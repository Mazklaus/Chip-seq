#!/usr/bin/env bash
# Ce script récupère les % de reads alignés, le nombre de reads récupérés et le nombre de gènes transcrits.
# Il se base sur le parcours du fichier conditionSheet et procède donc sample par sample.
# paramètre : il correspond au type de fichiers fastq (pairedEnd ou singleEnd)
# $1 correspond au fait que le fichier soit single ou pairedEnd
# $2 correspond au chemin des résultats

# DEFINITION DES FONCTIONS

# Calcul du Q30 depuis l'extrait du fastQC
Q30_calculation () {
  # Retire les 2 premières et les 2 dernières lignes de Q30.txt afin de ne conserver que les lignes indiquant la qualité de lecture des bases par position
  cat Q30.txt | head -n -1 | tail -n +3 > Q30bis.txt
  # Calcul du Q30
  NBBASES=$(sed -n '$=' Q30bis.txt)
  NBQ30=0
  while read laBase leQ30 leReste
  do
    if [ $(echo "$leQ30 > 30" | bc) -eq 1 ];
    then
      let NBQ30=NBQ30+1
    fi
  done < Q30bis.txt
  let NBQ30=100*NBQ30/NBBASES
  # Retourne le Q30
  echo "Q30 ecrit: ${NBQ30}" | cat >> $LOG
  echo ${NBQ30}
}


# Extraction des données depuis les fichiers bam et counts
data_extraction () {
  # Retire l'extension .bam afin d'obtenir le répertoire tophat
  BAMREP=$(echo ${THEFILE} | sed "s/\(.*\)\/alignment.bam/\1/")
  if [ ! -z ${THEFILE} ]
  then
    # Extraction du pourcentage de reads correctement alignés
    MAPPEDREADS=$(cat ${BAMREP}/align_summary.txt | tail -1 | sed "s/^\ *\(.*\)/\1/" | cut -d' ' -f1) # le sed permet de supprimer les éventuels espaces en début de la dernière ligne
    MAPPEDREADS=$(echo $MAPPEDREADS | sed "s/\%//g")  # Retire les %
    # Extraction du nombre de reads analysés
    NBREADS=$(cat ${BAMREP}/align_summary.txt | head -2 | tail -1)
    NBREADS=$(echo ${NBREADS} | sed "s/\ //g" | sed "s/Input:\(.*\)/\1/") # le 1er sed supprime le 1er espace en début de la dernière ligne
    if [[ "$FASTQ_TYPE" =~ pairedEnd ]]
    then
      if [ ! -z ${FASTQC_PATH}/${THEFASTQC1}.zip ] && [ ! -z ${FASTQC_PATH}/${THEFASTQC2}.zip ]
      then
        # Extraction du Q30 depuis les fastQC
        unzip -o ${FASTQC_PATH}/${THEFASTQC1}.zip -d ${FASTQC_PATH}
        unzip -o ${FASTQC_PATH}/${THEFASTQC2}.zip -d ${FASTQC_PATH}
        # Récupération de la partie des fichiers fastqc_data.txt permettant le calcul du Q30 puis calcul des Q30
        sed -n '/^>>Per\ base\ sequence\ quality/,/^>>END_MODULE/w Q30.txt' ${FASTQC_PATH}/${THEFASTQC1}/fastqc_data.txt
        Q301=$(Q30_calculation)
        sed -n '/^>>Per\ base\ sequence\ quality/,/^>>END_MODULE/w Q30.txt' ${FASTQC_PATH}/${THEFASTQC2}/fastqc_data.txt
        Q302=$(Q30_calculation)
        let Q30=(Q301+Q302)/2
      else
        echo "    FastQC ${THEFASTQC1} or ${THEFASTQC2} missing" | cat >> $LOG
      fi
      # Extraction du pourcentage de reads avec alignement multiple
      MULTIREADS=$(cat ${BAMREP}/align_summary.txt | tail -3 | head -1)
      MULTIREADS=$(echo ${MULTIREADS} | sed 's/[^\(]*.\([^%]*\).*/\1/')
      # Extraction du pourcentage de reads discordant
      DISCORDANTREADS=$(cat ${BAMREP}/align_summary.txt | tail -2 | head -1)
      DISCORDANTREADS=$(echo ${DISCORDANTREADS} | sed 's/[^\(]*.\([^%]*\).*/\1/')
    else
      if [ ! -z ${FASTQC_PATH}/${THEFASTQC}.zip ]
      then
        # Extraction du Q30 depuis les fastQC
        unzip -o ${FASTQC_PATH}/${THEFASTQC}.zip -d ${FASTQC_PATH}
        # Récupération de la partie des fichiers fastqc_data.txt permettant le calcul du Q30 puis calcul des Q30
        sed -n '/^>>Per\ base\ sequence\ quality/,/^>>END_MODULE/w Q30.txt' ${FASTQC_PATH}/${THEFASTQC}/fastqc_data.txt
        Q30=$(Q30_calculation)
      else
        echo "    FastQC ${THEFASTQC} missing" | cat >> $LOG
      fi
      # Extraction du pourcentage de reads avec alignement multiple
      MULTIREADS=$(cat ${BAMREP}/align_summary.txt | tail -2 | head -1)
      MULTIREADS=$(echo ${MULTIREADS} | sed 's/[^\(]*.\([^%]*\).*/\1/')
      # Extraction du pourcentage de reads discordant --> pas possible en singleEnd
      DISCORDANTREADS=""
    fi
    # Calcul du nombre de gènes exprimés
    COUNTSFILE=${COUNTS_PATH}/${THECOUNTS}
    if [ ! -z ${COUNTSFILE} ]
    then
      head -n -5 ${COUNTSFILE} |  (while read count1 count2
      do
        if [ "$count2" -gt 0 ]
        then
          let NBGENES=NBGENES+1
        fi
      done
      echo "${col1},${col2},${leReste},${MAPPEDREADS},${MULTIREADS},${DISCORDANTREADS},${NBREADS},${NBGENES},${Q30}" >> $OUTPUT)
    else
      echo "    Count ${COUNTSFILE} missing" | cat >> $LOG
    fi
  else
    echo "    Bam ${col2} missing" | cat >> $LOG
  fi
}


if [ "$#" -ne 2 ];
then

  echo "Argument manquant pour l'exécution du script QC.sh"

else

  # INITIALISATION DES VARIABLES
  FASTQ_TYPE=$1
  BAM_PATH=$2/aligmnent_res
  COUNTS_PATH=$2/aligmnent_res
  FASTQC_PATH=$2/quality_control
  LOG=$2/QC/log_QC.txt #fichier de sortie contenant ce qui c'est passée
  OUTPUT=$2/QC/QC.csv #fichier contenant les résultats
  ENTETE=$(head -1 $2/conditionSheet.csv)
  ENTETE="${ENTETE},Mapped_reads,Multiple_mapped,Discordant_reads,Reads_number,Expressed_genes,Q30"

  # Initialisation des fichiers QC.csv et log_QC.txt
  echo $ENTETE > $OUTPUT
  echo -n "" > $LOG


  # PROGRAMME PRINCIPAL
  echo "QC started"

  while IFS=',' read col1 col2 leReste
  do
    # Réinitialisation des variables
    MAPPEDREADS=""
    MULTIREADS=""
    DISCORDANTREADS=""
    NBREADS=""
    NBGENES=0

    THESAMPLE="alignment.bam"
    THESAMPLE=$(echo $THESAMPLE | sed "s/\"//g")  # Retire les guillemets
    THECOUNTS="${col2}.counts"
    THECOUNTS=$(echo $THECOUNTS | sed "s/\"//g")  # Retire les guillemets
    if [[ "$FASTQ_TYPE" =~ pairedEnd ]]
    then
      THEFASTQC1="${col2}.R1_fastqc"
      THEFASTQC1=$(echo $THEFASTQC1 | sed "s/\"//g")  # Retire les guillemets
      THEFASTQC2="${col2}.R2_fastqc"
      THEFASTQC2=$(echo $THEFASTQC2 | sed "s/\"//g")  # Retire les guillemets
    else
      THEFASTQC="${col2}_fastqc"
      THEFASTQC=$(echo $THEFASTQC | sed "s/\"//g")  # Retire les guillemets
    fi
    THEFILE=$(find ${BAM_PATH}/${col2} -maxdepth 1 -type f -name "${THESAMPLE}")
    echo "QC du sample ${col2}" | cat >> $LOG
    data_extraction
    echo "--"

  done < $2/conditionSheet.csv

  rm ./Q30.txt
  rm ./Q30bis.txt
  echo "QC ended"

fi
