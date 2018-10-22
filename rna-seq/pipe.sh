#!/usr/bin/env bash

# This a aligmnent test pipeline in order to manage RNA-seq data
# This script need to be adapt for both single and pair-end
# it will need to let the people choose
# add an option to add snp file for index building
# the option for p in hisat2 correspond to the number of
# core to use (more = faster)
# It will be better if i put a conditional here allowing to do both
# relative or absolute path to $1 $2

#Definition of the options
#: means an option is required
TEMP=`getopt -o AGNpis:kr:v:S:h --long AT,GC,dupliPos,dupliID,subsetSNP:,keepIntermediate,randomIndividuals:,splitIndiv:,randomSNP:,help -n 'qualityControl.sh' -- "$@"`
TEMP=`getopt -o hi:f:p:s:c: --long help,index,fastq,paired,snps,core -n 'pipe.sh' -- "$@"`
eval set -- "$TEMP"

#Initialise the flags that will be modified by the different options.
argIndex=""
argFastq=""
argPaired=""
argSnps=""
argCore=1

# function for help
usage() {
	clear

	echo -e "#######################################################################\n"
	echo -e "pipe.sh is a shell script taking as input a fileset of genomic data, namely SNP genotypes in the standard Plink format (BED/BIM/FAM). Adding different options allow subsetting of your initial dataset, however relying on the Plink software for most of the operations. Operations can be inputted simultaneously.\n"
	echo -e "\e[1m./pipe.sh [OPTION(S)] [INDEX FILE]\e[0m\n"
	echo -e "\t\e[1m-A/--AT\e[0m\t\t\t\t\tRemoves any SNP for which the alleles are A/T or T/A.\n\t\t\t\t\t\tThese alleles are troublesome because there may be ambiguities in the strands.\n"
	echo -e "\t\e[1m-G/--GC\e[0m\t\t\t\t\tRemoves any SNP for which the alleles are G/C or C/G for the same reasons as above.\n"
	echo -e "\t\e[1m-N\e[0m\t\t\t\t\tRemoves any SNP for which one of the two alleles is N.\n\t\t\t\t\t\tThe joker nucleotide 'N' can be found when quality of typing is low and most of the alleles\n\t\t\t\t\t\tin individuals were not clearly defined.\n"
	echo -e "\t\e[1m-p/--dupliPos\e[0m\t\t\t\tRemoves every SNP sharing the same position.\n\t\t\t\t\t\t\e[31m/!\ \e[0m Do not add this option if you allow multiallelic SNP.\n"
	echo -e "\t\e[1m-i/--dupliID\e[0m\t\t\t\tRemoves every SNP sharing their rsID with one or more SNP.\n"
	echo -e "\t\e[1m-s/--subsetSNP [.bim FILE]\e[0m\t\tSubset the SNP of your fileset by the list of SNP of another fileset.\n"
	echo -e "\t\e[1m-r/--randomIndividuals [INTEGER]\e[0m\tRandomly pick a defined number of individuals and generates the corresponding subset.\n"
	echo -e "\t\e[1m-v/--splitIndividuals\e[0m\t\t\tRandomly picks a defined number of individuals and generates the corresponding dataset.\n\t\t\t\t\t\tAdditionally provides a dataset with the remaining individuals of the original dataset.\n"
	echo -e "\t\e[1m-S/--randomSNP\e[0m\t\t\t\tRandomly picks a defined number of SNP.\n"
	echo -e "\t\e[1m-k/--keepIntermediate\e[0m\t\t\tSince operations can be simultaneous (but are always executed in the same order),\n\t\t\t\t\t\tyou can choose to keep the datasets at each step.\n\t\t\t\t\t\t\e[31m/!\ \e[0m This can flood your folders with numerous and potentially heavy files.\n"
	echo -e "\n#######################################################################"

	exit 1
}


# Parsing input parameters

#According to the inputted options, flags are set to 1 (or the value of the option) to apply only the wanted transformations of data
while true ;
do
	case "$1" in
		-h|help)
			usage ; shift ;;
		-i|index)
			case "$2" in
				"")
					echo "You need to give the fasta file used for index the genome of the studied species" ;
          exit 1;;
				*)
					if [ -e $2 ]; then
						argIndex=$2
					else
						echo "The fasta file for the index does not exist !";
						exit 1;
					fi;
					shift 2 ;;
			esac ;;
		-f|fastq)
			case "$2" in
				"")
          echo "You need to give the fastaq file you want to align" ;
          exit 1;;
				*)
          if [ -e $2 ]; then
            argFastq=$2
          else
            echo "The fastaq file you are trying to align does not exist !";
            exit 1;
          fi;
          shift 2 ;;
			esac ;;
		-p|paired)
			case "$2" in
				"")
		      shift 2 ;;
				*)
          if [ -e $2 ]; then
            argPaired=$2
          else
            echo "Your paired file does not exist !";
            exit 1;
          fi;
          shift 2 ;;
			esac ;;
		-s|snps)
			case "$2" in
				"")
				  shift 2 ;;
				*)
          if [ -e $2 ]; then
            argPaired=$2
          else
            echo "The snps file you want to use does not exist !";
            exit 1;
          fi;
          shift 2 ;;
			esac ;;
      -c|core)
        case "$2" in
          "")
            shift 2 ;;
          *)
            re='^[0-9]+$'
    				if ! [[ $2 =~ $re ]] ; then
    					echo "error: $2 is not a number" >&2; exit 1
    				fi;
    				argRandom=$2
    				shift 2;;
        esac ;;
		--) shift ; break ;;
		*) echo "The option $1 does not exist." ; exit 1 ;;
	esac
done


# Creation of the output folder
mkdir ./built_index
mkdir ./aligmnent_res

# Production of the index file
if [ -n "$argPaired" ]; then
  hisat2-build --snp ./$argSnps ./$argIndex ./built_index/refer
elif
  hisat2-build ./$argIndex ./built_index/refer
fi


# Alignement of the file
hisat2 -p $argCore -x ./built_index/refer -U ./$argFastq -S ./aligmnent_res/aligmnent.sam

samtools view -bS aligmnent.sam > alignement.bam
