#!/usr/bin/env bash

# This a aligmnent test pipeline in order to manage RNA-seq data
# This script need to be adapt for both single and pair-end
# It will be better if i put a conditional here allowing to do both
# relative or absolute path to $1 $2

#Definition of the options
#: means an option is required
TEMP=`getopt -o hqdo:p:s:c: --long help,quality,strand-specific,output,paired,snps,core -n 'pipe.sh' -- "$@"`
eval set -- "$TEMP"

#Initialise the flags that will be modified by the different options.
argIndex=""
argFastq=""
argPaired=""
argSnps=""
argOutput="./"
argCore=1
argQualityControl=0
argStrandSpecific=0

# function for help
usage() {
	#clear

	echo -e "#######################################################################\n"
	echo -e "\e[1mNAME\e[0m\n"
	echo -e "pipe.sh is a shell pipeline based on the hisat2 package which allow to easely align reads from\nRNA-seq experiment to a reference genome. It's final goal is to cover the entire analysis.\n"
	echo -e "\e[1mDESCRIPTION\e[0m\n"
	echo -e "\t\e[1m./pipe.sh [OPTION(S)] [INDEX FILE] [FASTQ FILE] [ANNOTATION FILE]\e[0m\n"
	echo -e "\t\e[1m[INDEX FILE]\e[0m\t\t\tThe fasta file that will be used to build the index.\n"
	echo -e "\t\e[1m[FASTQ FILE]\e[0m\t\t\tThe fastq file containing the reads you want to align.\n"
	echo -e "\t\e[1m[ANNOTATION FILE]\e[0m\t\tThe gtf file containing the position of genes for your\n\t\t\t\t\treferecence genome.\n"
	echo -e "\e[1mOPTION\e[0m\n"
	echo -e "\t\e[1m-q|--quality\e[0m\t\t\tGive a html quality summary using fastqc.\n"
	echo -e "\t\e[1m-d|--strand-specific\e[0m\t\tUse this argument if your RNA-seq experiment is\n\t\t\t\t\tstrand specific.\n"
	echo -e "\t\e[1m-o|--output [DIRECTORY]\e[0m\t\tThe output folder where you want to place result.\n\t\t\t\t\t\e[31m/!\ \e[0m If it does not exist it will be create.\n"
	echo -e "\t\e[1m-p|--paired [.fastq FILE]\e[0m\tUsed this argument with the paired file if you have pair-end\n\t\t\t\t\tdata (without it asume that you are managing single end data).\n"
	echo -e "\t\e[1m-s|--snps [.fa FILE]\e[0m\t\tGive a file adding snps information to the index.\n"
	echo -e "\t\e[1m-c/--core [INTEGER]\e[0m\t\tNumber of CPU core used for alignement.\n\t\t\t\t\t\e[31m/!\ \e[0m A high number can improve speed but can also reduce\n\t\t\t\t\toverall computer speed during the computing time.\n"
	echo -e "\t\e[1m-h|--help\e[0m\t\t\tPrint the help.\n"
	echo -e "\n#######################################################################"

	exit 1
}

# print help if no argument
if [[ "$@" = "--" ]]; then
  usage
	exit 1
fi


# Parsing input parameters
#According to the inputed options, flags are set to 1 (or the value of the option) to apply only the wanted transformations of data
while true ;
do
	case "$1" in
		-h|help)
			usage ; shift ;;
		-q|quality)
			argQualityControl=1 ; shift ;;
		-d|strand-specific)
			argStrandSpecific=1 ; shift ;;
		-o|output)
			case "$2" in
				"")
					shift 2 ;;
				*)
					argOutput=$2
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
            argSnps=$2
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
    				argCore=$2
    				shift 2;;
        esac ;;
		--) shift ; break ;;
		*) echo "The option $1 does not exist." ; exit 1 ;;
	esac
done


# Creation of the output folder
if [ ! "$argOutput" = "./" ]; then
	mkdir $argOutput
	argOutput="${argOutput}/"
fi
if [ $argQualityControl = 1 ];then
	mkdir ${argOutput}quality_control
fi
mkdir ${argOutput}built_index
mkdir ${argOutput}aligmnent_res
mkdir ${argOutput}QC

# Production of the index file
if [ ! -n "$argPaired" ]; then
	if [ -n "$argSnps" ]; then
	  hisat2-build -p $argCore --snp $argSnps $1 ${argOutput}built_index/refer
	else
	  hisat2-build -p $argCore $1 ${argOutput}built_index/refer
	fi
else
	if [ -n "$argSnps" ]; then
		hisat2-build -p $argCore --snp $argSnps -1 $1 -2 $argSnps ${argOutput}built_index/refer
	else
		hisat2-build -p $argCore -1 $1 -2 $argSnps ${argOutput}built_index/refer
	fi
fi

# Alignement of the file
hisat2 -p $argCore -x ${argOutput}built_index/refer -U $2 -S ${argOutput}aligmnent_res/aligmnent.sam

# turn sam file into bam one (binaries) which are much more faster in compute
samtools view -bS ${argOutput}aligmnent_res/aligmnent.sam > ${argOutput}aligmnent_res/alignment.bam

# Fastqc creation of output
if [ $argQualityControl = 1 ]; then
	fastqc -o ${argOutput}quality_control ${argOutput}aligmnent_res/alignment.bam
fi

# Generation of the count table
if [ $argStrandSpecific = 1 ]; then
	htseq-count -f bam -r pos -s yes ${argOutput}aligmnent_res/alignment.bam $3 > ${argOutput}aligmnent_res/samples.counts
else
	htseq-count -f bam -r pos -s no ${argOutput}aligmnent_res/alignment.bam $3 > ${argOutput}aligmnent_res/samples.counts
fi

# Count control quality QC
# in this part the zip file from fastqc will be taken and put in used
# with the count file from htsq_count
