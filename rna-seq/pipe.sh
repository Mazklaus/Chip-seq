#!/usr/bin/env bash

# This a aligmnent test pipeline in order to manage RNA-seq data
# This script need to be adapt for both single and pair-end
# it will need to let the people choose
# add an option to add snp file for index building
# It will be better if i put a conditional here allowing to do both
# relative or absolute path to $1 $2
# need to add the possibility to choose the output directory

#Definition of the options
#: means an option is required
TEMP=`getopt -o hqp:s:c: --long help,quality,paired,snps,core -n 'pipe.sh' -- "$@"`
eval set -- "$TEMP"

#Initialise the flags that will be modified by the different options.
argIndex=""
argFastq=""
argPaired=""
argSnps=""
argCore=1
argQualityControl=0

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
	echo -e "\t\e[1m[ANNOTATION FILE]\e[0m\t\tThe gtf file containing the position of genes for your\n\t\t\t\t\treferecence genome(WORK IN PROGRESS not really needed).\n"
	echo -e "\e[1mOPTION\e[0m\n"
	echo -e "\t\e[1m-h|--help\e[0m\t\t\tPrint the help.\n"
	echo -e "\t\e[1m-q|--quality\e[0m\t\t\tGive a html quality summary using fastqc.\n"
	echo -e "\t\e[1m-p|--paired [.fastq FILE]\e[0m\tUsed this argument with the paired file if you have pair-end\n\t\t\t\t\tdata (without it asume that you are managing single end data).\n"
	echo -e "\t\e[1m-s|--snps [.fa FILE]\e[0m\t\tGive a file adding snps information to the index.\n"
	echo -e "\t\e[1m-c/--core [INTEGER]\e[0m\t\tNumber of CPU core used for alignement.\n\t\t\t\t\t\e[31m/!\ \e[0m A high number can improve speed but can also reduce\n\t\t\t\t\toverall computer speed during the computing time.\n"
	echo -e "\n#######################################################################"

	exit 1
}

# print help if no argument
if [[ "$@" = "--" ]]; then
  usage
	exit 1
fi


# Parsing input parameters

#According to the inputted options, flags are set to 1 (or the value of the option) to apply only the wanted transformations of data
while true ;
do
	case "$1" in
		-h|help)
			usage ; shift ;;
		-q|quality)
			argQualityControl=1 ; shift ;;
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
mkdir ./built_index
mkdir ./aligmnent_res
if [ $argQualityControl = 1 ];then
	mkdir ./quality_control
fi

# Production of the index file
if [ -n "$argPaired" ]; then
	if [ -n "$argSnps" ]; then
	  hisat2-build -p $argCore --snp ./$argSnps ./$1 ./built_index/refer
	else
	  hisat2-build -p $argCore ./$1 ./built_index/refer
	fi
else
	if [ -n "$argSnps" ]; then
		hisat2-build -p $argCore --snp ./$argSnps -1 ./$1 -2 $argSnps ./built_index/refer
	else
		hisat2-build -p $argCore -1 ./$1 -2 $argSnps ./built_index/refer
	fi
fi
# Alignement of the file
hisat2 -p $argCore -x ./built_index/refer -U ./$2 -S ./aligmnent_res/aligmnent.sam

# turn sam file into bam one (binaries) which are much more faster in compute
samtools view -bS ./aligmnent_res/aligmnent.sam > ./aligmnent_res/alignment.bam

# Fastqc creation of output
if [ $argQualityControl = 1 ]; then
	fastqc -o ./quality_control ./aligmnent_res/alignment.bam
fi

# Generation of the count table
# the -s parameters must be put on yes if the alignement is strand specific
htseq-count -f bam -r pos -s no ./aligmnent_res/alignment.bam $3 > ./aligmnent_res/samples.counts
