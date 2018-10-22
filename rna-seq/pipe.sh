#!/usr/bin/env bash

# This a aligmnent test pipeline in order to manage RNA-seq data
# This script need to be adapt for both single and pair-end
# it will need to let the people choose
# add an option to add snp file for index building
# the option for p in hisat2 correspond to the number of
# core to use (more = faster)
# It will be better if i put a conditional here allowing to do both
# relative or absolute path to $1 $2
# add an option to allow users to get a fastqc output giving a quality control view

#Definition of the options
#: means an option is required
TEMP=`getopt -o hp:s:c: --long help,paired,snps,core -n 'pipe.sh' -- "$@"`
eval set -- "$TEMP"

#Initialise the flags that will be modified by the different options.
argIndex=""
argFastq=""
argPaired=""
argSnps=""
argCore=1

# function for help
usage() {
	#clear

	echo -e "#######################################################################\n"
	echo -e "\e[1mNAME\e[0m\n"
	echo -e "pipe.sh is a shell pipeline based on the hisat2 package which allow to easely align reads from\nRNA-seq experiment to a reference genome. It's final goal is to cover the entire analysis.\n"
	echo -e "\e[1mDESCRIPTION\e[0m\n"
	echo -e "\t\e[1m./pipe.sh [OPTION(S)] [INDEX FILE] [FASTQ file]\e[0m\n"
	echo -e "\t\e[1m[INDEX FILE]\e[0m\t\t\tThe fasta file that will be used to build the index.\n"
	echo -e "\t\e[1m[FASTQ file]\e[0m\t\t\tThe fastq file containing the reads you want to align.\n"
	echo -e "\e[1mOPTION\e[0m\n"
	echo -e "\t\e[1m-h|--help\e[0m\t\t\tPrint the help.\n"
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
else
  hisat2-build ./$argIndex ./built_index/refer
fi


# Alignement of the file
hisat2 -p $argCore -x ./built_index/refer -U ./$argFastq -S ./aligmnent_res/aligmnent.sam

# turn sam file into bam one (binaries) which are much more faster in compute
samtools view -bS aligmnent.sam > alignement.bam
