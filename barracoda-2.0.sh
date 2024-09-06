#!/bin/bash

# 2022-01-19
# Maintained by Eva Kjær and Kamilla K. Munk
# Based on program by Andrea Marquard 2015-03-20

#####################################################################
# Define paths to R
#R="/home/local/tuba-nobackup/shared/R/R-3.2.0/bin/R"
R="/opt/R-4.3.1/bin/R"
# Define paths to bowtie tools
bowtie2="/home/local/tuba-nobackup/shared/bin/bowtie2-align"
bowtie2Build="/home/local/tuba-nobackup/shared/bin/bowtie2-build"

# Define path to barracoda scripts
barracoda_script_dir="/home/local/barracoda/tools/barracoda-2.0/scripts"

# Default storage dir
default_storage_dir="/home/local/barracoda/archive"


#####################################################################
# Help message <3
help() {
	echo
	echo
	echo -e '  B A R R A C O D A'
	echo 
	echo "  Original version by Andrea Marquard, March 2015"
	echo 
	echo "  Description: Analyse FASTQ file from sequencing of DNA barcodes"
	echo "  Usage: $0 -f <FASTQ-file> -t <tag-dir> -m <sample-map-file>"
	echo
	echo
	echo -e '  Required arguments:'
	echo "  Input files:"
	echo "  -f   Sequencing data file (fasta or fastq). Example: seq_data.fastq"
	echo "  -m   Sample identification table - a table mapping sample ID tags to sample names. Example: sample_id_table.xlsx"
	echo "  -a   Barcode annotations table- a table with barcode annotations. Example: barcode_annotations_small.xlsx"
	echo
	echo "  Barcode information:"
	echo "  -A   Sample identification tag - a FASTA file with sample tags. Example: sample_id_tags.fasta"
	echo "  -B   Forward primer A sequence. Example: GAAGTTCCAGCCAGCGTCACAGTTT"
	echo "  -C   N-sequence length (A end). Example: 6"
	echo "  -D   Epitope tag A FASTA file.  Example: a.fasta"
	echo "  -E   Annealing region sequence. Example: GGTCAGCATCATTTCC"
	echo "  -F   Epitope tag B FASTA file.  Excample: b.fasta"
	echo "  -G   N-sequence length (B end). Esample: 6"
	echo "  -H   Forward primer B sequence. Example: CAATCTTGAGCGTGACTTAAG"
	echo
	echo
	echo -e '  Optional arguments:'
	echo "  -p   Barcode plate setup table - table with barcode plate layouts. Example: barcode_plate_setup.xlsx"
	echo "  -s   Storage directory in which jobid folder exist (Default:" $default_storage_dir")."
	echo "       By the end of a Barracoda run the jobid folder will contain:"
	echo "					- Input directory to which input files are moved and stored."
	echo "					- Intermediate directory in which intermediate file is stored"
	echo " 					- Output directory in which user (web) output is stored"
	echo "  -o   Optional path that output directory is copied to (e.g. path to web server)"
        echo "  -c   Sum counts for duplicated sampes. Can be TRUE or FALSE (Default:FALSE)" # NEW     
	echo "  -k   Keep all intermediate files (Default:Off)"           
        echo "  -w   Webserver mode (Default:Off)"  
	echo "  -h   Print this help information."
	echo 
	echo
}

#####################################################################
# Define variables based on arguments
[ $# -eq 0 ] && help
# Add ":" after characters that needs arguments (not k and h)
sum_of_counts="FALSE" # NEW Used when dealing with duplicated sampes

while getopts "f:m:a:A:B:C:D:E:F:G:H:p:s:o:c:wkh" arg; do # NEW
  case $arg in
	  # REQUIRED ARGUMENTS:
	  f) seq_data_fastq=${OPTARG} ;; # Could be .fasta
	  m) sample_id_table_file=${OPTARG} ;; # Could be .txt
		a) barcode_annotations_xlsx=${OPTARG} ;; 

		A) sample_id_tags_fasta=${OPTARG} ;; 
		B) a_forward_primer_seq=${OPTARG} ;; 
		C) a_end_n_seq_length=${OPTARG} ;; 
		D) a_epitope_tag_fasta=${OPTARG} ;; 
		E) annealing_seq=${OPTARG} ;; 
		F) b_epitope_tag_fasta=${OPTARG} ;; 
		G) b_end_n_seq_length=${OPTARG} ;; 
		H) b_forward_primer_seq=${OPTARG} ;; 

    # OPTIONAL ARGUMENTS:
    p) barcode_plate_setup_xlsx=${OPTARG} ;; 
		s) storage_dir=${OPTARG} ;;
		o) extra_output_dir=${OPTARG} ;; 
                c) sum_of_counts=${OPTARG} ;; # NEW
		w) web_mode=1 ;;
		k) keep_all=1 ;;
    h) help; exit 0 ;;
  esac
done

#####################################################################
#####################################################################
## Very initial checks

# Check if paths to R exists
if [ ! -f "$R" ]; then echo "$R does not exists."; fi

# Check if paths to bowtie tools exist
if [ ! -f "$bowtie2" ]; then echo "$bowtie2 does not exists."; fi
if [ ! -f "$bowtie2Build" ]; then echo "$bowtie2Build does not exists."; fi

# Check if path to barracoda scripts exists
if [ ! -d "$barracoda_script_dir" ]; then echo "Barracoda script directory $barracoda_script_dir does not exists."; fi

# Check if default storage dir exists
if [ ! -d "$default_storage_dir" ]; then echo "Default storage directory $default_storage_dir does not exists."; fi


#####################################################################
# Job id folder
if [ -d $default_storage_dir/store.barracoda_$(date +"%Y-%m-%d").1 ]; then last_num=$(printf "%s\n" $default_storage_dir/store.barracoda_$(date +"%Y-%m-%d").* | sort -Vr | head -1 | rev  | cut -d. -f1 | rev); else last_num=0; fi
jobid_dir="$default_storage_dir/store.barracoda_$(date +"%Y-%m-%d").$((last_num+1))"
mkdir $jobid_dir

# Redefine jobid_dir if -s was used
if [ -z $storage_dir ] ; then : ; else jobid_dir=$storage_dir/store.barracoda_${jobid_dir##*_} ; fi

# Autogenerated directories
input_dir="$jobid_dir/input"
output_dir="$jobid_dir/output" 
intermediate_dir="$jobid_dir/files_intermediate"
log_dir=$intermediate_dir/logs
mkdir -p $input_dir $output_dir $intermediate_dir $log_dir

# Make suer that we can still write to stdout. This is done by duping stdout to FD3 so we keep connection to CGI
exec 3>&1

# Save stdout and stderr to log file in output_dir
LogFile=`echo $output_dir/log_file_barracoda_${jobid_dir##*_}.txt`
exec >> $LogFile


#####################################################################
#####################################################################

# To the terminal
echo "The intermediate files and results can be found on TUBA: ${jobid_dir}" >&3

# Welcome message (to log file)
echo -e $'\n\n                     \n B A R R A C O D A \n                     '
echo $'\n You are now running Barracoda 2.0....\n In case of problems, please contact kamkj@dtu.dk, and provide the process id of your job:' barracoda_${jobid_dir##*_} $'\n The job folder is:' $jobid_dir

# Header function
Header() { 
	echo -e $'\n\n' $1 '\n' $(date +'Timestamp: %D at %T') '\n'
}

# show warnings and errors function 
ShowErrors() {
   warnmsg=`grep -h WARNING $log_dir/*.log $LogFile` # -h hide filename path
   errormsg=`grep -h ERROR $log_dir/*.log $LogFile`
   msg="${warnmsg}${errormsg}"
   if [[ $msg != "" ]] ; then
      echo "<h2>Warnings and Errors</h2>" >&3
      echo "<font color=\"red\">" >&3
      echo $"${msg}" >&3
      echo "</font>" >&3
   fi
}

ExitIfErrors() {
    errormsg=`grep -h ERROR $1` # -h hide filename path
    if [[ $errormsg != "" ]] ; then
       echo "<h2>Warnings and Errors</h2>" >&3
       echo "<font color=\"red\">" >&3
       echo $"${errormsg}" >&3
       echo "</font>" >&3
       exit
   fi
}

#####################################################################
#####################################################################
# Required arguments
Header 'Required arguments parsed to Barracoda'

RequiredOpt() {
	var=$(echo $1 | sed 's/var=//g')
	if [ -z $var ] ; then
		echo $' [ERROR] Script terminated because option' $2 $'is missing. Please provide the' $3 $'for this to work.\n\nCheck out the help page by using the -h option.\n\n';
		ShowErrors ;
		exit
	else
		echo $' The' $2 $'option was succesfully parsed to barracoda. The' $3 $'is defined as:\n       ' $var
	fi
}

RequiredOpt "var=${seq_data_fastq}" "-f" "Sequencing data file (fasta or fastq)"
RequiredOpt "var=${sample_id_table_file}" "-m" "sample identification table"
RequiredOpt "var=${barcode_annotations_xlsx}" "-a" "barcode annotations table"
RequiredOpt "var=${sample_id_tags_fasta}" "-A" " sample identification tag"
RequiredOpt "var=${a_forward_primer_seq}" "-B" "forward primer A sequence"
RequiredOpt "var=${a_end_n_seq_length}" "-C" "N-sequence length (A end)"
RequiredOpt "var=${a_epitope_tag_fasta}" "-D" "epitope tag A FASTA file"
RequiredOpt "var=${annealing_seq}" "-E" "annealing region sequence"
RequiredOpt "var=${b_epitope_tag_fasta}" "-F" "epitope tag B FASTA file"
RequiredOpt "var=${b_end_n_seq_length}" "-G" "N-sequence length (B end)"
RequiredOpt "var=${b_forward_primer_seq}" "-H" "forward primer B sequence"


#####################################################################
#####################################################################
# Optional arguments
Header 'Optional arguments parsed to Barracoda'

OptionalOpt() {
	var=$(echo $1 | sed 's/var=//g')
	if [ -z $var ] ; then
		echo $' Barracoda will be run without' $3 $'because option' $2 $'was not used.\n  -> Check out the help page by using the -h option..'
	else
		echo $' The' $2 $'option was succesfully parsed to barracoda and the' $3 $'is defined:\n  -> Path:' $var$'.'
	fi
}
OptionalDir() {
	var=$(echo $1 | sed 's/var=//g')
	if [ -z $var ] ; then
		echo $'\n Default setting will be used to create the' $3 $'because option' $2 $'was not used.\n  -> Check out the help page by using the -h option..'
	else
		echo $'\n The' $2 $'option was succesfully parsed to barracoda and the' $3 $'is defined:\n  -> Path:' $var$'.'
	fi
} 

OptionalOpt "var=${barcode_plate_setup_xlsx}" "-p" "barcode plate setup table"
OptionalOpt "var=${extra_output_dir}" "-o" "extra output directory"
OptionalDir "var=${storage_dir}" "-p" "working directory"

echo $'\n The following directories were created:'
echo $'  -> Input directory:' $input_dir
echo $'  -> Intermediate files directory:' $intermediate_dir
echo $'  -> Output directory:' $output_dir


#####################################################################
#####################################################################
# Define paths to R scripts
Header 'Paths to R scripts and logs'

# Paths to scripts and logs
echo -e $'  -> Defining paths to R scripts'
function_script=$barracoda_script_dir/functions.R

checkinputdata_script=$barracoda_script_dir/r_checkinput.R
readlength_script=$barracoda_script_dir/plot-read-lengths.R
summarize_script=$barracoda_script_dir/summarize-barcodes.R
pvalscript=$barracoda_script_dir/collect-pvals-and-logfc.R
platesetup_script=$barracoda_script_dir/plot-barcodes-on-plates.R

# Path to log
r_logcheckinputdata=$log_dir/r_checkinput.R.log
r_log_plot=$log_dir/plot-read-lengths.R.log
r_log_sum=$log_dir/summarize-barcodes.R.log
r_log_pval=$log_dir/collect-pvals-and-logfc.R.log
r_log_plate=$log_dir/plot-barcodes-on-plates.R.log
echo -e $'  -> Log can be found in' $r_log


#####################################################################
#####################################################################
# Save files to input dir - not seq_data_fastq

file_type=$(file -b ${sample_id_table_file}) # find the file type of the sample_id_table_file, xlsx (Microsoft/Excel) or txt (text)
if [[ $file_type = *"Microsoft"* ]] || [[ $file_type = *"Excel"* ]] ; then cp $sample_id_table_file $input_dir/sample_id_table.xlsx ; sample_id_table_file=$input_dir/sample_id_table.xlsx ; fi
if [[ $file_type = *"ASCII"* ]] || [[ $file_type = *"text"* ]] ; then cp $sample_id_table_file $input_dir/sample_id_table.txt ; sample_id_table_file=$input_dir/sample_id_table.txt ; fi

cp $barcode_annotations_xlsx $input_dir/barcode_annotations.xlsx
barcode_annotations_xlsx=$input_dir/barcode_annotations.xlsx
cp $sample_id_tags_fasta $input_dir/sample_id_tags.fasta
sample_id_tags_fasta=$input_dir/sample_id_tags.fasta
cp $a_epitope_tag_fasta $input_dir/a_epitope_tag.fasta
a_epitope_tag_fasta=$input_dir/a_epitope_tag.fasta
cp $b_epitope_tag_fasta $input_dir/b_epitope_tag.fasta
b_epitope_tag_fasta=$input_dir/b_epitope_tag.fasta
if [ $barcode_plate_setup_xlsx ] ; then cp $barcode_plate_setup_xlsx $input_dir/barcode_plate_setup.xlsx ; fi

# Primer info
echo -e $'Forward primer A:' $a_forward_primer_seq $'(length:' $a_forward_primer_length $'nucleotides)\nN-sequence (A end):' $a_end_n_seq_length $'\nAnnealing region:' $annealing_seq $'(length:' $annealing_length $'nucleotides)\nN-sequence (B end):' $b_end_n_seq_length $'\nForward primer A:' $b_forward_primer_seq $'(length:' $b_forward_primer_length $'nucleotides)' > $input_dir/primer_info.txt

#####################################################################
# Primer seq to fasta
a_forward_primer_fasta=$input_dir/a_forward_primer.fasta
echo $'>a_forward_primer\n'$a_forward_primer_seq > $a_forward_primer_fasta
annealing_fasta=$input_dir/annealing.fasta
echo -e $'>annealing_region\n'$annealing_seq > $annealing_fasta
b_forward_primer_fasta=$input_dir/b_forward_primer.fasta
echo -e $'>b_forward_primer\n'$b_forward_primer_seq > $b_forward_primer_fasta


#####################################################################
#####################################################################
# Input data check
Header 'Checking input data'

#####################################################################
# Check that sequences only contain base pair letters (A, C, T, G)
echo -e $' Checking if sequences only contain base pair letters (A, C, T, G) ...'
for seq in $a_forward_primer_seq"_forward primer A" $annealing_seq"_annealing region" $b_forward_primer_seq"_forward primer B" ; do
	if [[ ${seq%_*} =~ ^[ACTG]+$ ]] ; then : ; else 
		echo -e $'  -> [ERROR] The' ${seq#*_} $'sequence contains other characters than A, C, T and G!!!!! The program was terminated!!\n Please provide sequences (options -B, -E and -H) only containing base pair letters (A, C, T, G)...\n\n' ; 
		ShowErrors ; 
		exit ; fi
done
echo -e $'  -> The forward A primer, annealing region and forward B primer sequences seem to be right..'

#####################################################################
# Make checks in R
echo -e $' Making checks in R '
echo -e $'  -> Running read checkinputdata_script with following command line:'
$R --vanilla --slave --args $sample_id_table_file $input_dir/sample_id_table.txt $barcode_annotations_xlsx < <( cat $function_script $checkinputdata_script ) > $r_logcheckinputdata 2>&1
echo -e $'  \t' $R '--vanilla --slave --args' $sample_id_table_file $input_dir/sample_id_table.txt $barcode_annotations_xlsx  '< <(cat' $function_script $checkinputdata_script ') >' $r_logcheckinputdata '2>&1'

##### EXIT IF checks went wrong..
ExitIfErrors $r_logcheckinputdata

#####################################################################
# Check that sample ids match in sample_id_table_file and sample_id_tags_fasta 
echo -e $'\n Checking if a-keys in sample_id_table file exist in sample FASTA file ...'
for a_key in $(cut -f1 $input_dir/sample_id_table.txt) ; do 
	if grep -q $a_key $sample_id_tags_fasta ; then : ; else 
		echo '  -> ' [ERROR] $a_key NOT found!!!!! ; 
		echo -e $'\n The program was terminated!!\n Please provide a sample id table with a-keys of sample FASTA file ...\n\n' ; 
		ShowErrors ;
		exit ; fi
done
echo -e $'  -> All a-keys in sample_id_table matches headers in sample_id_tags_fasta..'

#####################################################################

# Lengths of sequences in fasta files 
echo -e $'\n Lengths and positions of sequences '
CheckLengthFasta() {
        if [[ $(echo $1 | wc -w) == 1 ]] ; then : ; else 
		echo -e $'  -> [ERROR] The' $2 $'FASTA file contains reads of different lengths.. The script was terminated because of a problem with the sequences (options -B, -E and -H)...\n\n' ; 
		ShowErrors ;
		exit ; fi
}

# exstract only sequence with awk, then remove all special characters from text, get unique length of the sequences
sample_id_tags_length=$(awk '{if(NR%2==0) print $1}' $sample_id_tags_fasta | sed $'s/[^[:print:]\t]//g' | awk '{print length}' | uniq)
CheckLengthFasta "$sample_id_tags_length" "sample identification tag"
a_epitope_tag_length=$(awk '{if(NR%2==0) print $1}' $a_epitope_tag_fasta | sed $'s/[^[:print:]\t]//g' | awk '{print length}' | uniq)
CheckLengthFasta "$a_epitope_tag_length" "epitope tag A"
b_epitope_tag_length=$(awk '{if(NR%2==0) print $1}' $b_epitope_tag_fasta | sed $'s/[^[:print:]\t]//g' | awk '{print length}' | uniq)
CheckLengthFasta "$b_epitope_tag_length" "epitope tag B"
echo -e $'  -> All fasta files contained reads of same lengths..'

# Length of primers
a_forward_primer_length=${#a_forward_primer_seq} 
annealing_length=${#annealing_seq} 
b_forward_primer_length=${#b_forward_primer_seq}

# Sum of lengths (expected barcode length)
expected_barcode_length=$((sample_id_tags_length+a_forward_primer_length+a_end_n_seq_length+a_epitope_tag_length+annealing_length+b_epitope_tag_length+b_end_n_seq_length+b_forward_primer_length))
echo -e $'  -> The expected barcode length is' $expected_barcode_length

# Make lengths array for later use
declare -A lengths
lengths[sample_id_tags_fasta]=$sample_id_tags_length
lengths[a_epitope_tag_fasta]=$a_epitope_tag_length
lengths[b_epitope_tag_fasta]=$b_epitope_tag_length
lengths[a_forward_primer_fasta]=$a_forward_primer_length
lengths[annealing_fasta]=$annealing_length
lengths[b_forward_primer_fasta]=$b_forward_primer_length

# EXIT 1  - check 

#####################################################################
# Check if sequencing data file is zipped
echo -e $'\n Checking sequence file ...'

fastx_type=$(file -b --mime-type ${seq_data_fastq}) # find the file type of the sequencing data file, will include zip if the file is zipped 
if [[ $fastx_type = *"zip"* ]] ; then
        echo "  -> Sequencing file is zipped, unzipping using gunzip..." ;
        gunzip $seq_data_fastq ; seq_data_fastq=${seq_data_fastq%.*} ;
        echo -e $'  -> Unzipping was succesful! Sequencing file is now:\n      ' $seq_data_fastq ;
else : ; fi

# Check if sequencing data file is fastq or fasta
fastx_format=$(head -n 1 ${seq_data_fastq} | cut -c1-1) # this will extract the character of the first line in the sequencing data file if this is a "@" the file is a fastq file and if this is a ">" the file is a fasta file
fastx_suffix='fastq' # can be either 'fastq' or 'fasta'. Changed to 'fasta' if the first line in the sequencing data file if this is a ">" 
if [[ $fastx_format == "@" ]] ; then echo -e $'  -> The sequencing file is FASTQ file, we can carry on :)' ;
        elif [[ $fastx_format == ">" ]] ; then fastx_suffix='fasta' ; echo -e $'  -> The sequencing file is a FASTA file, we can carry on :)' ;
        else echo -e $'  -> [ERROR] The sequencing file is neither FASTQ nor FASTA format!!! The program was terminated!\n\n Please input a sequencing file that is either FASTQ or FASTA!!!\n' ; 
		ShowErrors ; 
		exit ;
fi

# Check if sample id table file is xlsx or txt
if [ ${sample_id_table_file##*.} = "xlsx" ] ; then echo -e $'  -> The sample id table is XLSX, we can carry on :)' ;
	elif [ ${sample_id_table_file##*.} = "txt" ] ; then echo -e $'  -> The sample id table is TXT, we can carry on :)' ;
	else echo -e $'  -> [ERROR] The sample id table is neither XLSX nor TXT format!!! The program was terminated!\n\n Please input a sample id table i that is either XLSX or TXT!!!\n' ; 
		ShowErrors ;
		exit ;
fi

#####################################################################
# Settings based on file type

# Define 'divide_for_nr_reads' for calculating number of reads based on sequencing data file type
if [[ $fastx_format == "@" ]] ; then divide_for_nr_reads=4 ;  paste_for_sorting='- - - -' ; modulus_division_for_length=4 ; line_per_modulus_for_length=2 # if fastq
elif [[ $fastx_format == ">" ]] ; then divide_for_nr_reads=2 ; paste_for_sorting='- -' ; modulus_division_for_length=2 ; line_per_modulus_for_length=0 # if fasta
fi


#####################################################################
#####################################################################
# Positions of sequences
Header 'Determining positions of sequences and adding to fasta'
sample_id_tags_start=1
let sample_id_tags_end=$sample_id_tags_length
let a_forward_primer_start=$sample_id_tags_end+1
let a_forward_primer_end=$sample_id_tags_end+$a_forward_primer_length
let a_end_n_seq_start=$a_forward_primer_end+1
let a_end_n_seq_end=$a_forward_primer_end+$a_end_n_seq_length
let a_epitope_tag_start=$a_end_n_seq_end+1
let a_epitope_tag_end=$a_end_n_seq_end+$a_epitope_tag_length
let annealing_start=$a_epitope_tag_end+1
let annealing_end=$a_epitope_tag_end+$annealing_length
let b_epitope_tag_start=$annealing_end+1
let b_epitope_tag_end=$annealing_end+$b_epitope_tag_length
let b_end_n_seq_start=$b_epitope_tag_end+1
let b_end_n_seq_end=$b_epitope_tag_end+$b_end_n_seq_length
let b_forward_primer_start=$b_end_n_seq_end+1
let b_forward_primer_end=$b_end_n_seq_end+$b_forward_primer_length

# Add positons fasta files (0 indexed!)
fasta_dir_pos=$intermediate_dir/fasta_with_positions
mkdir $fasta_dir_pos
echo -e $' A new directory has been created to store fasta with positions:\n  ' $fasta_dir_pos $'\n\n Lengths and positions has been added to fasta files:'
AddPositionsToFasta() {
	fasta=$1 ; name=$2 ; length=$3 ; start=$4 ; end=$5 ; newfile=$6
	sed 's/>.*/&			'"${start} ${end}"'/' $fasta | sed 's/   / /g' >> $newfile
	echo -e $'  -> The' $name $'sequence length is' $length $'( Start:' $start 'End:' $end ")"
}
AddPositionsToFasta $sample_id_tags_fasta "sample identification tag" $sample_id_tags_length $sample_id_tags_start $sample_id_tags_end $fasta_dir_pos/sample_id_tags.fasta
AddPositionsToFasta	$a_forward_primer_fasta "forward primer A" $a_forward_primer_length $a_forward_primer_start $a_forward_primer_end $fasta_dir_pos/a_forward_primer.fasta
echo -e $'  -> The N-sequence (A end) length length is' $a_end_n_seq_length $'( Start:' $a_end_n_seq_start 'End:' $a_end_n_seq_end ")" 
AddPositionsToFasta $a_epitope_tag_fasta "epitope tag A" $a_epitope_tag_length $a_epitope_tag_start $a_epitope_tag_end $fasta_dir_pos/a_epitope_tag.fasta
AddPositionsToFasta $annealing_fasta "annealing region" $annealing_length $annealing_start $annealing_end $fasta_dir_pos/annealing.fasta
AddPositionsToFasta $b_epitope_tag_fasta "epitope tag B" $b_epitope_tag_length $b_epitope_tag_start $b_epitope_tag_end $fasta_dir_pos/b_epitope_tag.fasta
echo -e $'  -> The N-sequence (B end) length length is' $b_end_n_seq_length $'( Start:' $b_end_n_seq_start 'End:' $b_end_n_seq_end ")" 
AddPositionsToFasta	$b_forward_primer_fasta "forward primer B" $b_forward_primer_length $b_forward_primer_start $b_forward_primer_end $fasta_dir_pos/b_forward_primer.fasta


#####################################################################
#####################################################################
# Bowtie index databases
Header 'Bowtie index databases <3'
declare -A alignment_stats

# Making bowtie dir
echo -e $' A new directory has been created to store fasta with positions:\n  ' $fasta_dir_pos
bowtie_dir="$intermediate_dir/bowtie_alignment"
mkdir -p $bowtie_dir

#####################################################################
# Bowtie indices with Bowtie2build
echo -e $'\n Creating indices of sequence files with bowtie2Build  '
mkdir -p $bowtie_dir/databases
for file in a_forward_primer_fasta a_epitope_tag_fasta annealing_fasta b_epitope_tag_fasta b_forward_primer_fasta ; do
	echo -e $'  -> Creating index for' $file":" ${!file}
	mkdir -p $bowtie_dir/databases/$file
	python3 $bowtie2Build -f ${!file} $bowtie_dir/databases/$file/tag &> $bowtie_dir/databases/$file/build.log
	#python $bowtie2Build -f ${!file} $bowtie_dir/databases/$file/tag &> $bowtie_dir/databases/$file/build.log
done

#####################################################################
# Make seq_data dir
seq_data_dir=$intermediate_dir/seq_data
mkdir $seq_data_dir

# Sort sequencing data file by id
echo -e $'\n Sorting sequencing data file ... '
sorted_seq_data_fastq=$seq_data_dir/seq_data_sorted.${fastx_suffix}
cat $seq_data_fastq | paste $paste_for_sorting | sort -k1,1 -t ' ' | tr '\t' '\n' > $sorted_seq_data_fastq
echo -e $'  -> Sorted file saved as:' $sorted_seq_data_fastq

# Number of reads for alignment_stats array
alignment_stats[total]=$(($(wc -l < $sorted_seq_data_fastq)/$divide_for_nr_reads))
echo -e $'  -> Total number of reads:' ${alignment_stats[total]}

# Length of reads
read_lengths_dir=$intermediate_dir/read_lengths
mkdir $read_lengths_dir
read_lengths_all=$read_lengths_dir/read_lengths_all.txt
( echo -e "length\treads" ; awk -v modulus_division=$modulus_division_for_length -v line=$line_per_modulus_for_length 'NR%modulus_division==line {print length}' $sorted_seq_data_fastq | sort -nr | uniq -c | awk '{print $2"\t"$1}' ) > $read_lengths_all

echo -e $'  -> Length of reads saved as:' $read_lengths_all


#####################################################################
#####################################################################
# Bowtie2 alignment - all reads
Header 'Bowtie alignment - all reads <3'
	# --nofw    do not align forward (original) version of read 
	# --norc 		do not align reverse-complement version of read 

# Bowtie alignment - primers and annealing region (constant sequences)
echo -e $' Aligning reads to constant sequences (primers and annealing region) with bowtie2 '
echo -e $'COMMAND LINES TO RUN BOWTIE:\n' > $log_dir/bowtie.commands.log
for filesetting in a_forward_primer_fasta'.--norc' annealing_fasta'.--norc' b_forward_primer_fasta'.--nofw' ; do	
	file=$(echo $filesetting | cut -d'.' -f1)
	no_setting=$(echo $filesetting | cut -d'.' -f2)
	length=$(echo ${lengths[$file]})
	mkdir -p $bowtie_dir/$file
	echo -e $'  -> Aligning' $file 'to reads using the' $no_setting 'setting and length of' $length
	cat $sorted_seq_data_fastq | $bowtie2 -x $bowtie_dir/databases/$file/tag -t -p 10 -i L,1,0 -U - -S $bowtie_dir/$file/aligns.sam -N 1 -L 10 --un $bowtie_dir/$file/un.${fastx_suffix} --score-min C,$length --reorder $(echo $no_setting) --no-hd --local --al $bowtie_dir/$file/al.${fastx_suffix} 2> $bowtie_dir/$file/log
	echo -e 'cat' $sorted_seq_data_fastq '|' $bowtie2 '-x' $bowtie_dir'/databases/'$file/tag '-t -p 10 -i L,1,0 -U - -S' $bowtie_dir/$file'/aligns.sam -N 1 -L 10 --un' $bowtie_dir'/'$file'/un.'${fastx_suffix} '--score-min C,'$length '--reorder' $(echo $no_setting) '--no-hd --local --al' $bowtie_dir'/'$file'/'al'.'${fastx_suffix} $'\n' >> $log_dir/bowtie.commands.log
	alignment_stats[$file]=$(($(wc -l < $bowtie_dir/$file/al.${fastx_suffix})/$divide_for_nr_reads))
done

# Filtering fastq of alignment files to continue with reads that aligned to 2 out of 3
	# Figure 7.1 s. 57 i Andreas PhD afhandling
echo -e $'\n Filtering sequencing file to only keep reads that aligned to 2 of 3 constant sequences '
seq_data_2of3_fastq=$seq_data_dir/seq_data_2of3.${fastx_suffix}
echo -e $"  -> Filtering with 'make-fastq-from-sam.pl' script"
$barracoda_script_dir/make-fastq-from-sam.pl $bowtie_dir/a_forward_primer_fasta/aligns.sam $bowtie_dir/annealing_fasta/aligns.sam $bowtie_dir/b_forward_primer_fasta/aligns.sam 2 ${fastx_suffix} > $seq_data_2of3_fastq
echo -e $"  -> Calculating number of reads for alignment statistics"
alignment_stats[seq_data_2of3_fastq]=$(($(wc -l < $seq_data_2of3_fastq)/$divide_for_nr_reads))

# Output alignment statistics
echo -e $'\n Alignment statistics for constant sequences (number of reads): '
	echo $'  -> Total (all reads):\t'	${alignment_stats[total]}
	echo $'  -> Forward primer A:\t'	${alignment_stats[a_forward_primer_fasta]}
	echo $'  -> Annealing region:\t'	${alignment_stats[annealing_fasta]}
	echo $'  -> Forward primer B:\t'	${alignment_stats[b_forward_primer_fasta]}
	echo $'  -> Aligned to 2 of 3:\t'	${alignment_stats[seq_data_2of3_fastq]}

# Length of reads
read_lengths_2of3=$read_lengths_dir/read_lengths_2of3.txt
( echo -e "length\treads" ; awk -v modulus_division=$modulus_division_for_length -v line=$line_per_modulus_for_length 'NR%modulus_division==line {print length}' $seq_data_2of3_fastq | sort -nr | uniq -c | awk '{print $2"\t"$1}' ) > $read_lengths_2of3
echo -e $'  -> Length of reads saved as:' $read_lengths_2of3


#####################################################################
#####################################################################
# Bowtie2 alignment - 2 of 3
Header 'Bowtie alignment - 2 of 3 <3'
	# --nofw    do not align forward (original) version of read 
	# --norc 		do not align reverse-complement version of read 

# Bowtie aligment - epitope tags 
echo -e $' Aligning filtered reads to epitope tags with bowtie2 '
	# --nofw    do not align forward (original) version of read 
	# --norc 		do not align reverse-complement version of read 
for filesetting in a_epitope_tag_fasta'.--norc' b_epitope_tag_fasta'.--nofw' ; do	
	file=$(echo $filesetting | cut -d'.' -f1)
	no_setting=$(echo $filesetting | cut -d'.' -f2)
	length=$(echo ${lengths[$file]})
	mkdir -p $bowtie_dir/$file
	echo -e $'  -> Aligning' $file 'to filtered reads using the' $no_setting 'setting and length of' $length
	cat $seq_data_2of3_fastq | $bowtie2 -x $bowtie_dir/databases/$file/tag -t -p 10 -i L,1,0 -U - -S $bowtie_dir/$file/aligns.sam -N 1 -L 10 --score-min C,$length --reorder $(echo $no_setting) --no-unal --no-hd --local --al $bowtie_dir/$file/al.${fastx_suffix} 2> $bowtie_dir/$file/log
	alignment_stats[$file]=$(($(wc -l < $bowtie_dir/$file/al.${fastx_suffix})/$divide_for_nr_reads))
done

# Merge alignmets from epitope tags
echo -e $'\n Filter sequencing file to only keep reads that aligned to epitope tags A and B '
seq_data_AandB_sam=$seq_data_dir/seq_data_AandB_sam.sam
echo -e $"  -> Filtering with 'merge-sam-by-read-ID.pl' script"
$barracoda_script_dir/merge-sam-by-read-ID.pl $bowtie_dir/a_epitope_tag_fasta/aligns.sam $bowtie_dir/b_epitope_tag_fasta/aligns.sam > $seq_data_AandB_sam
echo -e $"  -> Calculating number of reads for alignment statistics"
alignment_stats[seq_data_AandB_sam]=$(($(wc -l < $seq_data_AandB_sam)))

# Output alignment statistics
echo -e $'\n Alignment statistics for epitope tags (number of reads): '
	echo $'  -> Total (Aligned to 2 of 3):\t'	${alignment_stats[seq_data_2of3_fastq]}
	echo $'  -> Epitope tag A:      \t'	${alignment_stats[a_epitope_tag_fasta]}
	echo $'  -> Epitope tag B:      \t'	${alignment_stats[b_epitope_tag_fasta]}
	echo $'  -> Aligned to A and B:\t'	${alignment_stats[seq_data_AandB_sam]}

# Length of reads in sam file 
read_lengths_AandB=$read_lengths_dir/read_lengths_AandB.txt
cut -f4 $seq_data_AandB_sam | $barracoda_script_dir/stdin-lengths.pl > $read_lengths_AandB
echo -e $'  -> Length of reads saved as:' $read_lengths_AandB


#####################################################################
#####################################################################
# MERGED MAPPED READS -> to locate the sample id and N6 sequence
Header 'IT IS TIME TO MAP SOME MERGED READS!'

# Map constant sequences (primers and annealing region) and epitope tags within each read
echo -e $' Mapping constant sequences (primers and annealing region) and epitope tags within each read '

# Chunks to run in parallel in order to save time
echo -e $'  -> Creating chunks to run in parallel!'
chunk_dir=$intermediate_dir/chunks
mkdir $chunk_dir
lines=$(cat $seq_data_AandB_sam | wc -l)
njobs=30
chunk_size=$((1 + (lines / njobs)))
split -l $chunk_size $seq_data_AandB_sam $chunk_dir/chunk

# Dissecting barcodes fast
echo -e $'  -> Dissecting barcodes fast'
ls -1 $chunk_dir | parallel --joblog ${chunk_dir}/log -j $njobs "${barracoda_script_dir}/dissect-barcodes--fast.pl ${chunk_dir}/{} ${fasta_dir_pos} > ${chunk_dir}/{}.mapped.reads.txt 2> ${chunk_dir}/{}.log" 2> $chunk_dir/log

# Map constant sequences (primers and annealing region) and epitope tags within each read
echo -e $'  -> Making merged.mapped.reads.txt file!'
merged_mapped_reads=$intermediate_dir/merged.mapped.reads.txt
awk 'NR==1 || FNR>1' $chunk_dir/*.mapped.reads.txt > $merged_mapped_reads

# -s == file exists and is not empty and has more than one headerline 
checkFile() {
        if [[ $(wc -l <$1) -ge 2 ]]; then
                echo $' The' $2 $' file was created succesfully.\n'
        else
                echo $' [ERROR] Script terminated because' $2 $'was not generated!!!\n';
                ShowErrors ;
                exit
        fi

}

checkFile "${merged_mapped_reads}" "merged.mapped.reads.txt"

#####################################################################
#####################################################################

# Remove created files to minimize footprint
Header 'Clean up to minimize footprint'

if [[ "$keep_all" = 1 ]] ; then 
	echo -e $'  -> Since the -k option was used, files will not be removed '
else
	echo -e $'  -> Removing seq_data_dir, bowtie_dir, and chunk_dir'
	rm -r $seq_data_dir
	rm -r $bowtie_dir
	rm -r $chunk_dir
fi


#####################################################################
#####################################################################
# Analysis in R incl. plotting read lengths and summarizing barcodes
Header 'Analysis in R'

# Carrying out analysis in R
echo -e $' Carrying out R analysis including summarization of barcodes '

# Run read lengths scripts
echo -e $'  -> Running read lengths scripts with following command line'
$R --vanilla --slave --args $read_lengths_all $read_lengths_2of3 $read_lengths_AandB $output_dir $expected_barcode_length < <(cat $function_script $readlength_script ) > $r_log_plot 2>&1
echo -e $'  \t' $R '--vanilla --slave --args' $read_lengths_all $read_lengths_2of3 $read_lengths_AandB $output_dir $expected_barcode_length '< <(cat' $function_script $readlength_script ') >' $r_log_plot '2>&1'

# Run summarize barcodes scripts
echo -e $'  -> Running summarize barcodes with following command line'
$R --vanilla --slave --args $merged_mapped_reads $sample_id_table_file $output_dir $a_epitope_tag_fasta $b_epitope_tag_fasta $a_end_n_seq_length $b_end_n_seq_length $barcode_annotations_xlsx $sum_of_counts < <(cat $function_script $summarize_script ) > $r_log_sum 2>&1
echo -e $'  \t' $R '--vanilla --slave --args' $merged_mapped_reads $sample_id_table_file $output_dir $a_epitope_tag_fasta $b_epitope_tag_fasta $a_end_n_seq_length $b_end_n_seq_length $barcode_annotations_xlsx '< <(cat' $function_script $summarize_script ') >' $r_log_sum '2>&1'

# Run pval results
echo -e $'  -> Running read lengths scripts with following command line'
$R --vanilla --slave --args $output_dir < <(cat $function_script $pvalscript ) > $r_log_pval 2>&1
echo -e $'  \t' $R '--vanilla --slave --args' $output_dir '< <(cat' $function_script $pvalscript ') >' $r_log_pval '2>&1'

# Run barcode plate 
echo -e $'\n Plot of barcodes on plates '
if [ -z $barcode_plate_setup_xlsx ] ; then
		echo -e $'  -> The plot was not made because the option -p was not used'
		echo -e $'  -> Check out the help page by using the -h option..'
	else
		echo -e $'  -> Plots are being made in R...'
		$R --vanilla --slave --args $output_dir $input_dir/barcode_plate_setup.xlsx $barcode_annotations_xlsx < <(cat $function_script $platesetup_script ) > $r_log_plate 2>&1
		echo -e $'  \t' $R '--vanilla --slave --args' $output_dir $barcode_plate_setup_xlsx $barcode_annotations_xlsx  '< <(cat' $function_script $platesetup_script ') >' $r_log_plate '2>&1'
fi

# make files accesiable for all  
chmod -R 777 $jobid_dir

#####################################################################
#####################################################################
# Last-minute copying
Header 'Last-minute copying'

# Copy log file to log directory
echo -e $' Log file should be kept in log_dir '
echo -e $'  -> Moving log file to' $log_dir
cp  $output_dir/log_file_barracoda_${jobid_dir##*_}.txt $log_dir

# Copy output to other location with -o argument
echo -e $'\n Output directory copied to other location'
if [ -z $extra_output_dir ] ; then 
		echo -e $'  -> The outout directory was not copied because the option -o was not used'
		echo -e $'  -> Check out the help page by using the -h option..'
else 
	echo -e $'  -> The outout directory copied to other location:'
	echo -e $'    ' $extra_output_dir/barracoda_${jobid_dir##*_}
        
        mkdir -p $extra_output_dir # make directory if it does not exists
        webresultdir="barracoda_${jobid_dir##*_}"
	cp -r $output_dir $extra_output_dir/$webresultdir # copy the output directory to new location defined by the option -o
        cd $extra_output_dir
        zip -rq ${webresultdir}.zip ${webresultdir}/*
        
fi

#####################################################################
#####################################################################
# Heartfelt message
Header 'Goodbye'
echo $' The script did not hault with error...... Well done. :)\n'


#####################################################################
#####################################################################
# Make html for the webserver output 

if [[ "$web_mode" = 1 ]] ; then

   # Create HTML output
   webdir=`echo $extra_output_dir | cut -d'/' -f5-`
   webdir="/tuba/${webdir}/"
   #webdir="$extra_output_dir/$webresultdir"  
   echo "In case of problems, please email kamkj (at) dtu.dk, and provide the process id of your job: ${jobid_dir##*_}" >&3

   # Show warnings with the ShowErrors function
   ShowErrors

   # Write a reference to the results
   echo "<h2>Download results</h2>" >&3
   echo "<a href=\"${webdir}${webresultdir}.zip\">Download results as .zip file</a>" >&3
   echo "<a href=\"${webdir}${webresultdir}\">Show entire results folder</a>" >&3
	
   echo "<h2>Summary of data and analysis</h2>" >&3
   
   # show read lengths plot
   echo "<h4>NGS read lengths</h4>" >&3
   echo "<img src='${webdir}${webresultdir}/read-lengths.png'>" >&3
   
   # show total reads per A-key plot	
   echo "<h4>Distribution of reads among sample keys</h4>" >&3
   echo "<img src='${webdir}${webresultdir}/total-reads-per-key.png'>" >&3

fi


# How to copy newest script from mac
# mount_tuba
# cp /Users/evakjr/Documents/DTU/Barracoda/2022/ny\ barracoda\ \<3/barracoda-2.0.sh /Users/evakjr/Documents/DTU/Tuba/tuba/barracoda/barracoda-2.0 

# How to run: 
	# 100k
		# ./barracoda-2.0.sh -f test_data/seq_data.fastq -m test_data/sample_id_table.txt -a test_data/barcode_annotations.xlsx -A test_data/sample_id_tags.fasta -B GAAGTTCCAGCCAGCGTCACAGTTT -C 6 -D test_data/a.fasta -E GGTCAGCATCATTTCC -F test_data/b.fasta -G 6 -H GTTATCGGCTCGTTCACACTCGA -p test_data/barcode_plate_setup.xlsx

	# 10K
		# ./barracoda-2.0.sh -f test_data/test_10k.fastq -m test_data/sample_id_table.txt -a test_data/barcode_annotations.xlsx -A test_data/sample_id_tags.fasta -B GAAGTTCCAGCCAGCGTCACAGTTT -C 6 -D test_data/a.fasta -E GGTCAGCATCATTTCC -F test_data/b.fasta -G 6 -H GTTATCGGCTCGTTCACACTCGA -p test_data/barcode_plate_setup.xlsx

	# 10K xlsx sample_id_table
		# /home/tuba/barracoda/barracoda-2.0/barracoda-2.0.sh -f test_data/test_10k.fastq -m test_data/sample_id_table.xlsx -a test_data/barcode_annotations.xlsx -A test_data/sample_id_tags.fasta -B GAAGTTCCAGCCAGCGTCACAGTTT -C 6 -D test_data/a.fasta -E GGTCAGCATCATTTCC -F test_data/b.fasta -G 6 -H GTTATCGGCTCGTTCACACTCGA -p test_data/barcode_plate_setup.xlsx

# stort eksempel 
# cd /home/people/kamkj/run_barracoda/Susana/Bseq188/
# /home/tuba/barracoda/barracoda-2.0/barracoda-2.0.sh -f data/Bseq188.fastq -m data/Sample_identification_table.xlsx  -a data/Annotations-removed-barcodes-with-no-inputs.xlsx -A data/2_SampleIentificationTag_71_to_259_nov_2021.fasta -B GAAGTTCCAGCCAGCGTCACAGTTT -C 6 -D data/3_EpitopeTagA_new_barcodes.fasta -E GGTCAGCATCATTTCC -F data/4_EpitopeTagB_plate_43-54.fasta -G 6 -H CAATCTTGAGCGTGACTTAAG



