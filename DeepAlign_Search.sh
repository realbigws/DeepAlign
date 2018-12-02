#!/bin/bash

# ============== global variables defined here ========= # start
declare RXTHREAD_HOME           #-> root directory
if [ -z "${RXTHREAD_HOME}" ]
then
	#echo "RXTHREAD_HOME not set. Use default value '~/GitBucket'"
	RXTHREAD_HOME=~/GitBucket
fi

# ============== global variables defined here ========= # start
declare MAXSIZE=3000            #-> maximal subjext length   (default: 3000)
declare FUNC_RET                #-> for the function echo value
# ============== global variables defined here ========= # end



#---------------------------------------------------------#
##### ===== All functions are defined here ====== #########
#---------------------------------------------------------#


# ----- usage ------ #
function usage()
{
	echo "DeepAlign_Search v1.00 [Nov-30-2018] "
	echo "    Search a template database to find structurally similar proteins for a query protein structure "
	echo ""
	echo "USAGE:  ./DeepAlign_Search.sh <-q query_pdb> [-l data_list] [-d data_db] [-L refer_list] [-D refer_db] "
	echo "                              [-t tmsco] [-p pval] [-k topK] [-s score_func] [-C cut_alignment] [-c CPU_num]"
	echo "                              [-o output_root] [-O output_file] [-s sort] [-S options] [-H home] "
	echo "Options:"
	echo ""
	echo "***** required arguments *****"
	echo "-q query_pdb         : Query protein file in PDB format. "
	echo ""
	echo "***** relevant directories ***"
	echo "-l data_list         : The list of protein database [default = databases/bc40_list]"
	echo ""
	echo "-d data_db           : The folder containing the database files (i.e., .pdb files) [default = databases/pdb_BC100/]"
	echo ""
	echo "-L refer_list        : The list of reference database [default = databases/reference_pdb_list]"
	echo ""
	echo "-D refer_db          : The folder containing the reference files (in .pdb format) [default = databases/CAL_PDB/]"
	echo ""
	echo "***** optional arguments *****"
	echo "-t tmsco             : Apply TMscore cutoff during searching process [default = 0.35]"
	echo ""
	echo "-p pval              : Keep the results for top proteins according to P-value cutoff [default = 0.001; set -1 to disable]"
	echo ""
	echo "-k topK              : Keep the results for top topK proteins [default = 100]"
	echo ""
	echo "-s score_func        : 1:distance-score,2:vector-score,4:evolution-score; these scores could be combined [default 7]"
	echo ""
	echo "-C cut_alignment     : If specified, then the final template structure will be cut according to the alignment [default 0]"
	echo ""
	echo "-c CPU_num           : Number of processors. [default = 1]"
	echo ""
	echo "***** output arguments *****"
	echo "-o output_root       : The root for output file. At least one brief summary 'query.SortedScore_pvalue' will be generated."
	echo "                       [default = './\${input_name}_DeepSearch'] "
	echo ""
	echo "-O output_file       : The file containing a complex summary of the topK prediction results."
	echo "                       if not specified, then only a brief summary of the results will be generated. [default = null]"
	echo ""
	echo "***** other arguments *****"
	echo "-s sort              : screen output is sorted with respect to a specific column. [default = 8 for DeepScore] "
	echo "   tnam qnam tlen qlen -> BLOSUM CLESUM DeepScore -> LALI RMSDval TMscore -> MAXSUB GDT_TS GDT_HA -> SeqID nLen dCut"
	echo "      1    2    3    4 (5)     6      7         8 (9)  10      11      12 (13)   14     15     16 (17)  18   19   20"
	echo ""
	echo "-S options           : the arguments for DeepAlign, such as '-n -1' [default = null]"
	echo ""
	echo "-H home              : home directory of DeepSearch (i.e., \$DeepSearchHome)."
	echo "                       [default = $RXTHREAD_HOME/DeepAlign_Package]"
	echo ""
	exit 1
}



#----- check file existence ------# -> this is a function
function file_exist()
{
	#-- input --#
	local file=${1}    #-> 1st input is the file content
	local name=${2}    #-> 2nd input is the file name
	#-- check --#
	if [ -z "$file" ]
	then
		echo "$name is null !!" >&2
		return 1
	fi
	if [ ! -s "$curdir/$file" ]
	then
		if [ ! -s "$file" ]
		then
			echo "$name $file not found !!" >&2
			return 1
		fi
	else
		file=$curdir/$file
	fi
	#-- return --#
	FUNC_RET=$file
	return 0
}


#----- check root existence ------# -> this is a function
function root_exist()
{
	FUNC_RET=""
	#-- input --#
	local root=${1}    #-> 1st input is the root content
	local name=${2}    #-> 2nd input is the root name
	#-- check --#
	if [ -z "$root" ]
	then
		echo "$name is null !!" >&2
		return 1
	fi
	if [ ! -d "$curdir/$root" ]
	then
		if [ ! -d "$root" ]
		then
			echo "$name $root not found !!" >&2
			return 1
		fi
	else
		root=$curdir/$root
	fi
	#-- return --#
	FUNC_RET=$root
	return 0
}




#-------------------------------------------------------------#
##### ===== get pwd and check DeepThreaderHome ====== #########
#-------------------------------------------------------------#


#------ current directory ------#
curdir="$(pwd)"


#-------- check usage -------#
if [ $# -lt 1 ];
then
	usage
fi


#---------------------------------------------------------#
##### ===== All arguments are defined here ====== #########
#---------------------------------------------------------#

# ----- get arguments ----- #
#-> required arguments
query_pdb=""

#-> optional arguments
#--| data related
data_list=""
data_db=""
refer_list=""
refer_db=""
#--| parameter related
tmsco=0.35                      #-> tmscore cutoff [0.35]
pval=0.001                      #-> pvalue cutoff [0.001]
topK=100                        #-> minimal output number [100]
score_func=7                    #-> score function [7: using all for score_func]
cut_alignment=0                 #-> cut alignment or not [0: we don't cut template]
CPU_num=1                       #-> default CPU_num is 1
#--| output related
output_root=""
output_file=""
#--| other arguments
sort_col=8                      #-> sort by DeepScore
#---- screen output format -----#
#   tnam qnam tlen qlen -> BLOSUM CLESUM DeepScore -> LALI RMSDval TMscore -> MAXSUB GDT_TS GDT_HA -> SeqID nLen dCut
#      1    2    3    4 (5)     6      7         8 (9)  10      11      12 (13)   14     15     16 (17)  18   19   20
options=""
#--| home directory
home=$RXTHREAD_HOME/DeepAlign_Package         #-> home directory


#-> parse arguments
while getopts ":q:l:d:L:D:t:p:k:s:C:c:o:O:a:s:S:H:" opt;
do
	case $opt in
	#-> required arguments
	q)
		query_pdb=$OPTARG
		;;
	#-> optional arguments
	#--| data related
	l)
		data_list=$OPTARG
		;;
	d)
		data_db=$OPTARG
		;;
	L)
		refer_list=$OPTARG
		;;
	D)
		refer_db=$OPTARG
		;;
	#--| parameter related
	t)
		tmsco=$OPTARG
		;;
	p)
		pval=$OPTARG
		;;
	k)
		topK=$OPTARG
		;;
	s)
		score_func=$OPTARG
		;;
	C)
		cut_alignment=$OPTARG
		;;
	c)
		CPU_num=$OPTARG
		;;
	#--| output related
	o)
		output_root=$OPTARG
		;;
	O)
		output_file=$OPTARG
		;;
	#--| other arguments
	s)
		sort_col=$OPTARG
		;;
	S)
		options=$OPTARG
		;;
	H)
		home=$OPTARG
		;;
	#-> default
	\?)
		echo "Invalid option: -$OPTARG" >&2
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		exit 1
		;;
	esac
done




#---------------------------------------------------------#
##### ===== Part 0: initial argument check ====== #########
#---------------------------------------------------------#


#========================== Part 0.1 check the relevant roots of input arguments =======================#

#----------- check home -----------------#
root_exist $home home
retval=$?
if [ "$retval" != "0"  ]
then
	exit 1
else
	home=$FUNC_RET
fi
LOCAL_HOME=${home}
#echo "LOCAL_HOME=$LOCAL_HOME"

#----------- check query  -----------#
if [ -z "$query_pdb" ]
then
	echo "input query_pdb is null !!" >&2
	exit 1
fi
fulnam=`basename $query_pdb`
relnam=${fulnam%.*}

#----------- check data  -----------#
if [ "$refer_list" == "" ]
then
	refer_list="${LOCAL_HOME}/databases/reference_pdb_list"
fi
if [ "$refer_db" == "" ]
then
	refer_db="${LOCAL_HOME}/databases/CAL_PDB/"
fi
if [ "$data_list" == "" ]
then
	data_list="${LOCAL_HOME}/databases/bc40_list"
fi
if [ "$data_db" == "" ]
then
	data_db="${LOCAL_HOME}/databases/pdb_BC100/"
fi

#-- get name ---#
refer_fulnam=`basename $refer_list`
refer_relnam=${refer_fulnam%.*}
data_fulnam=`basename $data_list`
data_relnam=${data_fulnam%.*}

#-- data length ---#
data_len=`wc $data_list | awk '{print $1}'`
if [ $topK -le 0 ] || [ $topK -ge $data_len ]
then
	topK=$data_len
fi

#========================== Part 0.2 check the relevant roots of output arguments ========================#

#----------- create output_root with absolute directory --------#
if [ "$output_root" == "" ]
then
	output_root=$curdir/${relnam}_DeepSearch
fi
dir_out_root=`dirname $output_root`
nam_out_root=`basename $output_root`
if [ "$dir_out_root" == "." ]
then
	if [ "$nam_out_root" == "." ]
	then
		output_root=$curdir
	else
		output_root=$curdir/$nam_out_root
	fi
fi
mkdir -p $output_root
#-> change to absolute
if [ ! -d "$curdir/$output_root" ]
then
	if [ ! -d "$output_root" ]
	then
		echo "outroot $output_root not found !!" >&2
		exit 1
	fi
else
	output_root=$curdir/$output_root
fi


#----------- assign output_file with absolute directory --------#
if [ "$output_file" != "" ]
then
	dir_output_file=`dirname $output_file`
	nam_output_file=`basename $output_file`
	if [ "$dir_output_file" == "." ]
	then
		output_file=$output_root/$output_file
	fi
else
	output_file="-1"
fi


#-----------------------------------------------------#
##### ===== Part 1: DeepSearch process ====== #########
#-----------------------------------------------------#

#-> create tmp root
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
tmp="tmp_${relnam}_${RANDOM}_${DATE}"
prot=${output_root}/${tmp}
mkdir -p $prot


#================== Part 1.1 estimate p-value from refer_list ==================# 

#---- default value ----#
MEAN=0
VARI=0.5
#---- calculate Z-score and E-value -------#
if [ 1 -eq "$(echo "$pval > 0" | bc)" ]    #-> calculate refer_list for p-value if pval > 0
then

	#--------- preliminary ----------#
	#--| cut refer_list into N threads
	${LOCAL_HOME}/util/List_Div_Shuf $refer_list $CPU_num $prot
	retval=$?
	if [ "$retval" != "0" ]
	then
		exit 1
	fi
	#--| run DeepAlign for these N-cut refer_list
	deepalign="DeepAlign $options -P 0 -u 0 -e $MAXSIZE "

	#--------- Z-score ----------#
	#--| calculate Z-score
	for ((i=0;i<$CPU_num;i++))
	do
		# screen output for 'proc=0'
		if [ $i -eq 0 ]   #-> for proc=0, we add option "-v"
		then
		  addi=" -v "
		else
		  addi=""
		fi
		# run program
		(${LOCAL_HOME}/${deepalign} ${addi} -f ${prot}/${refer_relnam}.$i -r $refer_db -q $query_pdb -m 3 -g 0 -h 1 -c 0 -p 0 -b ${prot}/${relnam}_${refer_relnam}_zscore.$i -w ${prot}/${relnam}_${refer_relnam}_tmpout.$i)&
		OUT=$?
		if [ $OUT -ne 0 ]
		then
			echo "${LOCAL_HOME}/${deepalign} ${addi} -f ${prot}/${refer_relnam}.$i -r $refer_db -q $query_pdb -m 3 -g 0 -h 1 -c 0 -p 0 -b ${prot}/${relnam}_${refer_relnam}_zscore.$i -w ${prot}/${relnam}_${refer_relnam}_tmpout.$i" >&2
			exit 1
		fi
	done
	wait
	#--| collect Z-score
	rm -f ${prot}/${relnam}_${refer_relnam}.Score_zsco
	for ((i=0;i<$CPU_num;i++))
	do
		rm -f ${prot}/${relnam}_${refer_relnam}_tmpout.$i
		cat ${prot}/${relnam}_${refer_relnam}_zscore.$i >> ${prot}/${relnam}_${refer_relnam}.Score_zsco
		rm -f ${prot}/${relnam}_${refer_relnam}_zscore.$i
	done
	#--| calculate mean/vari
	awk '{print $6}' ${prot}/${relnam}_${refer_relnam}.Score_zsco | sort -g -r | tail -n+4 > ${prot}/${relnam}_${refer_relnam}.rank_zscore_val
	reso=`${LOCAL_HOME}/util/Stat_List ${prot}/${relnam}_${refer_relnam}.rank_zscore_val`
	MEAN=`echo $reso | cut -d ' ' -f 3`
	VARI=`echo $reso | cut -d ' ' -f 7`
	rm -f ${prot}/${relnam}_${refer_relnam}.rank_zscore_val
	rm -f ${prot}/${relnam}_${refer_relnam}.Score_zsco

	#--------- E-value ----------#
	#--| calculate E-value
	for ((i=0;i<$CPU_num;i++))
	do
		# screen output for 'proc=0'
		if [ $i -eq 0 ]   #-> for proc=0, we add option "-v"
		then
		  addi=" -v "
		else
		  addi=""
		fi
		# run program
		(${LOCAL_HOME}/${deepalign} ${addi} -f ${prot}/${refer_relnam}.$i -r $refer_db -q $query_pdb -j 0 -m 0 -p 0 -w ${prot}/${relnam}_${refer_relnam}_evalue.$i -s $score_func)&
		OUT=$?
		if [ $OUT -ne 0 ]
		then
			echo "${LOCAL_HOME}/${deepalign} ${addi} -f ${prot}/${refer_relnam}.$i -r $refer_db -q $query_pdb -j 0 -m 0 -p 0 -w ${prot}/${relnam}_${refer_relnam}_evalue.$i -s $score_func" >&2
			exit 1
		fi
	done
	wait
	#--| collect E-value
	rm -f ${prot}/${relnam}_${refer_relnam}.Score_evd
	for ((i=0;i<$CPU_num;i++))
	do
		rm -f ${prot}/${refer_relnam}.$i
		cat ${prot}/${relnam}_${refer_relnam}_evalue.$i >> ${prot}/${relnam}_${refer_relnam}.Score_evd
		rm -f ${prot}/${relnam}_${refer_relnam}_evalue.$i
	done
	#--| calculate miu/beta
	awk '{print $a}' a=${sort_col} ${prot}/${relnam}_${refer_relnam}.Score_evd | sort -g -r | tail -n+4 > ${prot}/${relnam}_${refer_relnam}.rank_evalue_val
	${LOCAL_HOME}/util/Fitting_EVD ${prot}/${relnam}_${refer_relnam}.rank_evalue_val > ${output_root}/${relnam}_${refer_relnam}.pvalue_param
	rm -f ${prot}/${relnam}_${refer_relnam}.rank_evalue_val
	rm -f ${prot}/${relnam}_${refer_relnam}.Score_evd
fi


#================== Part 1.2 run main search for data_list ==================# 

if true 
then

	#--------- preliminary ----------#
	#--| cut refer_list into N threads
	${LOCAL_HOME}/util/List_Div_Shuf $data_list $CPU_num $prot
	retval=$?
	if [ "$retval" != "0" ]
	then
		exit 1
	fi
	#--| run DeepAlign for these N-cut refer_list
	deepalign="DeepAlign $options -P 0 -u 0 -e $MAXSIZE "

	#--------- main search ----------#
	#--| calculate main search
	for ((i=0;i<$CPU_num;i++))
	do
		# screen output for 'proc=0'
		if [ $i -eq 0 ]   #-> for proc=0, we add option "-v"
		then
		  addi=" -v "
		else
		  addi=""
		fi
		# run program
		(${LOCAL_HOME}/${deepalign} ${addi} -f ${prot}/${data_relnam}.$i -r $data_db -q $query_pdb -m 2 -g $MEAN -h $VARI -c $tmsco -p 0 -w ${prot}/${relnam}_${data_relnam}.$i -s $score_func)&
		OUT=$?
		if [ $OUT -ne 0 ]
		then
			echo "${LOCAL_HOME}/${deepalign} ${addi} -f ${prot}/${data_relnam}.$i -r $data_db -q $query_pdb -m 2 -g $MEAN -h $VARI -c $tmsco -p 0 -w ${prot}/${relnam}_${data_relnam}.$i -s $score_func" >&2
			exit 1
		fi
	done
	wait

	#--| collect main search
	rm -f ${prot}/${relnam}_${data_relnam}.Score
	for ((i=0;i<$CPU_num;i++))
	do
		rm -f ${prot}/${data_relnam}.$i
		cat ${prot}/${relnam}_${data_relnam}.$i >> ${prot}/${relnam}_${data_relnam}.Score
		rm -f ${prot}/${relnam}_${data_relnam}.$i
	done
	sort -g -r -k ${sort_col} ${prot}/${relnam}_${data_relnam}.Score > ${output_root}/${relnam}_${data_relnam}.SortedScore
	rm -f ${prot}/${relnam}_${data_relnam}.Score
fi


#================== Part 1.3 re-align topK ==================# 

#-- determine topK to output the detailed rank ----#
rel_topk=$topK
if [ 1 -eq "$(echo "$pval > 0" | bc)" ]  #--| consider topK by p-value 
then
	pval_para=${output_root}/${relnam}_${refer_relnam}.pvalue_param
	sorted_score=${output_root}/${relnam}_${data_relnam}.SortedScore
	sorted_file=${prot}/${relnam}_${data_relnam}.SortedScore_${sort_col}
	awk '{print $a}' a=${sort_col} $sorted_score > $sorted_file
	topk_from_pval=`${LOCAL_HOME}/util/TopK_by_EVD $sorted_file $pval_para $pval`
	if [ $topk_from_pval -gt $topK ]
	then
		rel_topk=$topk_from_pval
	fi
	rm -f $sorted_file
fi
#-> topK_list
topk_list=${output_root}/${relnam}_${data_relnam}_topKlist
head -n $rel_topk ${output_root}/${relnam}_${data_relnam}.SortedScore | awk '{print $1}' > $topk_list
topk_align=${output_root}/${relnam}_${data_relnam}_topKalign
mkdir -p $topk_align


#------------- re-align TopK -------------#
if true 
then

	#--------- preliminary ----------#
	#--| cut refer_list into N threads
	${LOCAL_HOME}/util/List_Div_Shuf $topk_list $CPU_num $prot
	retval=$?
	if [ "$retval" != "0" ]
	then
		exit 1
	fi
	#--| run DeepAlign for these N-cut refer_list
	deepalign="DeepAlign $options -P 0 -u 0 -e $MAXSIZE "

	#--------- main search ----------#
	#--| calculate main search
	for ((i=0;i<$CPU_num;i++))
	do
		# screen output for 'proc=0'
		if [ $i -eq 0 ]   #-> for proc=0, we add option "-v"
		then
		  addi=" -v "
		else
		  addi=""
		fi
		# run program
		(${LOCAL_HOME}/${deepalign} ${addi} -f ${prot}/${relnam}_${data_relnam}_topKlist.$i -r $data_db -q $query_pdb -m 0 -p 1 -d $topk_align -w ${prot}/${relnam}_${data_relnam}_topKsco.$i -s $score_func)&
		OUT=$?
		if [ $OUT -ne 0 ]
		then
			echo "${LOCAL_HOME}/${deepalign} ${addi} -f ${prot}/${relnam}_${data_relnam}_topKlist.$i -r $data_db -q $query_pdb -m 0 -p 1 -d $topk_align -w ${prot}/${relnam}_${data_relnam}_topKsco.$i -s $score_func" >&2
			exit 1
		fi
	done
	wait

	#--| cut_template or not
	if [ $cut_alignment -eq 1 ]
	then
		# cut templates
		for i in `cat $topk_list`
		do
			${LOCAL_HOME}/util/Domain_Proc $topk_align/${i}-${relnam}.fasta ${prot}/${i}-${relnam}.fasta_cut 0 0 0 1
			${LOCAL_HOME}/util/PDB_File_Cut $data_db/$i.pdb ${prot}/${i}-${relnam}.fasta_cut_seq1 ${prot}/$i.pdb 0
		done
		# realign templates
		for ((i=0;i<$CPU_num;i++))
		do
			## screen output for 'proc=0'
			if [ $i -eq 0 ]   #-> for proc=0, we add option "-v"
			then
			  addi=" -v "
			else
			  addi=""
			fi
			## run program
			(${LOCAL_HOME}/${deepalign} ${addi} -f ${prot}/${relnam}_${data_relnam}_topKlist.$i -r $prot -q $query_pdb -m 0 -p 1 -i 0 -d $topk_align -w ${prot}/${relnam}_${data_relnam}_topKsco.$i -s $score_func)&
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "${LOCAL_HOME}/${deepalign} ${addi} -f ${prot}/${relnam}_${data_relnam}_topKlist.$i -r $prot -q $query_pdb -m 0 -p 1 -i 0 -d $topk_align -w ${prot}/${relnam}_${data_relnam}_topKsco.$i -s $score_func" >&2
				exit 1
			fi
		done
		wait
		# clear templates
		for i in `cat $topk_list`
		do
			rm -f ${prot}/${i}-${relnam}.fasta_cut*
			rm -f ${prot}/$i.pdb
		done
	fi

	#--| collect main search
	rm -f ${prot}/${relnam}_${data_relnam}.TopKScore
	for ((i=0;i<$CPU_num;i++))
	do
		rm -f ${prot}/${relnam}_${data_relnam}_topKlist.$i
		cat ${prot}/${relnam}_${data_relnam}_topKsco.$i >> ${prot}/${relnam}_${data_relnam}.TopKScore
		rm -f ${prot}/${relnam}_${data_relnam}_topKsco.$i
	done
	sort -g -r -k ${sort_col} ${prot}/${relnam}_${data_relnam}.TopKScore > ${output_root}/${relnam}_${data_relnam}.SortedTopKScore
	rm -f ${prot}/${relnam}_${data_relnam}.TopKScore
fi


#================== Part 1.4 generate the detailed summary file ==================#

#-- output the result file ---#
if [ "$output_file" != "-1" ]     #-> need to output the detailed summary file
then
	rank_simp=${output_file}_simp
	rank_file=${output_root}/${relnam}_${data_relnam}.SortedTopKScore
	rank_root=${output_root}/${relnam}_${data_relnam}_topKalign
	pval_para=${output_root}/${relnam}_${refer_relnam}.pvalue_param
	run_command="$0 $@"
	${LOCAL_HOME}/util/DeepSearch_Rank $relnam $rank_file $rank_root $pval_para $output_file "${run_command}" $pval $data_len $rel_topk > $rank_simp
	retval=$?
	if [ "$retval" != "0" ]
	then
		exit 1
	fi
fi


#--- remove tmp -----#
rm -rf $prot


#======================= exit ====================#
exit 0


