#!/bin/bash

# ============== global variables defined here ========= # start
declare RXTHREAD_HOME           #-> root directory
if [ -z "${RXTHREAD_HOME}" ]
then
	#echo "RXTHREAD_HOME not set. Use default value '~/GitBucket'"
	RXTHREAD_HOME=~/GitBucket
fi


# ----- usage ------ #
usage()
{
	echo "USAGE: $0 <-i list1> <-r root1> <-I list2> <-R root2> <-p program> "
	echo "               [-o out_root] [-j job_id] [-z TMSorGDT] [-n NAMorNOT] [-s options] [-f suffix] "
	echo "               [-K remove_tmp] [-H home] " 
	echo "[note1]:  <listX> and <rootX> are relative or absolute directory. "
	echo "          <jobid> would be the name of output and temporary directory. "
	echo "          <program> should be either DeepScore or DeepAlign. "
	echo "          DeepScore, DeepAlign should exist. "
	echo "[note2]:  Default value of 'job_id' is struct, default value of 'out_root' is './\${job_id}_out "
	echo "          'TMSorGDT' (default = 0) for output TMscore[1], or uGDT[2], or both[0]. "
	echo "          'NAMorNOT' (default = 0) for output header[1], or not[0]. "
	echo "          'options' (default = null) is the arguments for <program>, such as '-n -1' "
	echo "          'suffix'  (default = null) is for the suffix name to align. "
	echo "[note3]:  default home directory of DeepAlign_package. "
	echo "          [default = $RXTHREAD_HOME/DeepAlign_package]"
	echo "[note4]:  remove_tmp is set to 1 by default to remove temporary folder. "
	exit 1
}

if [ $# -lt 1 ];
then
	usage
fi
curdir="$(pwd)"


# ----- get arguments ----- #
#-> optional arguments
out_root=""
jobid="struct"
TMSorGDT=0         #-> '0' for output both TMscore and uGDT; '1' for output TMscore; '2' for output uGDT
NAMorNOT=0         #-> '1' for output header; '0' for not output. 
options=""         #-> the arguments for <program>, such as '-n -1'
suffix=""          #-> suffix name 
#-> required arguments
list1=""
root1=""
list2=""
root2=""
program=""
kill_tmp=1      #-> default: kill temporary root
#-> home directory
home=$RXTHREAD_HOME/DeepAlign_package         #-> home directory


#-> parse arguments
while getopts ":i:r:I:R:p:j:z:n:s:f:K:H:" opt; 
do
	case $opt in
	#-> required arguments
	i)
		list1=$OPTARG
		;;
	r)
		root1=$OPTARG
		;;
	I)
		list2=$OPTARG
		;;
	R)
		root2=$OPTARG
		;;
	p)
		program=$OPTARG
		;;
	#-> optional arguments
	o)
		out_root=$OPTARG
		;;
	j)
		jobid=$OPTARG
		;;
	z)
		TMSorGDT=$OPTARG
		;;
	n)
		NAMorNOT=$OPTARG
		;;
	s)
		options=$OPTARG
		;;
	f)
		suffix=$OPTARG
		;;
	K)
		kill_tmp=$OPTARG
		;;
	#-> home directory
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




#=========================== initial check =================================#

# ------ check home directory ---------- #
if [ ! -d "$home" ]
then
	echo "home directory $home not exist " >&2
	exit 1
fi
home=`readlink -f $home`

# ------ check output directory -------- #
if [ "$out_root" == "" ]
then
	out_root=$curdir/${jobid}_out
fi
mkdir -p $out_root
out_root=`readlink -f $out_root`

# ------ check required arguments ------ #
if [ ! -f "$list1" ]
then
	echo "list1 $list1 not found !!" >&2
	exit 1
fi
list1=`readlink -f $list1`

if [ ! -f "$list2" ]
then
	echo "list2 $list2 not found !!" >&2
	exit 1
fi
list2=`readlink -f $list2`

if [ ! -d "$root1" ]
then
	echo "root1 $root1 not found !!" >&2
	exit 1
fi
root1=`readlink -f $root1`

if [ ! -d "$root2" ]
then
	echo "root2 $root2 not found !!" >&2
	exit 1
fi
root2=`readlink -f $root2`

if [ "$program" != "DeepScore" ] && [ "$program" != "DeepAlign" ]
then
	echo "program must be either DeepScore or DeepAlign !!" >&2
	exit 1
fi


# ========================== main process ============================#

# --- create temporary folder --#
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
tmp_root="${out_root}/TMP_STRUCT_${jobid}_${RANDOM}_${DATE}"
mkdir -p $tmp_root

#------ extract lists with the first column ----#
awk '{print $1}' $list1 > $tmp_root/$jobid.list1
awk '{print $1}' $list2 > $tmp_root/$jobid.list2
relist1=$tmp_root/$jobid.list1
relist2=$tmp_root/$jobid.list2

#----- mpi run matrix ------#
proclist=$tmp_root/$jobid.strucsco_matrix_proc
tmpdir=$tmp_root/tmp"_"$jobid"_"strucsco_matrix
mkdir -p $tmpdir
for f in `cat $relist1`
do
	for g in `cat $relist2`
	do
		$home/$program $root1/$f$suffix $root2/$g$suffix -P 0 $options > $tmpdir/$f-$g.tmp_sco &
	done
done
wait

#----- get TMscore ------#
if [ $TMSorGDT -eq 0 ] || [ $TMSorGDT -eq 1 ]
then
	TmscoOut=$out_root/$jobid.tmsco_out
	rm -f $TmscoOut
	#-> print header
	if [ $NAMorNOT -eq 1 ]
	then
		str="-----"
		for f in `cat $relist2`
		do
			str=$str" "${f:0:5}
		done
		echo $str >> $TmscoOut
	fi
	#-> print content
	for f in `cat $relist1`
	do
		if [ $NAMorNOT -eq 1 ]
		then
			str=${f:0:5}
		else
			str=""
		fi
		for g in `cat $relist2`
		do
			score=`awk '{print $12}' $tmpdir/$f-$g.tmp_sco`
			str=$str" "$score
		done
		echo $str >> $TmscoOut
	done
fi

#----- get uGDT ------#
if [ $TMSorGDT -eq 0 ] || [ $TMSorGDT -eq 2 ]
then
	uGDTout=$out_root/$jobid.ugdt_out
	rm -f $uGDTout
	#-> print header
	if [ $NAMorNOT -eq 1 ]
	then
		str="-----"
		for f in `cat $relist2`
		do
			str=$str" "${f:0:5}
		done
		echo $str >> $uGDTout
	fi
	#-> print content
	for f in `cat $relist1`
	do
		if [ $NAMorNOT -eq 1 ]
		then
			str=${f:0:5}
		else
			str=""
		fi
		for g in `cat $relist2`
		do
			score=`awk '{printf $NF}' $tmpdir/$f-$g.tmp_sco`
			sco=${score:0:5}
			str=$str" "$sco
		done
		echo $str >> $uGDTout
	done
fi

#---- delete temporary files -----#
if [ $kill_tmp -eq 1 ]
then
	rm -rf $tmpdir
else
	rm -rf $out_root/"TMP_STRUCT_"${jobid}
	mv $tmpdir $out_root/"TMP_STRUCT_"${jobid}
fi
rm -rf $tmp_root
rm -f $relist1
rm -f $relist2

#---- success exit ------#
exit 0

