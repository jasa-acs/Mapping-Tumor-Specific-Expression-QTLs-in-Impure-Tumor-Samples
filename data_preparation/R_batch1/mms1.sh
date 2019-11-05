#!/bin/bash


###############################################################################
# This script will make multiple jobs to be submitted based on the variables  #
# that the user "plugs in". Once  you have put in your variables and executed #
# this script make sure that you verify that the output is correct and then   #
# run submit.job.$OUTFILE and this will submit the jobs to the cluster        #
###############################################################################

###############################################################################
# If you make changes to this file you can only change it within a unix shell #
# and not wordpad, word or any other windows text editor!!! If you do make    #
# changes with *ANY* windows text editor, this script will not work!!!        #
###############################################################################



###--- The variable "OUTFILE" will create files named what you put in the variable "OUTFILE" ---###

OUTFILE="step11_ligate"




###--- OUTFILEEXT is the extension that you want to add to the end of your file. e.g. if you want the end of your file to be named .R then put in .R, if you want your file to be named .sas then put in .sas ---###

OUTFILEEXT=".R"




###--- MYCODE is the code that you want to run against your variables ---###

MYCODE="step11_ligate.R"
CODE=`cat $MYCODE`



###--- LIB is a list of libraries or that you may want to load before you setup your variables. This is a list of libraries that you want to insert in yor code before you define you variables. Each library that you specify you must add a retur character after the library is defined. If you do not use this then make sure you set USELIB="0". If you do use LIB then USELIB must be set to 1 ---###

LIB=""


###--- This will tell the script to use LIB or not ---###
USELIB="0"


###--- These are the names of the variables ---###

VARNAME1="i"
VARNAME2=""
VARNAME3=""
VARNAME4=""
VARNAME5=""
VARNAME6=""
VARNAME7=""

###--- This is the number of variables you intend to use... If you only have 5 variables then put in 5. If you do not put in this field then it may break the output file ---###

NUMVAR="1"



###--- These are the variables that you will put in to run against your code, each variable is seperated by a space, do not use tabs or any other delimiting character ---###

VAR1="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
VAR2=""
VAR3=""
VAR4=""
VAR5=""
VAR6=""
VAR7=""



### This is the command that you want to use to submit your jobs to the cluster. If you want to submit the job to a different queue or different application, you will have to modify this variable"

BSUB="bsub -q week R CMD BATCH"



## DO NOT MODIFY ANY CODE BELOW THIS LINE OTHERWISE YOU MAY BREAK THIS SCRIPT ##




###--- Initialize the arrays to generate output ---###
array1=( $VAR1 );
array2=( $VAR2 );
array3=( $VAR3 );
array4=( $VAR4 );
array5=( $VAR5 );
array6=( $VAR6 );
array7=( $VAR7 );

###--- Set the araay length so we know how many times to loop ---###
arraylength=${#array1[@]}

i=0;



while [ "$i" -lt "$arraylength" ]
do
	
	fromarray1=${array1["$i"]};
	fromarray2=${array2["$i"]};
	fromarray3=${array3["$i"]};
	fromarray4=${array4["$i"]};
	fromarray5=${array5["$i"]};
	fromarray6=${array6["$i"]};
	fromarray7=${array7["$i"]};


	if [ $USELIB -eq "1" ]; then
	echo "$LIB" >> $OUTFILE.$i$OUTFILEEXT
	fi

	echo "$VARNAME1 = $fromarray1" >>	$OUTFILE.$i$OUTFILEEXT

	if [ $NUMVAR -ge "2" ]; then
	echo "$VARNAME2 = $fromarray2" >>	$OUTFILE.$i$OUTFILEEXT
	fi

	if [ $NUMVAR -ge "3" ]; then
	echo "$VARNAME3 = $fromarray3" >>	$OUTFILE.$i$OUTFILEEXT
	fi

	if [ $NUMVAR -ge "4" ]; then
	echo "$VARNAME4 = $fromarray4" >>	$OUTFILE.$i$OUTFILEEXT
	fi

	if [ $NUMVAR -ge "5" ]; then
	echo "$VARNAME5 = $fromarray5" >>	$OUTFILE.$i$OUTFILEEXT
	fi

	if [ $NUMVAR -ge "6" ]; then
	echo "$VARNAME6 = $fromarray6" >>	$OUTFILE.$i$OUTFILEEXT
	fi

	if [ $NUMVAR -ge "7" ]; then
	echo "$VARNAME7 = $fromarray7" >>	$OUTFILE.$i$OUTFILEEXT
	fi

	echo "$BSUB $OUTFILE.$i$OUTFILEEXT" >> submit-jobs.sh	

	echo "$CODE" >> $OUTFILE.$i$OUTFILEEXT

echo "Creating file $OUTFILE.$i$OUTFILEEXT"

let "i = $i +1";
done

chmod 700 submit-jobs.sh

