#!/bin/bash

#$ -cwd   -l virtual_free=4G
#$ -o /dev/null
#$ -e /dev/null

date

OLD_IFS=IFS
export OLD_IFS

COUNTER=0
while read -d ^ LINE
do

   COUNTER=$(($COUNTER+1))
   if [ $COUNTER == $SGE_TASK_ID ]; then
       IFS="
"
       set -- $LINE
       array=($LINE)

       for COMMAND in ${array[@]}; do
           echo $COMMAND
           eval $COMMAND
       done
   fi
done < $INFILE

date

