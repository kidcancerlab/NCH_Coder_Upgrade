#!/bin/bash
# Author : Ali Snedden
# Date   : 5/31/19
# License: MIT
# Purpose: 
#
# Notes :
#
# Questions:
#
# References :  
#
## Get section of genome of interest. It contains a fragment of the trascript 
## spaced with introns
READ_LEN=50
START=$1
EXON1_BEG=$2
EXON1_END=$3
EXON2_BEG=$4
EXON2_END=$5
INSERT=\
`awk -v START="$START" -v FINISH=$EXON2_END\
    'BEGIN{
        LEN=60;
        READLEN=FINISH-START;
        LINE=int(START/60+2);
        CHAR=int(START%60);
        #print LINE1_BEG " "CHAR1_BEG " " REMAIN
    }
    {
        if(NR==LINE){
            printf substr($0,CHAR,60);
            REMAIN = READLEN-(LEN-CHAR);
            while(REMAIN > 0){
                getline;
                if(REMAIN >= 60)
                    printf substr($0,1,60) 
                if(REMAIN < 60)
                    printf substr($0,1,REMAIN) 
                REMAIN = REMAIN - 60
            }
        }
    }\
' ~/Code/Python/simulate_fastq_data/validate/data/chr1_short.fa`
#echo $INSERT

STR1=`python3 -c "import sys; insert=sys.argv[1]; print(insert[0:int(sys.argv[3]) - int(sys.argv[2])+1])" $INSERT $START $EXON1_END ` 
#python3 -c "import sys; print(sys.argv[1])" $INSERT
STR2=`python3 -c "import sys; insert=sys.argv[1]; print(insert[int(sys.argv[3]) - int(sys.argv[2]):int(sys.argv[4]) - int(sys.argv[2])+1])" $INSERT $START $EXON2_BEG $EXON2_END ` 

STR="$STR1$STR2"
echo ${STR:0:$READ_LEN}

