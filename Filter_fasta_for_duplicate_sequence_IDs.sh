{\rtf1\ansi\ansicpg1252\cocoartf2636
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fnil\fcharset0 AndaleMono;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;\csgray\c0\c0;}
\margl1440\margr1440\vieww19260\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs26 \cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 #!/bin/bash
\fs24 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\fs28 \cf2 \cb3 \kerning1\expnd0\expndtw0 \CocoaLigature0 \outl0\strokewidth0 \
INPUT_FASTA_FILE=\'93RNAseq_AEC_uninfected_DEG.fa\'94\
OUTPUT_FILE=\'93Filtered.fa\'94\
FINAL_FASTA=\'93RNAseq_AEC_uninfected_DEG.fa\'94\
\
awk 'BEGIN\{RS=">"\}NR>1\{sub("\\n","\\t"); gsub("\\n",""); print RS$0\}' $INPUT_FASTA_FILE | awk '!seen[$1]++' >$OUTPUT_FILE\
\
\
awk 'BEGIN\{RS=">"\}NR>1\{sub("\\n","\\t"); gsub("\\n",""); print RS$0\}' $OUTPUT_FILE | awk '!seen[$1]++' | awk -v OFS="\\n" '\{print $1,$2\}' >$FINAL_FASTA}