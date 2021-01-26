#!/bin/bash

if [ "$1" == "" ]
then
    echo "Need FASTA filename as first parameter"
	exit 1
fi

FileName="$1"

if [ "$2" == "" ]
then
    echo "Need sequence similarity threshold as second parameter"
	exit 1
fi

Threshold="$2"

if [ "$3" == "" ]
then
    echo "Result cluster filename is needed as third parameter"
	exit 1
fi

ResClustFileName="$3"

PrefixFileName=$ResClustFileName

TmpFileFastaWithShortIDs="${PrefixFileName}_TmpClustShortID.faa"
TmpFileTmpClustDB="${PrefixFileName}_TmpClustDB"
TmpClustFolder="${PrefixFileName}_TmpClustFolder"
TmpFileClusters="${PrefixFileName}_TmpClust"
TmpFileMMseqsClusters="${PrefixFileName}_TmpClustMMSEQS.tsv"

mkdir $TmpClustFolder
python3 RemoveFASTAIDRedundency.py -f $FileName > $TmpFileFastaWithShortIDs
/home/hazirliver/mmseqs/bin/mmseqs createdb $TmpFileFastaWithShortIDs $TmpFileTmpClustDB # не добавляется mmseqs в PATH, поэтому указан полный путь
/home/hazirliver/mmseqs/bin/mmseqs cluster $TmpFileTmpClustDB $TmpFileClusters $TmpClustFolder --cluster-mode 1 -c 0.1 --min-seq-id $Threshold -e 0.01
/home/hazirliver/mmseqs/bin/mmseqs createtsv $TmpFileTmpClustDB $TmpFileTmpClustDB $TmpFileClusters $TmpFileMMseqsClusters
python3 ConvertOutput.py -f $TmpFileMMseqsClusters > $ResClustFileName

rm -R ${PrefixFileName}_*

