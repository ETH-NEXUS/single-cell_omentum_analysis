#!/bin/bash

# prepare fastq files for novaseq preprocessing; create symlinks that correct the wrong novaseq sample name (change lane-specific extra tags)

outFastqDirec=$1
inSampleNames=$2
downloadInfo=$3

for name in $(less $inSampleNames)
do
	echo $name
	outPath=${outFastqDirec}${name}
	echo $outPath

	for i in {1..4}
	do
		echo $i
		patternR1=L00${i}_R1
		patternR2=L00${i}_R2
		patternI1=L00${i}_I1

		myFileR1=$(grep $name $downloadInfo | grep "fastq.gz" | grep $patternR1)
		myFileR2=$(grep $name $downloadInfo | grep "fastq.gz" | grep $patternR2)
		myFileI1=$(grep $name $downloadInfo | grep "fastq.gz" | grep $patternI1)

		sTag=$(ls $myFileR1 | rev | cut -d "_" -f 4 | rev)
		echo $sTag
		echo $myFileR1
		echo $myFileR2
		echo $myFileI1

		cp $myFileR1 ${outPath}/${name}_${sTag}_L00${i}_R1.fastq.gz
		cp $myFileR2 ${outPath}/${name}_${sTag}_L00${i}_R2.fastq.gz
		cp $myFileI1 ${outPath}/${name}_${sTag}_L00${i}_I1.fastq.gz
	done
done

	
