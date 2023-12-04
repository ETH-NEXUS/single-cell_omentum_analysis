#!/bin/bash
conda activate rna_velocity

myinputfolder=$1
mysamplename=$2


mkdir a-velocyto/${mysamplename}
velocyto run --bcfile ${myinputfolder}/${myinputfolder}.barcodes.tsv --outputfolder velocyto/${mysamplename} --sampleid ${mysamplename} --mask grch38_rmsk.gtf -@ 27 ${myinputfolder}/${myinputfolder}_possorted_genome_bam.bam refdata-cellranger-GRCh38-3.0.0/genes/genes.filtered.gtf 
