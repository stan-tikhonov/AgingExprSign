#!/bin/sh
#1) Prepare output directories
mkdir star_output
mkdir count_result_star

#2) run STAR and FeatureCounts for every sample
cd fastq.trim
for var1 in *_1.fastq.gz
do
  var2="$(echo $var1 | sed -e "s/_1.fastq.gz/_2.fastq.gz/g")"
  outcount="$(echo $var1 | sed -e "s/_1.fastq.gz/.count/g")"
  dirname="$(echo $var1 | sed -e "s/_1.fastq.gz//g")"
  mkdir ../star_output/$dirname
  STAR --runThreadN 10 --genomeDir ~/genomes/mus_musculus/star_indexed --readFilesIn $var1 $var2 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ../star_output/$dirname/
  featureCounts -F GTF -a ~/genomes/mus_musculus/gtf/Mus_musculus.GRCm38.98.gtf -g gene_id -t exon -o ../count_result_star/$outcount -s 0 -p -C -B -T 10 ../star_output/$dirname/Aligned.sortedByCoord.out.bam
done
