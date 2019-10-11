#0) look at several of the fastqc results, find out if there are any contaminants or duplicated sequences
#1) $1 should be *_1 adapter's sequence, $2 should be *_2 adapter's sequence
#2) cut adapters/other contaminants (all fastq files should be in ./fastq)
mkdir fastq.trim
# the ribo-seq way (only for unpaired-end fastq's):
ruby /home/smtb2019/bin/h/h5_gen_cutadapt_template.rb fastq fastq.cutadapt.sh $1
chmod +x fastq.cutadapt.sh
./fastq.cutadapt.sh
# for paired-end fastq's:
for name in ./fastq/*_1.fastq.gz
do
  name2="$(echo $name | sed -e "s/_1.fastq.gz/_2.fastq.gz/g")"
  cutadapt -a $1 -A $2 -j 4 --minimum-length 20 -q 20 -o ./fastq.trim/$(basename $name) -p ./fastq.trim/$(basename $name2) $name $name2 
done
# for unpaired-end fastq's:
for name in ./fastq/*.fastq.gz
do
  cutadapt -a $1 -j 4 --minimum-length 20 -q 20 -o ./fastq.trim/$(basename $name) $name 
done
#3) run fastqc again to ensure that there are no more adapters in the data
#4) deduplicate if needed
mkdir fastq.dedup
ruby /home/smtb2019/bin/h/h6b_gen_deduplicator1_template.rb fastq.trim deduplicator.sh
chmod +x deduplicator.sh
./deduplicator.sh
#5) now you have a trimmed and deduplicated fastq data, that can be mapped using STAR
