#0) look at several of the fastqc results, find out if there are any contaminants or duplicated sequences
#1) $1 should be adapter's sequence
#2) cut adapters/other contaminants (all fastq files should be in ./fastq)
ruby /home/smtb2019/bin/h/h5_gen_cutadapt_template.rb fastq fastq.cutadapt.sh $1
mkdir fastq.trim
chmod +x fastq.cutadapt.sh
./fastq.cutadapt.sh
#3) run fastqc again to ensure that there are no more adapters in the data
#4) deduplicate if needed
mkdir fastq.dedup
ruby /home/smtb2019/bin/h/h6b_gen_deduplicator1_template.rb fastq.trim deduplicator.sh
chmod +x deduplicator.sh
./deduplicator.sh
#5) now you have a trimmed and deduplicated fastq data, that can be mapped using STAR
