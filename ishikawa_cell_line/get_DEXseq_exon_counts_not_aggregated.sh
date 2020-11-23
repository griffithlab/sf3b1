#1. Create a GFF file for use with DEXSeq
cd /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/refs
wget ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz
gunzip Homo_sapiens.GRCh38.100.gtf.gz

isub -m 10 -i 'docker(malachig/htseq-0.12.4)'
python /usr/local/DEXSeq/dexseq_prepare_annotation.py --aggregate no
/storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/refs/Homo_sapiens.GRCh38.100.gtf
/storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/Homo_sapiens.GRCh38.100.gff
exit

#2. Create annotations from the GFF to help annotate the practically useless exon count file that will be produced in the next step
cd /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/
grep -v aggregate_gene Homo_sapiens.GRCh38.100.gff | perl -ne 'chomp; $gid="ZZZ"; $eid="ZZZ"; if ($_ =~ /exonic\_part\_number\s+\"(\w+)\".*gene\_id\s+\"(\w+)\"/){$eid=$1; $gid=$2;} print "$gid:$eid\t$_\n"' | sort > Homo_sapiens.GRCh38.100.annotations.tsv

#NOTE: Initially I was understanding that the DEXseq package in R will reassociate counts with annotations using the GFF file.  So the step above is probably not needed.

#3. Prepare exon counts using DEXseq-count and HTSeq-count
#you need to inform the script by specifying the option -s no. If your library preparation protocol reverses the strand (i.e., reads appear on the strand opposite to their gene of origin), use -s reverse. In case of paired-end data, the default (-s yes) means that the read from the first sequence pass is on the same strand as the gene and the read from the second pass on the opposite strand (“forward-reverse” or “fr” order in the parlance of the Bowtie/TopHat manual) and the options -s reverse specifies the opposite case.

#from the collaborator: The complete details of RNA sequencing procedure is included in the method section of manuscript (attached). In brief, the samples were prepared with the random-priming un-stranded Kappa RiboErase kit and sequenced as paired-end 2x150bp reads on an Illumina HiSeq 3000.  The bam files were aligned to GRCh38 with STAR and gene counts were generated with Subread:featureCounts.  Differential expression analysis was conducted with Limma using the voomWithQualityWeights method.

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/Control_Rep1.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/Homo_sapiens.GRCh38.100.gff /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/inputs/Control_Rep1.bam /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/Control_Rep1.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/Control_Rep2.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/Homo_sapiens.GRCh38.100.gff /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/inputs/Control_Rep2.bam /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/Control_Rep2.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/Control_Rep3.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/Homo_sapiens.GRCh38.100.gff /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/inputs/Control_Rep3.bam /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/Control_Rep3.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/SF3B1_Rep1.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/Homo_sapiens.GRCh38.100.gff /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/inputs/SF3B1_Rep1.bam /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/SF3B1_Rep1.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/SF3B1_Rep2.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/Homo_sapiens.GRCh38.100.gff /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/inputs/SF3B1_Rep2.bam /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/SF3B1_Rep2.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/SF3B1_Rep3.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/Homo_sapiens.GRCh38.100.gff /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/inputs/SF3B1_Rep3.bam /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/SF3B1_Rep3.exon_counts.tsv'

#Create clean versions of the DEXseq count files without the statistics lines at the bottom
cd /storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/exon_counts/not_aggregated/
grep -v "^_" Control_Rep1.exon_counts.tsv > Control_Rep1.exon_counts.clean.tsv
grep -v "^_" Control_Rep2.exon_counts.tsv > Control_Rep2.exon_counts.clean.tsv
grep -v "^_" Control_Rep3.exon_counts.tsv > Control_Rep3.exon_counts.clean.tsv
grep -v "^_" SF3B1_Rep1.exon_counts.tsv > SF3B1_Rep1.exon_counts.clean.tsv
grep -v "^_" SF3B1_Rep2.exon_counts.tsv > SF3B1_Rep2.exon_counts.clean.tsv
grep -v "^_" SF3B1_Rep3.exon_counts.tsv > SF3B1_Rep3.exon_counts.clean.tsv
wc -l *.clean.tsv *annotations.tsv


