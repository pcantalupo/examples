# bowtie2 align and save output as unsorted bam
bowtie2 ... | samtools view -Sb - > foo.bam

# bowtie2 align and save output as sorted bam (output file will be 'foo.bam')
bowtie2 ... | samtools view -Sb - | samtools sort - foo
bowtie2 ... | samtools view -Sb - | samtools sort -T foo -o foo.bam # (samtools v1.3; works without -T option on command line but not on SLURM)

# sort bam in-place and index it
samtools sort foo.bam foo && samtools index foo.bam
samtools sort -T foo -o foo.bam foo.bam && samtools index foo.bam  # (samtools v1.3; works without -T option on command line but not on SLURM)
