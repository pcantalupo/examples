# bowtie2 align and save output as unsorted bam
bowtie2 ... | samtools view -Sb - > foo.bam

# bowtie2 align and save output as sorted bam (output file will be 'foo.bam')
bowtie2 ... | samtools view -Sb - | samtools sort - foo
bowtie2 ... | samtools view -Sb - | samtools sort -o foo.bam # (samtools v1.3)

# convert sam to a sorted bam and index it
samtools view -Sb foo.sam | samtools sort - foo && samtools index foo.bam
samtools view -Sb foo.sam | samtools sort -o foo.bam && samtools index foo.bam # (samtools v1.3)

# sort bam in-place and index it
samtools sort foo.bam foo && samtools index foo.bam
samtools sort -o foo.bam foo.bam && samtools index foo.bam  # (samtools v1.3)


