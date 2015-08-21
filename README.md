# examples
A place to put all those one liners that I always forget

### samtools
convert sam to a sorted bam and index it
```
samtools view -b foo.sam | samtools sort - foo && samtools index foo.bam
```

sort bam in-place and index it
```
samtools sort foo.bam foo && samtools index foo.bam
```
