#!/bin/bash

ctr=1
fTotal=$(ls ../fastq/*.gz -1 | wc -l)

for file in ../fastq/*.gz
   do
   if test -f "$file"
   then
      fPrefix=$(echo "$file" | sed 's|.*/||' | cut -f1 -d ".")
      echo mapping and sorting individual "$ctr" of "$fTotal"
      bwa mem -t 24 hcomma "$file" | \
      samtools view -u -b - | \
      samtools sort -l0 -@24 -o "$fPrefix.sorted.bam"
      ((ctr++))
   fi
done
for sBam in *.sorted.bam
   do
   if test -f "$sBam"
   then
      samtools index -@24 "$sBam"
   fi
done
```