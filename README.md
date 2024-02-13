
```
 for a in $(basename  -a ../05_result/*) ; do echo $a ;  create_report ../05_result/$a/$a.bed  --genome hg38     --flanking 10000  --tracks  ../05_result/${a}/${a}.sorted.bam --output ${a}_bed.html; done
```
