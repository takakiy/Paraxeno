
# DATA PROCESSING & RdRP SEARCH

## Assembly

### 1. CLEANUP of reads
```
 rm sample_list.txt; \
 for i in ./fastq/UraH*-dsR_S*_L001_R1_001.fastq.gz; do a=$i; \
 j=${i##./*/}; \
 k=${j%%-dsR_S*_L001_R1_001.fastq.gz}; echo "$i $k"; \
 echo -e "$i\t${i/_R1/_R2}\t$k" >> sample_list.txt; \
 done

 ./Cleanup_FLDS.pl -lst sample_list.txt -lib FLDS -outdir cleanup

```

`OUTPUT:` ./cleanup/UraH02-A_PP_R1.fq


### 2. ASSEMBLY

> [!Note]
> Assembly using CLC genomics workbench

```
 read: UraH02-A_PP_R1.fq.gz

 WORD: 33
 BUBBLE: 300
 MAP SIMILARITY: 0.9
 MAP COVERATE: 0.9
 

```
> [!Note]
> Reform of output

```

 infile="UraH02-A_PP"; \
 iname="UraH02-A_PP"; \
 oname="UraH02_clc_A1"; \
 perl -pe 'if(/>'${iname}'_contig_(\d+)/){ $no=sprintf("%.4u",$1); s/>'${iname}'_contig_(\d+)/>'${oname}'_ctg${no}/ }' ${infile}_assembly.fa > ${oname}_assembly.fna; \
 perl -pe 'if(/'${iname}'_contig_(\d+)_mapping/){$no=sprintf("%.4u",$1); s/'${iname}'_contig_(\d+)_mapping/'${oname}'_ctg${no}/}' ${infile}_assembly.txt > ${oname}_assembly.txt; \
 gc_contentSkew.pl -if ${oname}_assembly.fna -p gc; \
 ~/Desktop/work_genome/bin/adding_info_list.pl -i ${oname}_assembly.txt -a outgc -key 0 0 -val 2 3; \
 mv out_${oname}_assembly.txt ${oname}_assembly_rd.txt; \
 echo ${oname}_assembly_rd.txt; \


```

`OUTPUT:` UraH02_clc_A1_assembly.fna

## RdRP SEARCH

### 1. SEARCH RdRP's contigs

> [!Note]
> Search ORFs in contigs

```
 ref="UraH02_clc_A1_assembly.fna"; \
 name="UraH02_clc_A1";

 $HOME/biotools/local/genome/Prodigal/prodigal -i $ref -m -a ${name}.orfs.faa -d ${name}.orfs.fna -o ${name}.gbk -f gbk -p meta -q -g 11

```
`OUTPUT:` UraH02_clc_A1.orfs.faa


> [!Note]
> Search RdRP

```
 query="UraH02_clc_A1.orfs.faa"

 ## Pfam
 HMM1="Pfam-A.hmm" 
 outname="Ura02_clc_A1_pfam"
 hmmscan --cpu 2 --cut_ga -o $outname --domtblout $outname $HMM1 $query

 ## Zayed et al 
 HMM2="all_virus_RdRP_profiles.hmm"  ## 
 outname="Ura02_clc_A1_rdrpsci"
 hmmscan --cpu 2 -o $outname --domtblout $outname $HMM1 $query

 ## NeoRdRp
 HMM3="NeoRdRp.hmm"  ## 
 outname="Ura02_clc_A1_neordrp"
 hmmscan --cpu 2 -o $outname --domtblout $outname $HMM1 $query
```

`OUTPUT:` Ura02_clc_A1_pfam_tbl.txt


### 2. MERGE CANDIDATE CONTIGS

   
```
 query="UraHall_rdrp.fa"
 outname="UraHall_rdrp_cdhit"

 cd-hit-est -d 100 -c 0.97 -aS 0.9 -G 0 -T 20 -M 20000 -i $query -o ${name}_cdhit.fa

```

`OUTPUT:` UraHall_rdrp_cdhit.fa


## ESTIMATE ABUNDANCE OF VIRAL CONTIGS

### 1. MAPPING

```
 read="UraH02-A_PP_R1.fq.gz"

 query="UraHall_rdrp_cdhit.fa"

 alname="UraH02-A"

 bbmap.sh -Xmx20g in1=$read in2=${read/_R1/_R2} ref=$ref out=${alname}.sam minid=0.97 maxindel=3

 samtools view -@ 20 -hb -o ${alname}.bam ${alname}.sam
 samtools sort -@ 20 -o ${alname}.sort.bam ${alname}.bam
 samtools index ${alname}.sort.bam

   # CALCULATION COVERGAE
 samtools coverage ${alname}.sort.bam > out_depth_${alname}.txt

```

`OUTPUT:`   out_depth_UraH02-A.txt


### 2. CALCULATE FPKM

```
$ R

 library(tidyverse)

 dat<-read.table("out_depth_UraH02-A.txt", header=T, sep="\t")
 colnames(dat)<-c("rname","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq")

 datS <- dat %>% dplyr::mutate(id=rname, Len=endpos, count=numreads) %>% dplyr::select(id,Len,count)
 Ssize<- sum(dat$numreads)


 countToFpkm <- function(counts, effLen, Ssize)
 {
     N <- Ssize
   # N <- sum(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
 }



 ## CALCULATION FPKM

 geneData<- datS %>% dplyr::select(id,Len)
 geneCount<- datS %>% dplyr::select(id,count)

 countDf<-data.frame(id=geneCount$id)
 fld<- 280

 countDf<-data.frame(count=geneCount$count,length=geneData$Len)
 colnames(countDf)<- c("count","length")

 countDf$effLength <- ifelse(countDf$length > fld, countDf$length - fld + 1,1)
 countDf$fpkm <- with(countDf, countToFpkm(count, effLength, Ssize))


 write.table(countDf,"out_count2fpkm.txt", quote=F, row.names=F, col.names=T, appen=F, sep="\t")


```

`OUTPUT:`  out_count2fpkm.txt



