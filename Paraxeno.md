
# DATA PROCESS

## Fastq format

![Fastq image](images/fastq_code.png)


#




### 1. Cleanup of reads
```
 rm sample_list.txt; \
 for i in ./fastq/UraH*-dsR_S*_L001_R1_001.fastq.gz; do a=$i; \
 j=${i##./*/}; \
 k=${j%%-dsR_S*_L001_R1_001.fastq.gz}; echo "$i $k"; \
 echo -e "$i\t${i/_R1/_R2}\t$k" >> sample_list.txt; \
 done

 ./Cleanup_FLDS_YT220418.pl -lst sample_list.txt -lib FLDS -outdir cleanup

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



### 3. SEARCH RdRP's contigs

```
 query="UraH02_clc_A1_assembly_cds.pep"

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


### 4. MERGE CANDIDATE CONTIGS
   
```
 query="UraHall_rdrp.fa"
 outname="UraHall_rdrp_cdhit"

 cd-hit-est -d 100 -c 0.97 -aS 0.9 -G 0 -T 20 -M 20000 -i $query -o ${name}_cdhit.fa

```

`OUTPUT:` UraHall_rdrp_cdhit.fa




## ESTIMATE ABUNDANCE OF VIAL CONTIGS



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



### 2 CALCULATE FPKM

```
$ R

 library(tidyverse)


##   RAW COUNT

 dat<-read.table("out_count_tbl.txt", header=T, sep="\t")
 colnames(dat)[1]<-"id"
 dmeta<-read.table("UraHall_clc_A1_ed1_meta.txt", header=T, sep="\t", stringsAsFactors=FALSE)

 dlib<-read.table("UraHall_clc_A1_ed1_data.txt", header=T, sep="\t",stringsAsFactors = FALSE)
 libsize<- dlib$readnum
 names(libsize)<- dlib$sample


 datS <- dat %>% dplyr::mutate( dmeta[match(id,dmeta$con),c("Len","Type")] )


countToFpkm <- function(counts, effLen, Ssize)
{
     N <- Ssize
   # N <- sum(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}



##   CAL FPKM IN clean reads


 geneData<- datS %>% dplyr::select(id,Len)
 geneCount<- datS %>% dplyr::select_if( grepl("Ura", names(.)) | grepl("id", names(.)) )

 countDf<-data.frame()

 j<-0
 countDf<-data.frame()
 countAll_fpkm<-geneData
 meanLen<- 280

 for (i in 2:19) {

   j<-j+1
#   print(i)
#   print(meanLen[j])
   libname<- colnames(geneCount)[i]
   Ssize<- libsize[libname]

   fld<-meanLen
   countDf<-data.frame(count=geneCount[,i],length=geneData[,2])
   colnames(countDf)<- c("count","length")

   countDf$effLength <- ifelse(countDf$length > fld, countDf$length - fld + 1,1)
   countDf$fpkm <- with(countDf, countToFpkm(count, effLength, Ssize))

   countDf_fpkm<- countDf %>% dplyr::select(fpkm)

   countE_fpkm<-data.frame(id=geneData$id,countDf_fpkm)

   countAll_fpkm<-dplyr::right_join(countAll_fpkm, countE_fpkm, by="id")

 }

 colnames(countAll_fpkm)[3:20]<- colnames(geneCount)[2:19]

 write.table(countAll_fpkm,"out_count2fpkm_all.txt", quote=F, row.names=F, col.names=T, appen=F, sep="\t")

    #==> out_count2fpkm_all.txt

```

`OUTPUT:`  out_count2fpkm_all.txt



`OUTPUT:` 
phylogeny-align-to-tree-mafft-fasttree/alignment.qza
phylogeny-align-to-tree-mafft-fasttree/masked_alignment.qza
phylogeny-align-to-tree-mafft-fasttree/tree.qza
phylogeny-align-to-tree-mafft-fasttree/rooted_tree.qza



### 8 EXPORT DATA (BIOM => COUNT TABLE)
 
+   **REPRESENT FAST**  

```
 qiime tools export --input-path rep-seqs-dada2-nochim.qza --output-path output
```

`OUTPUT:`  
 ./output/dna-sequences.fasta  

+   **COUNT TABLE**  
      sample-map.txtは、sample-metadata.tsvのヘッダーに#を付加したもの 
      
```
 perl -pe 's/sample-id/#sample-id/' sample-metadata.tsv > sample-map.txt
 qiime tools export --input-path table-dada2-nochim.qza --output-path output
 biom convert --to-tsv --table-type "OTU table" \
       -i ./output/feature-table.biom -o ./output/feature-count-table.txt \
       -m sample-map.txt
```

`OUTPUT:`  
 ./output/feature-table.biom   
 ./output/feature-count-table.txt






# 参照

Qiime2 を用いた 16S rRNA 菌叢解析
https://qiita.com/keisuke-ota/items/6399b2f2f7459cd9e418

Qiime2やRのPhyloseq、STAMPによる細菌叢解析
https://qiita.com/kuanl/items/a1f98f76ea5a651753f2

初心者の菌叢解析 Qiime2で解析(10) 多様性解析 ~α多様性~
https://note.com/nanaimo_/n/n8543a6008acd



