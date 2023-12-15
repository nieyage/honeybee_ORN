## get the chemoreceptors form the pfam scan results

1. get the protein sequences including 
* Odorant receptor domain: 7tm_6;
* Ionotropic receptor domain: Lig_chan
* Gustatory receptor domain: 7tm_7,Trehalose_recp

`grep "^#" GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_pfam_scan.csv|wc -l `

#back to R 
```
pfam_scan <- read.table("GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_pfam_scan.csv")
colnames(pfam_scan) <- c("seq_id", "alignment_start", "alignment_end", "envelope_start", "envelope_end", "hmm_acc", "hmm_name", "type", "hmm_start", "hmm_end", "hmm_length", "bit_score", "E-value", "significance", "clan")
OR_domain <- c("7tm_6");
GR_domain <- c("7tm_7","Trehalose_recp");
IR_domain <- c("Lig_chan");
OR_scan <- pfam_scan[which(pfam_scan$hmm_name %in% OR_domain),]
GR_scan <- pfam_scan[which(pfam_scan$hmm_name %in% GR_domain),]
IR_scan <- pfam_scan[which(pfam_scan$hmm_name %in% IR_domain),]
OR_transcript <- unique(OR_scan$seq_id)
GR_transcript <- unique(GR_scan$seq_id)
IR_transcript <- unique(IR_scan$seq_id)
```

2. the gene list of chemoreceptors by pfam scan

```
library(Biostrings)
gtf <- rtracklayer::import('/md01/nieyg/ref/10X/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid.sorted.gtf')

OR_gene_id <- unique(gtf[which(gtf$transcript_id %in% OR_transcript),]$gene_name)
GR_gene_id <- unique(gtf[which(gtf$transcript_id %in% GR_transcript),]$gene_name)
IR_gene_id <- unique(gtf[which(gtf$transcript_id %in% IR_transcript),]$gene_name)

```

3. the gene list of chemoreceptors before pfam scan

```
OR_info <- read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/odorant receptor gene_result.txt",sep="\t",header=T)
GR_info <- read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/Gustatory receptor gene_result.txt",sep="\t",header=T)
IR_info <- read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/Ionotropic receptor gene_result.txt",sep="\t",header=T)

diff_OR_id <- OR_gene_id[-grep(paste0(OR_info$Symbol,collapse="|"),OR_gene_id)]
diff_GR_id <- GR_gene_id[-grep(paste0(GR_info$Symbol,collapse="|"),GR_gene_id)]
diff_IR_id <- IR_gene_id[-grep(paste0(IR_info$Symbol,collapse="|"),IR_gene_id)]

diff_OR_id2 <- OR_info$Symbol[-grep(paste0(gsub("-a|-b|-c","",OR_gene_id),collapse="|"),OR_info$Symbol)]
diff_GR_id2 <- GR_info$Symbol[-grep(paste0(gsub("-a|-b|-c","",GR_gene_id),collapse="|"),GR_info$Symbol)]
diff_IR_id2 <- IR_info$Symbol[-grep(paste0(gsub("-a|-b|-c","",IR_gene_id),collapse="|"),IR_info$Symbol)]

pfam_scan[which(pfam_scan$seq_id%in%gtf[gtf$gene_id%in%diff_OR_id2,]$transcript_id),]
gtf[gtf$transcript_id %in% unique(gtf[gtf$gene_id%in%diff_OR_id2,]$transcript_id),]

gtf[gtf$transcript_id %in% unique(gtf[gtf$gene_id%in%diff_OR_id2,]$transcript_id),]


```
### the last chemoreceptor gene list is based on pfam scan results 

4. get the peptide sequence of chemoreceptors 

```
# gene id including the not change gene form gtf 

OR_gene_id <- c(OR_gene_id,"LOC725052","LOC102653615","LOC102653695","LOC102653738","LOC102653858","LOC408517","LOC102653979")
GR_gene_id <- c(GR_gene_id,"LOC413596","LOC726535")
IR_gene_id <- c(IR_gene_id,"LOC100578352")

OR_gtf <- gtf[which(gtf$gene_name %in% OR_gene_id & gtf$type=="transcript"),]
GR_gtf <- gtf[which(gtf$gene_name %in% GR_gene_id & gtf$type=="transcript"),]
IR_gtf <- gtf[which(gtf$gene_name %in% IR_gene_id & gtf$type=="transcript"),]

write.table(OR_transcript,"OR_transcript.txt",row.name=F,col.names=F)
write.table(GR_transcript,"GR_transcript.txt",row.name=F,col.names=F)
write.table(IR_transcript,"IR_transcript.txt",row.name=F,col.names=F)


```
#LOC413596       XM_016915853.2
#LOC726019       XM_026443708.1
LOC725052       GSAman00018
#LOC102653738    XM_026444706.1
#LOC102653695    XM_026444679.1
LOC102653615    XM_026444646.1

seqkit grep -f supply.txt GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_pep.fa -o supply.aa

```
sed -i 's/"//g' OR_transcript.txt
sed -i 's/"//g' GR_transcript.txt
sed -i 's/"//g' IR_transcript.txt
seqkit grep -f OR_transcript.txt GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_pep.fa -o OR_transcript_pep.aa 
seqkit grep -f GR_transcript.txt GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_pep.fa -o GR_transcript_pep.aa 
seqkit grep -f IR_transcript.txt GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_pep.fa -o IR_transcript_pep.aa 
```


5. Collate chemoreceptors information 

```
chemoreceptor_gtf <- c(OR_gtf,GR_gtf,IR_gtf)
chemoreceptor_gtf_data <- as.data.frame(chemoreceptor_gtf)
#add gene type 
chemoreceptor_gtf_data$gene_type<-"NA"
chemoreceptor_gtf_data$gene_type[which(chemoreceptor_gtf_data$gene_name %in% OR_gene_id)]="OR"
chemoreceptor_gtf_data$gene_type[which(chemoreceptor_gtf_data$gene_name %in% GR_gene_id)]="GR"
chemoreceptor_gtf_data$gene_type[which(chemoreceptor_gtf_data$gene_name %in% IR_gene_id)]="IR"
#add gene_domain 
chemoreceptor_gtf_data$gene_domain <- pfam_scan[match(chemoreceptor_gtf_data$transcript_id,pfam_scan$seq_id),]$hmm_name
write.table(chemoreceptor_gtf_data,"chemoreceptor_info_data.csv",row.names=F)

# Add suffixes to duplicate objects
ID_trans <- chemoreceptor_gtf_data[,c(10,12)]
ID_trans$gene_name <- make.unique(ID_trans$gene_name, sep = "_")
write.table(ID_trans,"chemoreceptor_IDtrans.list",row.names=F,col.names=F)

```

6. get the peptide sequence of chemoreceptors named by gene id

#ID trans 

```
sed -i 's/"//g' chemoreceptor_IDtrans.list

#!/bin/bash
nLine=`wc -l chemoreceptor_IDtrans.list | cut -d " " -f 1`
for i in `seq 1 ${nLine}`;
do
        transcript=`awk '{print $1}' chemoreceptor_IDtrans.list | sed -n "${i}p"`
        gene_name=`awk '{print $2}' chemoreceptor_IDtrans.list | sed -n "${i}p"`
        sed -i "s/${transcript}/${gene_name}/g" OR_transcript_pep.aa 
        sed -i "s/${transcript}/${gene_name}/g" GR_transcript_pep.aa 
        sed -i "s/${transcript}/${gene_name}/g" IR_transcript_pep.aa 
done
``` 
