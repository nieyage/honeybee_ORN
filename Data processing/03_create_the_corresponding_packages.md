## create the corresponding packages: BSgenome and 

1. BSgenome

* trans fasta to 2bit 
faToTwoBit /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.2bit

* bulid seed file
```
Package: BSgenome.Amel.HAv3.1.update.chemoreceptor
Title: Full genome sequences for update_chemoreceptor
Description: Full genome sequences for honeybee and stored in Biostrings objects.
Version: 1.0.0
organism: Apis mellifera (Honey Bee)
common_name: Honey Bee
organism_biocview: Apis mellifera
provider: NCBI
genome:Amel_HAv3.1
release_date: 2021
BSgenomeObjname: Amel
seqfile_name: genome.2bit

```

* bulid genome packages
```
library(BSgenome)
library(Biostrings)
forgeBSgenomeDataPkg("genome.seed", seqs_srcdir="/md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta", 
  destdir="/md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/BSgenome_Amel_HAv3_1/", verbose=TRUE)
```

```
cd /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/BSgenome_Amel_HAv3_1
R CMD build BSgenome.Amel.HAv3.1.update.chemoreceptor
R CMD check BSgenome.Amel.HAv3.1.update.chemoreceptor_1.0.0.tar.gz
R CMD INSTALL BSgenome.Amel.HAv3.1.update.chemoreceptor_1.0.0.tar.gz
```

library(BSgenome.Amel.HAv3.1.update.chemoreceptor)


# make the Apis_mellifera.OrgDb packages
```
#####GO term######
library(AnnotationHub)
library(biomaRt)
library(dplyr)
library(goseq)
library(DOSE)
library(stringr)
library(clusterProfiler)
library(GenomicRanges)
hub <- AnnotationHub::AnnotationHub() 
query(hub, "Apis mellifera") 
hub
Apis_mellifera.OrgDb <- hub[["AH102515"]]
columns(Apis_mellifera.OrgDb)
Apis_mellifera.OrgDb
head(keys(Apis_mellifera.OrgDb,keytype = "SYMBOL"),5) 
keytypes(Apis_mellifera.OrgDb)
keys(Apis_mellifera.OrgDb)[1:10] 
bitr(keys(Apis_mellifera.OrgDb)[1:10], 'ENTREZID', "SYMBOL", Apis_mellifera.OrgDb) 
bitr(keys(Apis_mellifera.OrgDb)[1:2], 'ENTREZID', c("SYMBOL","REFSEQ", "GO", "ONTOLOGY"), Apis_mellifera.OrgDb) 
saveDb(Apis_mellifera.OrgDb, "/md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/Apis_mellifera_AH102515.OrgDb")
```
