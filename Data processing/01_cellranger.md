## editing ref
1. get fasta and gff3 file from NCBI (Amel_HAv3.1)
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Apis_mellifera/latest_assembly_versions/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Apis_mellifera/latest_assembly_versions/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz
```


2. change the chr name to group 

* the chromsome ID relationship 
```
# trans NCBI_ID to Group123,etc
sed -n '35,121p' GCF_003254395.2_Amel_HAv3.1_assembly_report.txt | awk '{print $1"\t"$8}' > Amel_HAv3.1_assembly_report_1.txt
sed -n '122,211p' GCF_003254395.2_Amel_HAv3.1_assembly_report.txt | awk '{print $1"\t"$7}' > Amel_HAv3.1_assembly_report_2.txt
cat Amel_HAv3.1_assembly_report_1.txt Amel_HAv3.1_assembly_report_2.txt > Amel_HAv3.1_accession_convert.txt
```
```
#!/bin/bash
nLine=`wc -l Amel_HAv3.1_accession_convert.txt | cut -d " " -f 1`
for i in `seq 1 ${nLine}`;
do
        name=`awk '{print $1}' Amel_HAv3.1_accession_convert.txt | sed -n "${i}p"`
        accession=`awk '{print $2}' Amel_HAv3.1_accession_convert.txt | sed -n "${i}p"`
        sed -i "s/${accession}/${name}/g" GCF_003254395.2_Amel_HAv3.1_genomic.fna
done
```
```
#!/bin/bash
nLine=`wc -l Amel_HAv3.1_accession_convert.txt | cut -d " " -f 1`
for i in `seq 1 ${nLine}`;
do
        name=`awk '{print $1}' Amel_HAv3.1_accession_convert.txt | sed -n "${i}p"`
        accession=`awk '{print $2}' Amel_HAv3.1_accession_convert.txt | sed -n "${i}p"`
        sed -i "s/${accession}/${name}/g" GCF_003254395.2_Amel_HAv3.1_genomic.gtf
done
``` 

3. manually annotated chemoreceptors

* back to R
```
library(Biostrings)
gtf <- rtracklayer::import('/md01/nieyg/ref/10X/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.gtf')
```
* get all chemoreceptor gene # get from NCBI 

```
# unique(gtf[grep("odorant receptor",gtf$product),]$gene_id) # not all 
OR_info <- read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/odorant receptor gene_result.txt",sep="\t",header=T)
GR_info <- read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/Gustatory receptor gene_result.txt",sep="\t",header=T)
IR_info <- read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/Ionotropic receptor gene_result.txt",sep="\t",header=T)
chemoreceptor <- c(OR_info$Symbol,GR_info$Symbol,IR_info$Symbol)
chemoreceptor <- chemoreceptor[which(chemoreceptor%in% gtf$gene_id)]
chemoreceptor_gtf <- gtf[which(gtf$gene_id%in%chemoreceptor),]
chemoreceptor_gtf_exon <- chemoreceptor_gtf[which(chemoreceptor_gtf$type=="exon"),]
chemoreceptor_exon <- as.data.frame(table(chemoreceptor_gtf_exon$gene_id))
chemoreceptor_exon <- chemoreceptor_exon[order(chemoreceptor_exon$Freq,decreasing=T),]
chemoreceptor_exon$type <- "receptor";
chemoreceptor_exon$type[which(chemoreceptor_exon$Var1%in% OR_info$Symbol)]="OR";
chemoreceptor_exon$type[which(chemoreceptor_exon$Var1%in% IR_info$Symbol)]="IR";
chemoreceptor_exon$type[which(chemoreceptor_exon$Var1%in% GR_info$Symbol)]="GR";
write.csv(chemoreceptor_exon,"chemoreceptor_exon_number.csv")

```
* annotate chemoreceptor manually by bulk RNAseq data and ATACseq data 

4. modify the gtf file 

```
# check the integrality of gff3 file
agat_convert_sp_gxf2gxf.pl -gff GCF_003254395.2_Amel_HAv3.1_genomic_chemoreceptor_update.gff3  -o GCF_003254395.2_Amel_HAv3.1_genomic_chemoreceptor_update_fixed.gff3
# gff3 to gtf 
gffread GCF_003254395.2_Amel_HAv3.1_genomic_chemoreceptor_update_fixed.gff3 -T -o GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED.gtf

# gene_id = gene_name
library(Biostrings)
gtf <- rtracklayer::import('/md01/nieyg/ref/10X/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED.gtf')
gene<-c("LOC726834","LOC100579011","Or30","Gr10","LOC100576940","LOC100576992","LOC724202","LOC724462","LOC724590","Or63","LOC102656567","LOC100578751","LOC724673","LOC100577888")
for(i in gene){
        print(i);
        gtf[grep(i,gtf$gene_name),]$gene_id=gtf[grep(i,gtf$gene_name),]$gene_name
}
gtf[which(is.na(gtf$gene_id))]$gene_id=gtf[which(is.na(gtf$gene_id))]$gene_name
rtracklayer::export(gtf,"GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid.gtf" ,format = "gtf")
```

5. get the protein sequence 

* get cds (all gene)

`gffread GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid.gtf -g GCF_003254395.2_Amel_HAv3.1_genomic.fna -x GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_cds.fa`

* get pep (all gene)

`gffread GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid.gtf -g GCF_003254395.2_Amel_HAv3.1_genomic.fna -y GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_pep.fa`

* protein domain scan by pfam_scan

`nohup pfam_scan.pl -outfile GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_pfam_scan.csv -fasta GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_pep.fa -dir /data/R02/nieyg/software/pfamdb &`


[Tutorial pfam_scan](https://github.com/aziele/pfam_scan)


## get cell * gene/peak matrix by CellRanger 
1. create a cellranger-arc-compatible reference by cellranger

* "records for gene_id Mir2-1 are not contiguous in the file"  This is probably because the gtf is sorted by location, and if a gene overlaps another, then the gene records aren't going to be contiguous in the gtf.

(ref link) [https://github.com/10XGenomics/cellranger/issues/133]

* fix_gtf.py
```
#!/usr/bin/env python3

import argparse
import csv
import natsort
import pandas

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument(
  '-i',
  dest = 'input',
  required = True,
  help = 'Input gtf'
)
parser.add_argument(
  '-o',
  dest = 'output',
  required = True,
  help = 'Output sorted gtf'
)
args = parser.parse_args()

# Read gtf input
gtf_columns = {
  'chromosome': 'str',
  'source': 'str',
  'feature': 'str',
  'start': 'uint64',
  'end': 'uint64',
  'score': 'str',
  'strand': 'str',
  'frame': 'str',
  'attribute': 'str'
}

gtf = pandas.read_csv(
  args.input,
  sep = '\t',
  comment = '#',
  names = gtf_columns.keys(),
  dtype = gtf_columns
)

# Get gene id and transcript id for each row
gtf['gene'] = gtf['attribute'].str.extract(r'gene_id "(.+?)"')
gtf['transcript'] = gtf['attribute'].str.extract(r'transcript_id "(.+?)"')

# Get genes start and end positions
genes = gtf[gtf['feature'] == 'gene'][
  ['chromosome', 'strand', 'gene', 'start', 'end']
].rename(
  columns = {
    'start': 'gene_start',
    'end': 'gene_end'
  }
).set_index(
  ['chromosome', 'strand', 'gene']
)

# Get transcripts start and end positions
transcripts = gtf[gtf['feature'] == 'transcript'][
  ['chromosome', 'strand', 'transcript', 'start', 'end']
].rename(
  columns = {
    'start': 'transcript_start',
    'end': 'transcript_end',
  }
).set_index(
  ['chromosome', 'strand', 'transcript']
)

# Add gene and transcript start and end positions to each row
gtf = gtf.set_index(
  ['chromosome', 'strand', 'gene']
).merge(
  genes,
  how = 'left',
  on = ['chromosome', 'strand', 'gene']
).reset_index().set_index(
  ['chromosome', 'strand', 'transcript']
).merge(
  transcripts,
  how = 'left',
  on = ['chromosome', 'strand', 'transcript']
).reset_index()

# Sort rows
gtf = gtf.sort_values(
  by = [
    'chromosome',
    'strand',
    'gene_start',
    'gene_end',
    'gene',
    'transcript_start',
    'transcript_end',
    'transcript',
    'feature',
    'start',
    'end'
  ],
  key = lambda x: (
    [0 if i == 'gene' else 1 if i == 'transcript' else 2 for i in x]
    if x.name == 'feature'
    else natsort.natsort_key(x)
  )
)

# Write gtf to output
gtf.to_csv(
  args.output,
  sep = '\t',
  columns = gtf_columns.keys(),
  header = False,
  index = False,
  quoting = csv.QUOTE_NONE,
  float_format = '%.10g'
)
```
`./fixgtf.py -i /md01/nieyg/ref/10X/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid.gtf -o /md01/nieyg/ref/10X/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid.sorted.gtf`

```
organism: "Apis mellifera"
genome: ["Amel_HAv3_1"]
input_fasta: ["/md01/nieyg/ref/10X/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna"]
input_gtf: ["/md01/nieyg/ref/10X/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid.sorted.gtf"]
input_motifs: "/md01/nieyg/ref/10X/Amel_HAv3.1/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt"
non_nuclear_contigs: ["MT"]
```
`nohup cellranger-arc mkref --config=mkref.contig &`

2.  Run cellranger-arc count
* Forager-NCBI-manually.pbs
```
#PBS -N Forager
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=30000m

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/honeybee/data/cellranger/Forager
cellranger-arc count --id=Forager-NCBI-manually\
                       --reference=/md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1 \
                       --libraries=/md01/nieyg/project/honeybee/data/cellranger/Forager/libraries-add.csv \
                       --localcores=24 \
                       --localmem=30

```

* NE-NCBI-manually.pbs
```
#PBS -N NE
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=30000m

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/honeybee/data/cellranger/NE
cellranger-arc count --id=NE-NCBI-manually\
                       --reference=/md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1 \
                       --libraries=/md01/nieyg/project/honeybee/data/cellranger/NE/libraries-add.csv \
                       --localcores=24 \
                       --localmem=30

```


* Nurse-NCBI-manually.pbs
```
#PBS -N Nurse
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=30000m

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/honeybee/data/cellranger/Nurse
cellranger-arc count --id=Nurse-NCBI-manually\
                       --reference=/md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1 \
                       --libraries=/md01/nieyg/project/honeybee/data/cellranger/Nurse/libraries-add.csv \
                       --localcores=24 \
                       --localmem=30

```
