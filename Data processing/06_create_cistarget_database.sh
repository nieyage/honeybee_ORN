# create honeybee cistarget database 
# ref: https://github.com/aertslab/create_cisTarget_databases/issues/4

conda activate create_cistarget_databases


# Step1: create the fasta file 
# 1: upstream 5kb region of all genes
library("rtracklayer")
tss =  import("/md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/regions/tss.bed")
gene_up5k <- promoters(tss, upstream=5000, downstream=0, use.names=TRUE)
start(gene_up5k)[start(gene_up5k)<0]<- 0

write.table(gene_up5k, file="/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/gene_up5k.bed", sep="\t", quote=F, row.names=F, col.names=F)
awk '{print $1"\t"$2"\t"$3"\t"$6}' gene_up5k.bed > gene_up5k_update.bed

bedtools getfasta -name -fi /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa -bed gene_up5k_update.bed -fo gene_up5k.fasta


# 2: upstream 500bp region + transcript region of all genes
library("rtracklayer")
tss =  import("/md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/regions/tss.bed")
transcripts =  import("/md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/regions/transcripts.bed")
gene_up500bp <- promoters(tss, upstream=500, downstream=0, use.names=TRUE)
last<- c(transcripts,gene_up500bp)
start(last)[start(last)<0]<- 0
write.table(last, file="/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/gene_up500bp_transcript.bed", sep="\t", quote=F, row.names=F, col.names=F)
awk '{print $1"\t"$2"\t"$3"\t"$6}' gene_up500bp_transcript.bed > gene_up500bp_transcript_update.bed
bedtools getfasta -name -fi /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa -bed gene_up500bp_transcript_update.bed -fo gene_up500bp_transcript_update.fasta

# different region in same gene need to add annotation 
#!/usr/bin/env python

from Bio import SeqIO

records = set()
of = open("Numbered_gene_up500bp_transcript.fasta", "w")
for record in SeqIO.parse("gene_up500bp_transcript_update.fasta", "fasta"):
    ID = record.id
    num = 1
    while ID in records:
        ID = "{}#{}".format(record.id, num)
        num += 1
    records.add(ID)
    record.id = ID
    record.name = ID
    record.description = ID
    SeqIO.write(record, of, "fasta")
of.close()


# Step2: create the motif file 

# 1. use honeybee motif (we collected)
# motif file1: /data/R04/liyh526/project/Honeybee/03_CisBP2JASPAR/honeybee_jaspar/CisBP-honeybee.jaspar
# motif file2: /md01/nieyg/ref/10X/Amel_HAv3.1/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt
/data/R04/liyh526/project/Honeybee/03_CisBP2JASPAR/honeybee_jaspar/CisBP-honeybee.jaspar 

# trans GB ID to gene_symbol in our data 
#!/bin/bash
nLine=`wc -l /md01/nieyg/project/honeybee/joint/IDtrans/ID.list | cut -d " " -f 1`
for i in `seq 1 ${nLine}`;
do
        GB_id=`awk '{print $1}' /md01/nieyg/project/honeybee/joint/IDtrans/ID.list | sed -n "${i}p"`
        gene_symbol=`awk '{print $2}' /md01/nieyg/project/honeybee/joint/IDtrans/ID.list | sed -n "${i}p"`
        sed -i "s/${GB_id}/${gene_symbol}/g" CisBP-honeybee.jaspar
done

# trans PFM to CB 
from Bio import motifs
import os
from pathlib import Path
with open("CisBP-honeybee.jaspar", 'rt') as fh:jms = motifs.parse(fh, "JASPAR")
print(jms[0].matrix_id)
print(jms[0].name)
print(jms[0].counts)
print(format(jms[0], "clusterbuster"))

outdir = Path('01_motif_we_collect')
os.mkdir(outdir)
for jm in jms:
        motiffile = outdir.joinpath(str(jm.matrix_id + '.cb'))
        with open(motiffile, 'wt') as fh:
                fh.write(format(jm, "clusterbuster"))
ls 01_motif_we_collect/ | sed 's/\.cb//g' > 01_we_collect_motifs.lst 

# 2. use the public motif database 
wget https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip
# mapping human ID with honeybee gene ID 

# ref:https://github.com/JoGraesslin/Zebrafish_SCENIC

create_cistarget_databases_dir='/md01/nieyg/software/create_cisTarget_databases'
fasta_filename="gene_up500bp_transcript_update.fasta"
motifs_dir="01_motif_we_collect"
motifs_list_filename="01_we_collect_motifs.lst"
db_prefix="up500_trans_cisbp"
nbr_threads=10
genes="#[0-9]+$"

"${create_cistarget_databases_dir}/create_cistarget_motif_databases.py" \
    -f "${fasta_filename}" \
    -M "${motifs_dir}" \
    -m "${motifs_list_filename}" \
    -o "${db_prefix}" \
    -t "${nbr_threads}" \
    -g "${genes}"

ls singletons/ | sed 's/\.cb//g' > 02_v10nr_clust_public.lst 

conda activate create_cistarget_databases

create_cistarget_databases_dir='/md01/nieyg/software/create_cisTarget_databases'
fasta_filename="/data/R02/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/gene_up500bp_transcript_update.fasta"
motifs_dir="/data/R02/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/02_v10nr_clust_public/singletons"
motifs_list_filename="/data/R02/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/02_v10nr_clust_public/02_v10nr_clust_public.lst"
db_prefix="02_v10nr_clust_public_dmel6_to_honeybee"
nbr_threads=5
genes="#[0-9]+$"

nohup  "${create_cistarget_databases_dir}/create_cistarget_motif_databases.py" \
    -f "${fasta_filename}" \
    -M "${motifs_dir}" \
    -m "${motifs_list_filename}" \
    -o "${db_prefix}" \
    -t "${nbr_threads}" \
    -g "${genes}" &

# region vs motif 
awk '{print $1"\t"$2"\t"$3}' /md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/04_scenicplus/02_scATAC/consensus_peak_calling/consensus_regions.bed > consensus_regions.bed
bedtools getfasta  -fi /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa -bed consensus_regions.bed -fo consensus_regions.fasta

# honeybee cisBP motif 
create_cistarget_databases_dir='/md01/nieyg/software/create_cisTarget_databases'
fasta_filename="/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/03_region/consensus_regions.fasta"
motifs_dir="/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/01_motif_we_collect"
motifs_list_filename="/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/01_we_collect_motifs.lst"
db_prefix="01_honeybee_cisBp"
nbr_threads=10
genes="#[0-9]+$"

nohup  "${create_cistarget_databases_dir}/create_cistarget_motif_databases.py" \
    -f "${fasta_filename}" \
    -M "${motifs_dir}" \
    -m "${motifs_list_filename}" \
    -o "${db_prefix}" \
    -t "${nbr_threads}" \
    -g "${genes}" &


## PBS configure 
#PBS -N interproscan 
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=12
#PBS -l mem=16G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/04_scenicplus
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate create_cistarget_databases
create_cistarget_databases_dir='/md01/nieyg/software/create_cisTarget_databases'
fasta_filename="/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/03_region/consensus_regions.fasta"
motifs_dir="/data/R02/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/02_v10nr_clust_public/singletons"
motifs_list_filename="/data/R02/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/02_v10nr_clust_public/02_v10nr_clust_public.lst"
db_prefix="02_v10nr_clust_public_dmel6_to_honeybee"
nbr_threads=12
genes="#[0-9]+$"
 "${create_cistarget_databases_dir}/create_cistarget_motif_databases.py" \
    -f "${fasta_filename}" \
    -M "${motifs_dir}" \
    -m "${motifs_list_filename}" \
    -o "${db_prefix}" \
    -t "${nbr_threads}" \
    -g "${genes}" 

