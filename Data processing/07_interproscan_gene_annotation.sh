# Part1: interproscan annotate function of honeybee protein sequence 
conda create -n interproscan
conda activate interproscan


 JAVA_HOME=/md01/nieyg/my_interproscan/jdk-11.0.19+7
 PATH=$JAVA_HOME/bin:$PATH
 CLASSPATH=.:$JAVA_HOME/lib/dt.jar:$JAVA_HOME/lib/tools.jar
 export JAVA_HOME
 export PATH
 export CLASSPATH
 source ~/.profile
 java -version


# get the md5 of the databases
version_major=5.59
version_minor=91.0
wget -c http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${version_major}-${version_minor}/interproscan-${version_major}-${version_minor}-64-bit.tar.gz.md5
# get the databases (with core because much faster to download)
wget -c http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${version_major}-${version_minor}/interproscan-${version_major}-${version_minor}-64-bit.tar.gz
# checksum
md5sum -c interproscan-${version_major}-${version_minor}-64-bit.tar.gz.md5
# untar gz
tar xvzf interproscan-${version_major}-${version_minor}-64-bit.tar.gz

## PBS configure 
#PBS -N interproscan 
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=12
#PBS -l mem=16G
source /public/home/nieyg/.bash_profile
cd project/honeybee/honebee-latest-Version/10_fly_honeybee_ID_trans/
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate interproscan
 JAVA_HOME=/md01/nieyg/my_interproscan/jdk-11.0.19+7
 PATH=$JAVA_HOME/bin:$PATH
 CLASSPATH=.:$JAVA_HOME/lib/dt.jar:$JAVA_HOME/lib/tools.jar
 export JAVA_HOME
 export PATH
 export CLASSPATH
 source ~/.profile
 java -version
/md01/nieyg/my_interproscan/interproscan-5.59-91.0/interproscan.sh -i /md01/nieyg/ref/10X/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_cds.fa \
 -f tsv --goterms -dp -cpu 12



qsub /md01/nieyg/project/honeybee/honebee-latest-Version/10_fly_honeybee_ID_trans/interproscan.pbs
# Part2: honyebee to fly ID trans 
# blast: link fly ID and honeybee ID 
# Step1: split the fly nr database 
nohup wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz &
wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz.md5
md5sum -c prot.accession2taxid.FULL.gz.md5 prot.accession2taxid.FULL.gz
gzip -d prot.accession2taxid.FULL.gz

conda create -n split_nr_database
conda activate split_nr_database
conda install taxonkit -c bioconda -y
conda install csvtk -c bioconda -y

tar -zxvf taxdump.tar.gz
mkdir ~/.taxonkit
cp names.dmp ~/.taxonkit
cp nodes.dmp ~/.taxonkit
cp delnodes.dmp ~/.taxonkit
cp merged.dmp ~/.taxonkit

taxonkit list --ids 7227 --indent "" > Dmel6_7227.taxid.txt
cat prot.accession2taxid.FULL |csvtk -t grep -f taxid -P Dmel6_7227.taxid.txt |csvtk -t cut -f accession.version > Dmel6_7227.taxid.acc.txt
wc -l Dmel6_7227.taxid.acc.txt

conda deactivate 
blastdbcmd -db /md01/nieyg/ncbi_database/nr_Dmel6 -entry all -dbtype prot -out nr_Dmel.fa
  makeblastdb [-h] [-help] [-in input_file] [-input_type type]
    -dbtype molecule_type [-title database_title] [-parse_seqids]
    [-hash_index] [-mask_data mask_data_files] [-mask_id mask_algo_ids]
    [-mask_desc mask_algo_descriptions] [-gi_mask]
    [-gi_mask_name gi_based_mask_names] [-out database_name]
    [-max_file_sz number_of_bytes] [-logfile File_Name] [-taxid TaxID]
    [-taxid_map TaxIDMapFile] [-version]



blastdb_aliastool -seqidlist Dmel6_7227.taxid.acc.txt -db /md01/nieyg/ncbi_database/all_species/nr -out nr_Dmel6 -title nr_Dmel6

# Step2: blast honeyebee protein sequence in fly nr database 

####blastp NR database by diamond 
## PBS configure 
#PBS -N diamond 
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=16G
source /public/home/nieyg/.bash_profile
cd project/honeybee/honebee-latest-Version/10_fly_honeybee_ID_trans/
diamond blastx --db /md01/nieyg/ncbi_database/nr_Dmel6 --out diamond_fly_database.out2 \
--query /md01/nieyg/ref/10X/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_cds.fa \
 --max-target-seqs 5 --evalue 1e-5 --threads 24 \
 --block-size 10 --index-chunks 10 --more-sensitive -k 5 
#--outfmt 6 qseqid qlen sseqid qstart qend sstart send evalue pident stitle salltitles \

qsub diamond_blastp.pbs

# Step3: trans the transcript ID 2 gene symbol
# output the third col in file 01_fly_protein_ID.list
awk '{print $3}' diamond_fly_database.out > 01_fly_protein_ID.list

for i in `cat ../01_fly_protein_ID_uniq.list`
do
wget -nd -r -l1 --no-parent https://www.ncbi.nlm.nih.gov/gene/?term=$i
done

for i in `ls 02_fly_protein_html/* `
do
grep -A1 "Symbol" $i | sed -n '2p'| awk -v i=$i 'BEGIN{FS="<|>";OFS="\t"}{print i,$3}' >>ID.list
done

sed -e 's/index.html?term=//g' ID.list > 03_fly_proteinID_2_genesymbol.list
sed -i 's/02_fly_protein_html//g'  03_fly_proteinID_2_genesymbol.list
sed -i 's/^\///' 03_fly_proteinID_2_genesymbol.list

# honeybee transcript ID 2 gene symbol in diamond_fly_database.out 
# in R 
library(Biostrings)
gtf <- rtracklayer::import('/md01/nieyg/ref/10X/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid.sorted.gtf')
fly_proteinID2gene<- read.table("03_fly_proteinID_2_genesymbol.list")
diamond_out<- read.table("diamond_fly_database.out2")
diamond_out$honeybee_gene_name<- gtf[match(diamond_out$V1,gtf$transcript_id),]$gene_name
diamond_out$honeybee_gene_id<- gtf[match(diamond_out$V1,gtf$transcript_id),]$gene_id
diamond_out$fly_gene<- fly_proteinID2gene[match(diamond_out$V2,fly_proteinID2gene$V1),]$V2
diamond_out<- na.omit(diamond_out)
# get the multi mapping row 
duplicated_gene<- diamond_out$honeybee_gene_name[which(duplicated(diamond_out$honeybee_gene_name))]
duplicated_gene_data<- data.frame()
for (i in unique(duplicated_gene)){
    data_tmp<- diamond_out[diamond_out$honeybee_gene_name==i,]
    data_tmp<- data_tmp[order(data_tmp$V3,decreasing=T),]
    duplicated_gene_data<- rbind(duplicated_gene_data,data_tmp[1,])
}
unique_gene_data<- diamond_out[-which(diamond_out$honeybee_gene_name%in%duplicated_gene),]
last<- rbind(duplicated_gene_data,unique_gene_data)
write.csv(last,"04_addgenesymbol_diamond_out.csv")

gene_trans<- last[,c(13,15)]
gene_trans<- gene_trans[!duplicated(gene_trans),]
write.csv(gene_trans,"05_fly2honeybee.csv")

tbl_data<- read.table("/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl",sep="\t")
filter<- tbl_data[which(tbl_data$V6 %in% gene_trans[,2]),]
filter$V6<- gene_trans[match(filter$V6,gene_trans$fly_gene),]$honeybee_gene_name
colnames(filter)<- c("motif_id","motif_name","motif_description","source_name","source_version",
    "gene_name","motif_similarity_qvalue","similar_motif_id","similar_motif_description",
    "orthologous_identity","orthologous_gene_name","orthologous_species","description")
write.table(filter,"/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/motifs-v10nr_clust-nr.flybase-m0.001-o0.0-trans2honeybee.tbl",row.names=F,sep="\t")

sed -i 's/"//g' motifs-v10nr_clust-nr.flybase-m0.001-o0.0-trans2honeybee.tbl










