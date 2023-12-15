

# published data
/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_nameing/Supp_AmOrGrs.fasta
# our OR aa 
/md01/nieyg/ref/10X/Amel_HAv3.1/OR_transcript_pep2.aa

blastp -query /md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_nameing/Supp_AmOrGrs.fasta -subject /md01/nieyg/ref/10X/Amel_HAv3.1/all_OR_pep.aa \
-out OR_gene_mapping.txt -outfmt 6 


conda activate r4-base
R

data<- read.table("OR_gene_mapping.txt")

filtered_data<- data[data$V3>90,1:3]

chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)

OR_gene_name<- data.frame(OR_gene= OR_gene,GR2006= OR_gene, NCBI=OR_gene)
OR_gene_name$GR2006<-"NA"
OR_gene_name$NCBI<- "NA"
OR_gene_name$pident<- 0
tmp_data_set<- data.frame()
for(i in 1:nrow(OR_gene_name)){
	gene<- OR_gene_name$OR_gene[i]
	temp<- data[data$V2==gene,];
	max_pident<- max(temp$V3);
	if(max_pident>80){
	max_gene<- temp[temp$V3==max_pident,]$V1;
	if(length(max_gene)==1){
	OR_gene_name[i,]$GR2006=max_gene
	OR_gene_name[i,]$pident=max_pident
      }else{
        OR_gene_name[i,]$GR2006=max_gene[1]
	    OR_gene_name[i,]$pident=max_pident;
	    tmp_data<- OR_gene_name[i,]
	    tmp_data$GR2006=max_gene[2]
	    tmp_data_set<- rbind(tmp_data_set,tmp_data)
      }
	}

}

OR_gene_name<- rbind(OR_gene_name,tmp_data_set)

write.csv(OR_gene_name,"OR_gene_name.csv")

OR_gene_name<- read.csv("OR_gene_name.csv")

gtf <- rtracklayer::import('/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz')
gtf$gene_biotype<-"protein_coding"
gene.coords <- gtf[gtf$gene_biotype == 'protein_coding']

gtf_data<- as.data.frame(gene.coords[match(OR_gene_name$OR_gene,gene.coords$gene_name),])
last<- cbind(OR_gene_name[,c(2,4,3,5,6)],gtf_data[,c(1:7,10:13)])

write.csv(last,"OR_gene_naming_result.csv")

cp /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz .

awk -F',' '{print $2, $6}' OR_gene_naming_result.csv > gene_replacements.txt

sed -i 's/"//g' gene_replacements.txt

while read -r old_gene new_gene; do
  sed -i "s/${old_gene}/${new_gene}/g" genes.gtf
done < gene_replacements.txt




/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_nameing/OR_gene_naming_result.csv

