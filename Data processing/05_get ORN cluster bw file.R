### RNA track 

ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak.rds")
metadata <- ORN@meta.data
library(tidyr)
for(i in levels(ORN)){
	print(i);
	barcode <- rownames(metadata[metadata$subcluster==i,])
	barcode_sample<-separate(as.data.frame(barcode),"barcode",c("sample","barcode"),"_")
	barcode <- barcode_sample$barcode
	write.table(barcode,paste0("./08_ORN_cluster_bw/cell_barcode/",i,".txt",sep=""),row.name=F,col.names=F)
}

sed -i 's/"//g' *.txt

# merge bam 
samtools merge -@ 10 merged.bam ~/project/honeybee/data/cellranger/Forager/Forager-NCBI-manually/outs/gex_possorted_bam.unique.bam \
 ~/project/honeybee/data/cellranger/NE/NE-NCBI-manually/outs/gex_possorted_bam.unique.bam \
 ~/project/honeybee/data/cellranger/Nurse/Nurse-NCBI-manually/outs/gex_possorted_bam.unique.bam

samtools index merged.bam

# barcode bam 
ls *.txt> ORN_barcode.list
for file1 in $(<ORN_barcode.list)
do
  /md01/nieyg/software/subset-bam_linux --bam ../merged.bam --cell-barcodes $file1 --cores 10 --out-bam ../barcode_bam/$file1.bam
  # bam to bw 
  echo "make bedGraph"
  # make bedGraph by using bedtools;
  bedtools  genomecov  -bg -split -ibam ../barcode_bam/$file1.bam  > ../barcode_bedGraph/$file1.bedGraph
  echo "make normalized bedGraph"
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl ../barcode_bedGraph/$file1.bedGraph ../barcode_bedGraph/$file1.norm.bedGraph &> ../barcode_bedGraph/$file1.norm.bedGraph.log
  sort -k1,1 -k2,2n ../barcode_bedGraph/$file1.norm.bedGraph > ../barcode_bedGraph/$file1.norm.sorted.bedGraph
  bedGraphToBigWig ../barcode_bedGraph/$file1.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt ../barcode_bw/$file1.norm.bw 
done

#/md01/nieyg/software/subset-bam_linux --bam <FILE> --bam-tag <bam_tag> --cell-barcodes <FILE> --cores <INTEGER> --log-level <log_level> --out-bam <OUTPUT_FILE>
conda deactivate 
nohup bash get_bw.sh &



# merge bam 
nohup samtools merge -@ 10 ATAC-merged.bam ~/project/honeybee/data/cellranger/Forager/Forager-NCBI-manually/outs/atac_possorted_bam.bam \
 ~/project/honeybee/data/cellranger/NE/NE-NCBI-manually/outs/atac_possorted_bam.bam \
 ~/project/honeybee/data/cellranger/Nurse/Nurse-NCBI-manually/outs/atac_possorted_bam.bam &
samtools index ATAC-merged.bam

# barcode bam 
ls *.txt> ORN_barcode.list
for file1 in $(<ORN_barcode.list)
do
  /md01/nieyg/software/subset-bam_linux --bam ../ATAC-merged.bam --cell-barcodes $file1 --cores 10 --out-bam ../barcode_bam/ATAC/$file1.bam
  bedtools  genomecov  -bg -split -ibam ../barcode_bam/ATAC/$file1.bam  > ../barcode_bedGraph/ATAC/$file1.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl ../barcode_bedGraph/ATAC/$file1.bedGraph ../barcode_bedGraph/ATAC/$file1.norm.bedGraph &> ../barcode_bedGraph/ATAC/$file1.norm.bedGraph.log
  sort -k1,1 -k2,2n ../barcode_bedGraph/ATAC/$file1.norm.bedGraph > ../barcode_bedGraph/ATAC/$file1.norm.sorted.bedGraph
  bedGraphToBigWig ../barcode_bedGraph/ATAC/$file1.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt ../barcode_bw/ATAC/$file1.norm.bw 
done
cd cell_barcode
nohup bash ATAC_get_bw.sh &
