# Tutorial 1: Training an AtacWorks model

!wget https://api.ngc.nvidia.com/v2/models/nvidia/atac_dsc_atac_lowcellcount_1m_48m_50_2400/versions/0.3/files/train_data/noisy_data/dsc.1.Mono.50.cutsites.smoothed.200.bw
!wget https://api.ngc.nvidia.com/v2/models/nvidia/atac_dsc_atac_lowcellcount_1m_48m_50_2400/versions/0.3/files/train_data/clean_data/dsc.Mono.2400.cutsites.smoothed.200.bw
!wget https://api.ngc.nvidia.com/v2/models/nvidia/atac_dsc_atac_lowcellcount_1m_48m_50_2400/versions/0.3/files/train_data/clean_data/dsc.Mono.2400.cutsites.smoothed.200.3.narrowPeak

# Step 3: Train and validate a model using the parameters in the given config files
# training set
# validation set:chr2
# holdout set :chr10
atacworks=/md01/nieyg/software/AtacWorks

python $atacworks/scripts/main.py train \
    --config /md01/nieyg/software/AtacWorks/configs/train_config.yaml \
    --noisybw /md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/0_tutorial/dsc.1.Mono.50.cutsites.smoothed.200.bw \
    --cleanbw /md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/0_tutorial/dsc.Mono.2400.cutsites.smoothed.200.bw \
    --cleanpeakfile /md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/0_tutorial/dsc.Mono.2400.cutsites.smoothed.200.3.narrowPeak \
    --genome $atacworks/data/reference/hg19.auto.sizes \
    --val_chrom chr2 \
    --holdout_chrom chr10 \
    --out_home "./" \
    --exp_name "atacworks_train" \
    --distributed 

# Tutorial 2: Using a trained AtacWorks model to denoise ATAC-seq data and call peaks.
mkdir models
wget -P 1_models https://api.ngc.nvidia.com/v2/models/nvidia/atac_dsc_atac_lowcellcount_1m_48m_50_2400/versions/0.3/files/models/model.pth.tar
wget https://api.ngc.nvidia.com/v2/models/nvidia/atac_dsc_atac_lowcellcount_1m_48m_50_2400/versions/0.3/files/test_data/noisy_data/dsc.1.NK.50.cutsites.smoothed.200.bw
conda activate python37

python $atacworks/scripts/main.py denoise \
    --noisybw dsc.1.NK.50.cutsites.smoothed.200.bw \
    --genome $atacworks/data/reference/hg19.auto.sizes \
    --weights_path ../1_models/model.pth.tar \
    --out_home "./" \
    --exp_name "atacworks_denoise" \
    --distributed \
    --num_workers 0 \
    --config /md01/nieyg/software/AtacWorks/configs/model_structure.yaml




# in our honeybee data 
# all ORN data filtered 
# in R:
# all cells with Orco 
library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Amel.antenan)
library(patchwork)
set.seed(1234);
library(pheatmap)
ORN <- readRDS("./05_ORN_cluster2/01_first_cluster/ORN_integrated_antenna_first_cluster.rds")
ORN_last <- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak.rds")
ORN_filtered_cells<- setdiff(colnames(ORN),colnames(ORN_last))
library(tidyr)
barcode <- ORN_filtered_cells
barcode_sample<-separate(as.data.frame(barcode),"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
write.table(barcode,"/md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/2_honeybee/ORN_filtered_cells.txt",row.name=F,col.names=F)
# random select 50 cell to get the bw file form filtered ORN 
random_50cell<- sample(barcode,50)
write.table(random_50cell,"/md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/2_honeybee/random_50cell.txt",row.name=F,col.names=F)
sed -i 's/"//g' *.txt

# barcode bam 
ls *.txt> ORN_barcode.list
for file1 in $(<ORN_barcode.list)
do
  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/ATAC-merged.bam --cell-barcodes $file1 --cores 30 --out-bam ./$file1.bam
  bedtools  genomecov  -bg -split -ibam ./$file1.bam  > ./$file1.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl ./$file1.bedGraph ./$file1.norm.bedGraph &> ./$file1.norm.bedGraph.log
  sort -k1,1 -k2,2n ./$file1.norm.bedGraph > ./$file1.norm.sorted.bedGraph
  bedGraphToBigWig ./$file1.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt ./$file1.norm.bw 
done
nohup bash ATAC_get_bw.sh &

# bam call peak 
nohup macs2 callpeak -t ORN_filtered_cells.txt.bam -f BAMPE -g 2.7e+09 --keep-dup all -n ORN_filtered_cells --outdir . &



## PBS configure 
#PBS -N Atacworks_train
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=4
#PBS -l mem=32G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/2_honeybee
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate python37

honebee_input=/md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/2_honeybee
atacworks=/md01/nieyg/software/AtacWorks
python $atacworks/scripts/main.py train \
    --config $atacworks/configs/train_config.yaml \
    --noisybw $honebee_input/random_50cell.txt.norm.bw \
    --cleanbw $honebee_input/ORN_filtered_cells.txt.norm.bw \
    --cleanpeakfile $honebee_input/ORN_filtered_cells.narrowPeak \
    --genome /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt \
    --val_chrom Group2 \
    --holdout_chrom Group10 \
    --out_home "./" \
    --exp_name "atacworks_train_OSN_filtered" \
    --distributed 

## PBS configure 
#PBS -N Atacworks_denoise
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=6
#PBS -l mem=32G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/2_honeybee
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate python37

honebee_input=/md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/2_honeybee
atacworks=/md01/nieyg/software/AtacWorks
ATAC_bw=/md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/barcode_bw/ATAC

ls $ATAC_bw/*.bw |cut -d "/" -f 10 > ORN_ATAC_bw.list
for file1 in $(<ORN_ATAC_bw.list)
do
echo $file1
python $atacworks/scripts/main.py denoise \
    --noisybw $ATAC_bw/$file1 \
    --genome /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt \
    --weights_path $honebee_input/atacworks_train_OSN_filtered_latest/model_best.pth.tar \
    --out_home "./" \
    --exp_name "atacworks_denoise" \
    --distributed \
    --batch_size 10 \
    --num_workers 6 \
    --config /md01/nieyg/software/AtacWorks/configs/model_structure.yaml
done;


# plot the track plot by refer bw:
library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
set.seed(1234)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
# remove nopower cluster,then plot UMAP
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak.rds")
cluster="p2:21" #p2:21
obj<- subset(ORN,idents=cluster)
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot
DefaultAssay(obj)<-"peaks_ORN_subcluster"
## plot region 
start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
ranges.show <- paste(seq,start,end,sep="-")

listbw<-list("p2:21" = "/md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/3_OSN_denoise/atacworks_denoise_all_ORN_cluster/01_bw/p2:21_infer.track.bw", 
    "other" = "/md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/2_honeybee/random_50cell.txt.norm.bw")

pdf("./test_data.pdf")
BigwigTrack(
  StringToGRanges(ranges.show),
  listbw,
  smooth = 200,
  extend.upstream = 500,
  extend.downstream = 500,
  type = "coverage",
  y_label = "bigWig",
  bigwig.scale = "common",
  ymax = NULL,
  max.downsample = 30,
  downsample.rate = 0.8
)
dev.off()


# RT OR cluster bw and random 100 cell for other bw 
library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Amel.antenan)
library(patchwork)
set.seed(1234);
library(pheatmap)
ORN <- readRDS("./05_ORN_cluster2/01_first_cluster/ORN_integrated_antenna_first_cluster.rds")
ORN_last <- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak.rds")
ORN_filtered_cells<- setdiff(colnames(ORN),colnames(ORN_last))
library(tidyr)
barcode <- ORN_filtered_cells
barcode_sample<-separate(as.data.frame(barcode),"barcode",c("sample","barcode"),"_")
barcode <- barcode_sample$barcode
# random select 50 cell to get the bw file form filtered ORN 
random_50cell<- sample(barcode,50)
write.table(random_50cell,"/md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/5_RT_cluster/cell_barcode/random_50cell.txt",row.name=F,col.names=F)

# random select 100 cell to get the bw file form filtered ORN 
random_100cell<- sample(barcode,100)
write.table(random_100cell,"/md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/5_RT_cluster/cell_barcode/random_100cell.txt",row.name=F,col.names=F)
sed -i 's/"//g' *.txt

ls *.txt> ORN_barcode.list
for file1 in $(<ORN_barcode.list)
do
  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/merged.bam --cell-barcodes $file1 --cores 10 --out-bam $file1.bam
  bedtools  genomecov  -bg -split -ibam $file1.bam  > $file1.bedGraph
  echo "make normalized bedGraph"
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl $file1.bedGraph $file1.norm.bedGraph 
  sort -k1,1 -k2,2n $file1.norm.bedGraph > $file1.norm.sorted.bedGraph
  bedGraphToBigWig $file1.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt ./$file1.norm.bw 
done

# /md01/nieyg/project/honeybee/honebee-latest-Version/05_ORN_cluster2/05_combination_group_recluster/cell_barcode

## PBS configure 
#PBS -N Atacworks_denoise
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=6
#PBS -l mem=32G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/5_RT_cluster/cell_barcode
ls *.txt> ORN_barcode.list
for file1 in $(<ORN_barcode.list)
do
  /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/merged.bam --cell-barcodes $file1 --cores 10 --out-bam $file1.bam
  bedtools  genomecov  -bg -split -ibam $file1.bam  > $file1.bedGraph
  echo "make normalized bedGraph"
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl $file1.bedGraph $file1.norm.bedGraph 
  sort -k1,1 -k2,2n $file1.norm.bedGraph > $file1.norm.sorted.bedGraph
  bedGraphToBigWig $file1.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt ./$file1.norm.bw 
done

cd /md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/5_RT_cluster
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate python37

honebee_input=/md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/2_honeybee
atacworks=/md01/nieyg/software/AtacWorks
ATAC_bw=/md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/5_RT_cluster/cell_barcode

ls $ATAC_bw/*.bw |cut -d "/" -f 10 > ORN_ATAC_bw.list
for file1 in $(<ORN_ATAC_bw.list)
do
echo $file1
python $atacworks/scripts/main.py denoise \
    --noisybw $ATAC_bw/$file1 \
    --genome /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt \
    --weights_path $honebee_input/atacworks_train_OSN_filtered_latest/model_best.pth.tar \
    --out_home "./" \
    --regions C234.genome_intervals.bed \
    --exp_name "atacworks_denoise" \
    --distributed \
    --batch_size 10 \
    --num_workers 6 \
    --config /md01/nieyg/software/AtacWorks/configs/model_structure.yaml
done;





honebee_input=/md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/2_honeybee
atacworks=/md01/nieyg/software/AtacWorks
ATAC_bw=/md01/nieyg/project/honeybee/honebee-latest-Version/12_Atacworks/5_RT_cluster/cell_barcode

file1="C4.txt.norm.bw"

grep -v "GroupUN88" C234.genome_intervals.bed > temp.txt && mv temp.txt C234.genome_intervals.bed
python $atacworks/scripts/main.py denoise \
    --noisybw $ATAC_bw/$file1 \
    --genome /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt \
    --weights_path $honebee_input/atacworks_train_OSN_filtered_latest/model_best.pth.tar \
    --out_home "./" \
    --exp_name "atacworks_denoise" \
    --distributed \
    --batch_size 10 \
    --regions C234.genome_intervals.bed \
    --num_workers 6 \
    --config /md01/nieyg/software/AtacWorks/configs/model_structure.yaml




file1="random_50cell.txt.norm.bw"
grep -v "GroupUN11" C234.genome_intervals.bed > temp.txt && mv temp.txt C234.genome_intervals.bed
python $atacworks/scripts/main.py denoise \
    --noisybw $ATAC_bw/$file1 \
    --genome /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt \
    --weights_path $honebee_input/atacworks_train_OSN_filtered_latest/model_best.pth.tar \
    --out_home "./" \
    --exp_name "atacworks_denoise" \
    --distributed \
    --batch_size 10 \
    --regions C234.genome_intervals.bed \
    --num_workers 6 \
    --config /md01/nieyg/software/AtacWorks/configs/model_structure.yaml







