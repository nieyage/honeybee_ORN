# TOBIAS_TFBScan
# motif.scan overlap FIND

# The single-cell peaks merge the peaks of bulk ATAC
# motif scan of all Or promoter sequences of honeybees
cd /data/R04/liyh526/project/Honeybee/05_honeybee_bulkATAC_macs2_callpeak/2023.5.15/bulk_sc_associate/Or_promoter_fasta
for i in *.fasta
do
echo "$i"
TOBIAS TFBScan --motifs /data/R04/liyh526/project/Honeybee/03_CisBP2JASPAR/honeybee_jaspar/CisBP-honeybee.jaspar \
               --fasta "$i" --cores 8 --pvalue 0.0005 \
               --outfile /data/R04/liyh526/project/Honeybee/05_honeybee_bulkATAC_macs2_callpeak/2023.5.15/bulk_sc_associate/TOBIAS_motif/Or-0.0005/${i%.*}_motif.bed
done

#!/bin/bash

# Gets all the bed files in the folder
# Loop over the bed files for comparison

bedfiles=(/data/R04/liyh526/project/Honeybee/05_honeybee_bulkATAC_macs2_callpeak/2023.5.15/bulk_sc_associate/TOBIAS_motif/Or_pair_motif_wehave_0.0005/*.bed)
resultfolder=/data/R04/liyh526/project/Honeybee/05_honeybee_bulkATAC_macs2_callpeak/2023.5.15/bulk_sc_associate/TOBIAS_motif/Or_motif_scan_overlap_0.0005/

for ((i=0; i<${#bedfiles[@]}-1; i++))
do
  for ((j=i+1; j<${#bedfiles[@]}; j++))
  do
    # Output the motif common to the two Ors
    awk 'FNR==NR{a[$4];next}$4 in a' "${bedfiles[$i]}" "${bedfiles[$j]}" > "${resultfolder}/${bedfiles[$i]##*/}_${bedfiles[$j]##*/}.txt"
  done
done

cd /data/R04/liyh526/project/Honeybee/05_honeybee_bulkATAC_macs2_callpeak/2023.5.15/bulk_sc_associate/TOBIAS_motif/Or_motif_scan_overlap_0.0005
# Batch change file names
for file in *; do
    if [[ "$file" =~ Or_([[:alnum:]-]+)_promoter_motif\.bed_Or_([[:alnum:]-]+)_promoter_motif\.bed\.txt ]]; then
        newname="${BASH_REMATCH[1]}_${BASH_REMATCH[2]}_promoter_motif.bed"
        cp "$file" "$newname"
    else
        echo "No rename for $file"
    fi
done

# Take a look at the motif that we scan in the shell with TOBIAS
library(dplyr)
library(readr)
library(tidyr)

# Creates an empty result data frame
result_df <- data.frame()

# Process each.bed file one by one
setwd("/data/R04/liyh526/project/Honeybee/16_TF/TF_candidate/Or-0.0005")
bed_files <- dir()

for (bed_file in bed_files) {
  # Read the.bed file and create the data frame
  df <- read_delim(bed_file, delim = "\t", col_names = FALSE)
  colnames(df) <- c("Region", "Column2", "Column3", "MotifID", "Bindscore", "Column6")
  df <- df[, c("MotifID", "Bindscore")]
  
  file_name <- sub(".*_(.*?)_.*", "\\1", bed_file)
  
  # For the same MotifID, the Bindscore are summed
  df <- df %>%
    group_by(MotifID) %>%
    summarise(Bindscore = sum(Bindscore))
  
  df$symbol <- file_name
  
  result_df <- bind_rows(result_df, df)
}

head(result_df)
dim(result_df)

for (i in 1:length(result_df$MotifID)) {
  result_df$GBID[i] <- unlist(strsplit(result_df$MotifID[i], "_"))[1]
}

# MotifID
new_rows <- data.frame()

# Traverse the raw data frame of each line
for (i in 1:nrow(result_df)) {
  # Gets information about the current row
  row_info <- result_df[i, ]
  
  # Split GBID with " " as the separator
  ids <- unique(unlist(str_match_all(row_info$GBID, "(GB[0-9]+)")))
  
  # Create a new line
  for (id in ids) {
    new_row <- row_info
    new_row$GBID <- id
    new_rows <- rbind(new_rows, new_row)
  }
}

dim(new_rows)
result_df_annotated <- merge(new_rows, TF_symbol, by.x="GBID", by.y ="Gene.Symbol.and.Synonyms")

head(result_df_annotated)
result_df_annotated <- result_df_annotated[,c(3,4,5)]
dim(result_df_annotated)
result_df_annotated <- result_df_annotated[order(result_df_annotated$symbol),]
dim(result_df_annotated)

result_df_annotated_dedup <- result_df_annotated %>%
  group_by(symbol, Bindscore) %>%
  distinct(Gene.Symbol, .keep_all = TRUE)

# Using dcast function converts long matrix to matrix wide
matrix_df <- reshape2::dcast(result_df_annotated_dedup,
                             symbol ~ Gene.Symbol,
                             drop = T,
                             value.var = "Bindscore")
head(matrix_df)
write.csv(matrix_df,file = "matrix_df.csv")
# In this way, the motif obtained by scan in each Or cluster and its corresponding bindscore are obtained

setwd("/data/R04/liyh526/project/Honeybee/16_TF/TF_candidate/Foldchange")
ORN <- readRDS("/data/R02/nieyg/project/honeybee/honebee-latest-Version/05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_latest.rds")
Idents(ORN)
for (i in levels(ORN)) {
  # get foldchange
  Idents(ORN)
  # First, create a new vector clustering cluster labels, will all clustering cluster is initialized to "Others"
  new_clusters <- rep("Others", length(length(ORN$celltype)))
  # Store the cluster names that need to be retained in a variable
  cluster_to_keep <- i
  # Find the index of the cluster that needs to be retained
  cluster_to_keep_index <- which(ORN$celltype == cluster_to_keep)
  # Will need to keep the clustering cluster label set to its original label
  new_clusters[cluster_to_keep_index] <- cluster_to_keep
  # Sets the other clustering cluster labels to "Others"
  other_clusters_index <- which(!(ORN$celltype %in% cluster_to_keep))
  new_clusters[other_clusters_index] <- "Others"
  ORN$MP_cluster <- new_clusters
  Idents(ORN)<- ORN$MP_cluster
  DefaultAssay(ORN) <- "raw_RNA"
  markers_df <- FindMarkers(object = ORN, 
                            ident.1 = i, 
                            ident.2 = "Others", 
                            logfc.threshold = 0, 
                            test.use = "wilcox", 
                            min.pct = 0, # filter features
                            only.pos = FALSE)
  write.csv(markers_df,file = paste0(i,"_","Foldchange",".csv"),
            row.names =T, col.names =F, quote =F, sep = "\t")
  Idents(ORN) <- ORN$celltype
  print(i)
}

# Define the cell types
levels(ORN)
new.cluster.ids <- c("1"="LOC410603,LOC107963999,Or63-b", 
                     "2"="LOC410603,LOC107963999",
                     "3"="LOC100578045,LOC107963999",
                     "4"="LOC100578045", 
                     "5"="LOC413596",
                     "6"="LOC102653615,LOC102653695", 
                     "7"="Or130", 
                     "8"="LOC102656221,LOC102656904", 
                     "9"="LOC102653782,LOC102656904", 
                     "10"="LOC100577068,LOC100577101,Or41", 
                     "11"="LOC100576881",
                     "12"="LOC102656907,Or55,Or56", 
                     "13"="Or56",
                     "14"="LOC102653637,LOC102653703,Or57,Or58", 
                     "15"="LOC102653637,LOC102653703,LOC102655435",
                     "16"="LOC102653703,LOC102655435",
                     "17"="LOC725205,LOC100577787", 
                     "18"="Or4,Or5",
                     "19"="LOC100577715,Or12",
                     "20"="LOC100577671,LOC100577715", 
                     "21"="LOC100577590", 
                     "22"="Or26", 
                     "23"="LOC725861",
                     "24"="LOC726097,Or25,Or26,Or27,Or30-a", 
                     "25"="Or30-b,Or35",
                     "26"="LOC724763",
                     "27"="LOC100576522",
                     "28"="LOC107965760",
                     "29"="LOC102654841,Or170",
                     "30"="LOC100576153",
                     "31"="LOC107966050",
                     "32"="LOC102655559",
                     "33"="LOC100578724,LOC100579011-b",
                     "34"="LOC100578842,LOC102656567-a,LOC102656567-b,LOC107965011",
                     "35"="LOC107965011",
                     "36"="LOC102655285,LOC107965761",
                     "37"="LOC102679224",
                     "38"="LOC724673,LOC102655218,LOC102679224",
                     "39"="LOC100578751-b",
                     "40"="LOC100578617",
                     "41"="Or115",
                     "42"="LOC100577068",
                     "43"="LOC100577334",
                     "44"="LOC724911,LOC725052",
                     "45"="LOC725052,LOC100576839",
                     "46"="LOC100576839",
                     "47"="LOC100577955",
                     "48"="LOC102655147,Or105",
                     "49"="LOC102655147,LOC102655180,LOC102656805,Or105",
                     "50"="LOC102655180,LOC102656805,Or98,Or105",
                     "51"="LOC100576212",
                     "52"="LOC100576282",
                     "53"="LOC100576246",
                     "54"="LOC102653716",
                     "55"="Or109",
                     "56"="Or107",
                     "57"="LOC100577888-a,LOC102653865",
                     "58"="LOC100577888-a",
                     "59"="LOC102653810"
)
ORN <- RenameIdents(ORN, new.cluster.ids)
ORN$celltype <- ORN@active.ident

#setwd("/data/R04/liyh526/project/Honeybee/16_TF/TF_candidate/bubbleplot")
for (i in levels(ORN)) {
  levels(ORN)
  new_clusters <- rep("Others", length(length(ORN$celltype)))
  cluster_to_keep <- i
  cluster_to_keep_index <- which(ORN$celltype == cluster_to_keep)
  new_clusters[cluster_to_keep_index] <- cluster_to_keep
  other_clusters_index <- which(!(ORN$celltype %in% cluster_to_keep))
  new_clusters[other_clusters_index] <- "Others"
  ORN$MP_cluster <- new_clusters
  Idents(ORN)<- ORN$MP_cluster
  elements <- strsplit(i, ",")[[1]]
  for (j in elements) {
    if (j %in% unique(result_df_annotated_dedup$symbol) == T) {
      TF <- unique(result_df_annotated_dedup[result_df_annotated_dedup$symbol == j,]$Gene.Symbol)
      DefaultAssay(ORN) <- "raw_RNA"
      Avg_matrix <- AverageExpression(ORN,features=TF,assays = "raw_RNA",slot = "counts")
      Avg_matrix <- as.data.frame(Avg_matrix)
      for (k in 1:length(colnames(Avg_matrix))) {
        colnames(Avg_matrix)[k] <- gsub("[.]", ",", colnames(Avg_matrix)[k])
        colnames(Avg_matrix)[k] <- gsub("raw_RNA,", "", colnames(Avg_matrix)[k])
      }
      Avg_matrix_rowsum <- cbind(Avg_matrix,rowsum=rowSums(Avg_matrix))
      dim(Avg_matrix_rowsum)
      # Add a list of TAU
      Avg_matrix_rowsum <- cbind(Avg_matrix_rowsum,TAU=rowSums(Avg_matrix))
      # Calculate the TAU
      for (l in 1:nrow(Avg_matrix)) {
        Avg_matrix_rowsum[l,4] <- 2/1 - Avg_matrix_rowsum[l,3]/(1*max(Avg_matrix[l,]))
      }
      # Transcription factors with zero expression were removed
      dim(Avg_matrix_rowsum)
      Avg_matrix_rowsum_TAU_rmNaN <- Avg_matrix_rowsum[!is.nan(Avg_matrix_rowsum[,4]),]
      dim(Avg_matrix_rowsum_TAU_rmNaN)
      colnames(Avg_matrix_rowsum_TAU_rmNaN)[2] <- "exp"
      exp_tau <- Avg_matrix_rowsum_TAU_rmNaN[,c("exp","TAU")]
      promoter_motif <- result_df_annotated_dedup[result_df_annotated_dedup$symbol == j,]
      promoter_motif <- promoter_motif[!duplicated(promoter_motif$Gene.Symbol),]
      dim(promoter_motif)
      head(promoter_motif)
      exp_tau$symbol <- rownames(exp_tau)
      exp_tau <- merge(exp_tau,promoter_motif,by.x="symbol", by.y="Gene.Symbol")
      head(exp_tau)
      dim(exp_tau)
      markers_df <- read.csv(paste0(i,"_","Foldchange",".csv"))
      head(markers_df)
      colnames(markers_df)[1] <- "symbol"
      exp_tau <- merge(exp_tau,markers_df,by.x="symbol", by.y="symbol")
      
      write.csv(exp_tau,file = paste0(j,"_","bubbleplot",".csv"))
      print(j)
    } else {
      next
    }
  }
  Idents(ORN) <- ORN$celltype
  print(i)
}

# Draw a bubble plot
library(ggplot2)
library(ggrepel)
library(pheatmap)
setwd("/data/R04/liyh526/project/Honeybee/16_TF/TF_candidate/bubbleplot")
bubbleplot_data <- dir(pattern = "*_bubbleplot.csv", full.names = F)

fileNames <- gsub("_bubbleplot.csv","",bubbleplot_data)

for (i in fileNames) {
  data <- read.csv(paste0(i,"_","bubbleplot",".csv"))
  data <- data[,-1]
  max(data$exp)
  min(data$exp)
  median(data$exp)
  # if the range(data$exp) too large: data$exp<- log2(data$exp+1); (min(data$exp)==0,log2时+1，如果min(data$exp)> 0,log2(data$exp))
  max(data$Bindscore)
  min(data$Bindscore)
  median(data$Bindscore)
  # Filter meet the threshold condition of genes
  head(data)
  data <- data[which(data$exp > 0),]
  threshold_x <- median(data$Bindscore) # x阈值
  threshold_y <- 0.75   # y阈值
  filtered_data <- data[data$Bindscore > threshold_x & data$TAU > threshold_y, ]
  # Genes with 0 expression were removed
  filtered_data <- filtered_data[which(filtered_data$exp > 0),]
  dim(filtered_data)
  head(filtered_data)
  write.csv(filtered_data,file = paste0(i,"_","filtered_candidate",".csv"))
  # Take a look at the expression of the TF candidates we screened
  # cluster
  TF <- unique(filtered_data$symbol)
  Avg <- AverageExpression(ORN,features=TF,assays = "raw_RNA", slot = "counts")
  count=t(scale(t(Avg$raw_RNA),scale = T,center = F))
  pheatmap <- pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),show_rownames=T,show_colnames=T)
  pdf(paste0("./result/",i,"_filtered_TF_candidate_exp_pheatmap.pdf"),width=20,height=10)
  print(pheatmap)
  dev.off()
  
  # Bubble Plot
  p <- ggplot(data, aes(x = Bindscore, y = TAU, size = exp, fill = avg_log2FC)) +
    geom_point(shape = 21, stroke = 0, alpha = 2) +
    geom_text_repel(data = filtered_data, aes(label = symbol), color = "black", box.padding = 0.5, point.padding = 0, size = 3) +
    scale_size_identity() +
    scale_size(range = c(5, 15), name="Exp") +
    scale_fill_viridis(discrete=FALSE, 
                       option="F", 
                       direction = -1,
                       end = 0.9,
                       begin = 0.3) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    xlim(0,max(data$Bindscore)) +
    ylim(0.25,1) +
    labs(fill = "Log2FC") +
    xlab("BindScore") +
    ylab("Tau")
  
  # Add threshold line
  p <- p + geom_vline(xintercept = threshold_x, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = threshold_y, linetype = "dashed", color = "gray")
  
  pdf(paste0("./result/Bindscore_TAU_",i,"_bubbleplot.pdf"),width=12,height=10)
  print(p)
  dev.off()  
  print(i)
}

data.plot_V2 <- read.csv("/data/R04/liyh526/project/Honeybee/16_TF/TF_candidate/bubbleplot/Honeybee_Or_TF_Dotplot_V2.csv")
data.plot_V2 <- data.plot_V2[,c(-6,-7,-8,-9,-10,-11)]
data.plot_V2$id <- factor(data.plot_V2$id, levels=unique(data.plot_V2$id))
data.plot_V2 <- data.plot_V2[-which(data.plot_V2$id == "Or26"), ]
data.plot_V2 <- data.plot_V2[-which(data.plot_V2$id == "Or27"), ]

library(ggplot2)
p1<-ggplot(data.plot_V2,aes(x=id,
                            y=features.plot))+
  geom_point(aes(size=`Bindscore`,
                 color=`Average.Expression`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="lightgrey",high="blue")+
  labs(x = 'TF', y = 'Or')+
  guides(size=guide_legend(order=3))+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))

data.plot_V2$state<-"No";
for(i in 1:nrow(data.plot_V2)){
  if(data.plot_V2[i,]$Bindscore>4&&data.plot_V2[i,]$Average.Expression >= 0.1){data.plot_V2[i,]$state="Yes"};
}
data.plot_V3 <- data.plot_V2[which(data.plot_V2$state=="Yes"),]
unique(data.plot_V3$features.plot)

# 创建示例向量A和向量B
vector_A <- unique(data.plot_V2$features.plot)  # 示例向量A
vector_B <- unique(data.plot_V3$features.plot)  # 已排序的部分

# 使用order()函数对向量A的未排序部分进行排序，并获取排序后的索引
unsorted_indices <- setdiff(seq_along(vector_A), match(vector_B, vector_A))
sorted_unsorted_part <- vector_A[unsorted_indices[order(vector_A[unsorted_indices])]]

# 使用向量B的排序和未排序部分的排序合并得到最终结果
final_sorted_vector <- c(vector_B, sorted_unsorted_part)
final_sorted_vector <- rev(final_sorted_vector)

data.plot_V2$features.plot <- factor(data.plot_V2$features.plot, levels=final_sorted_vector)

data.plot_V3 <- data.plot_V2[-which(data.plot_V2$features.plot %in% sorted_unsorted_part),]
unique(data.plot_V3$features.plot)

data.plot_V3$id <- factor(data.plot_V3$id, levels=c("Or63-b","Or132","Or133N","Or130","Or102","Or128","Or129F","Or40","Or41",
                                                    "Or108","Or109","Or107","Or55","Or59","Or57","Or58","Or1","Or3",
                                                    "Or12","Or106","Or20","Or22","Or31","Or25","Or30-a","Or30-b","Or35","Or160",
                                                    "Or13a","Or164F","Or169","Or152","Or151","Or150","Or146","Or154",
                                                    "Or163","Or167PAR","Or165C","Or115","Or33","Or91","Or94","Or96",    
                                                    "Or161","Or101","Or98","Or83F","Or79","Or39"))

p1<-ggplot(data.plot_V3,aes(x=id,
                            y=features.plot))+
  geom_point(aes(size=`Bindscore`,
                 color=`Average.Expression`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="lightgrey",high="blue")+
  labs(x = 'TF', y = 'Or')+
  guides(size=guide_legend(order=3))+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))

pdf("Fig5_Honeybee_Or_TF_Dotplot.pdf",width=20,height=16)
p1
dev.off()

# Single example
# 128 & 129 Or pair
# Co-expression of multiple promoters / Motif scan
ORN <- readRDS("/data/R02/nieyg/project/honeybee/honebee-latest-Version/05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_latest.rds")
Idents(ORN)
new_clusters <- rep("Others", length(length(ORN$cell_group)))

cluster_to_keep <- "8"

cluster_to_keep_index <- which(ORN$cell_group == cluster_to_keep)

new_clusters[cluster_to_keep_index] <- cluster_to_keep

other_clusters_index <- which(!(ORN$cell_group %in% cluster_to_keep))
new_clusters[other_clusters_index] <- "Others"

ORN$MP_cluster <- new_clusters
Idents(ORN) <- ORN$MP_cluster

DefaultAssay(ORN) <- "raw_RNA"

markers_df <- FindMarkers(object = ORN,
                          ident.1 = "8",
                          ident.2 = "Others",
                          logfc.threshold = 0,
                          test.use = "wilcox",
                          min.pct = 0, # filter features
                          only.pos = FALSE)

# get foldchange and other related information
gene_names <- rownames(markers_df)
fold_changes <- markers_df$avg_log2FC

# OR128_OR129(LOC102656904_LOC102656221)
OR128_OR129_promoter_motif <- read.table("/data/R04/liyh526/project/Honeybee/05_honeybee_bulkATAC_macs2_callpeak/2023.5.15/bulk_sc_associate/TOBIAS_motif/Or_motif_scan_overlap_0.001/LOC102656221_LOC102656904_promoter_motif.bed")
for (i in 1:length(OR128_OR129_promoter_motif$V4)) {
  OR128_OR129_promoter_motif$V7[i] <- unlist(strsplit(OR128_OR129_promoter_motif$V4[i], "_"))[1]
}

OR128_OR129_promoter_motif <- merge(OR128_OR129_promoter_motif, TF_symbol, by.x="V7", by.y ="Gene.Symbol.and.Synonyms")
OR128_OR129_promoter_motif <- OR128_OR129_promoter_motif[,c(-2,-3,-4,-6,-7)]
dim(OR128_OR129_promoter_motif)
unique(OR128_OR129_promoter_motif$Gene.Symbol)

TF <- unique(OR128_OR129_promoter_motif$Gene.Symbol)
DefaultAssay(ORN) <- "raw_RNA"
Avg_OR128_OR129 <- AverageExpression(ORN,features=TF,assays = "raw_RNA",slot = "counts")

Avg_OR128_OR129 <- as.data.frame(Avg_OR128_OR129)
dim(Avg_OR128_OR129)
for (i in 1:length(colnames(Avg_OR128_OR129))) {
  colnames(Avg_OR128_OR129)[i] <- gsub("raw_RNA.", "", colnames(Avg_OR128_OR129)[i])
  print(i)
}

class(Avg_OR128_OR129)
dim(Avg_OR128_OR129)
Avg_OR128_OR129_rowsum <- cbind(Avg_OR128_OR129,rowsum=rowSums(Avg_OR128_OR129))
dim(Avg_OR128_OR129_rowsum)
Avg_OR128_OR129_rowsum <- cbind(Avg_OR128_OR129_rowsum,TAU=rowSums(Avg_OR128_OR129))
ncol(Avg_OR128_OR129_rowsum)
head(Avg_OR128_OR129_rowsum)

# calculation TAU
for (i in 1:nrow(Avg_OR128_OR129)) {
  Avg_OR128_OR129_rowsum[i,4] <- 2/1 - Avg_OR128_OR129_rowsum[i,3]/(1*max(Avg_OR128_OR129[i,]))
  print(i)
}

# rm Transcription factors with zero expression
dim(Avg_OR128_OR129_rowsum)
Avg_OR128_OR129_rowsum_TAU_rmNaN <- Avg_OR128_OR129_rowsum[!is.nan(Avg_OR128_OR129_rowsum[,4]),]
dim(Avg_OR128_OR129_rowsum_TAU_rmNaN)

exp_tau_OR128_OR129 <- Avg_OR128_OR129_rowsum_TAU_rmNaN[,c("8","TAU")]
OR128_OR129_promoter_motif <- OR128_OR129_promoter_motif[!duplicated(OR128_OR129_promoter_motif),]
dim(OR128_OR129_promoter_motif)
head(OR128_OR129_promoter_motif)
exp_tau_OR128_OR129$symbol <- rownames(exp_tau_OR128_OR129)
exp_tau_OR128_OR129 <- merge(exp_tau_OR128_OR129,OR128_OR129_promoter_motif,by.x="symbol", by.y="Gene.Symbol")
colnames(exp_tau_OR128_OR129) <- c("symbol","exp","tau","ID","Motifid")
head(exp_tau_OR128_OR129)
dim(exp_tau_OR128_OR129)

head(markers_df)
markers_df$symbol <- rownames(markers_df)
exp_tau_OR128_OR129 <- merge(exp_tau_OR128_OR129,markers_df,by.x="symbol", by.y="symbol")

setwd("/data/R04/liyh526/project/Honeybee/16_TF")
write.csv(exp_tau_OR128_OR129,"./exp_tau_OR128_OR129.csv")

# calculated the average binding score
csv_file="/data/R04/liyh526/project/Honeybee/16_TF/exp_tau_OR128_OR129.csv"

bed_file1="/data/R04/liyh526/project/Honeybee/16_TF/Or_LOC102656221_promoter_motif.bed"
bed_file2="/data/R04/liyh526/project/Honeybee/16_TF/Or_LOC102656904_promoter_motif.bed"

data <- read.csv(csv_file)

data$Bed1_Sum <- NA
data$Bed2_Sum <- NA

for (i in 1:nrow(data)) {
  motifid <- data[i, "Motifid"]
  
  bed1_data <- read.table(bed_file1)
  sum1 <- sum(bed1_data$V5[bed1_data$V4 == motifid])
  
  bed2_data <- read.table(bed_file2)
  sum2 <- sum(bed2_data$V5[bed2_data$V4 == motifid])
  
  data[i, "Bed1_Sum"] <- sum1
  data[i, "Bed2_Sum"] <- sum2
}

data$average_bindscore <- (data$Bed1_Sum + data$Bed2_Sum)/2
data <- data[,-1]
write.csv(data, "./exp_tau_bindscore_OR128_OR129.csv", row.names = FALSE)

setwd("/data/R04/liyh526/project/Honeybee/16_TF")
data <- read.csv("exp_tau_bindscore_OR128_OR129.csv")
# bubble plot
library(ggplot2)
library(ggrepel)

max(data$exp)
min(data$exp)
median(data$exp)
# if the range(data$exp) too large: data$exp<- log2(data$exp+1); (min(data$exp)==0,log2时+1，如果min(data$exp)> 0,log2(data$exp))

max(data$average_bindscore)
min(data$average_bindscore)
median(data$average_bindscore)
# Genes that met the threshold conditions were screened
threshold_x <- 2  # x阈值
threshold_y <- 0.8   # y阈值
filtered_data <- data[data$average_bindscore > threshold_x & data$tau > threshold_y, ]
# rm Transcription factors with zero expression
filtered_data <- filtered_data[which(filtered_data$exp > 0),]
dim(filtered_data)
head(filtered_data)
# rm Transcription factors with zero expression
data <- data[which(data$exp > 0),]
unique(data[data$average_bindscore > -10 & data$tau > 0.25, ]$symbol)

TF<- unique(filtered_data$symbol)
Avg <- AverageExpression(ORN,features=TF,assays = "raw_RNA", slot = "counts")
library(pheatmap)
count=t(scale(t(Avg$raw_RNA),scale = T,center = F))
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),show_rownames=T,show_colnames=T)
VlnPlot(ORN, TF, ncol = 3,pt.size = 0,assay = "raw_RNA", slot = "counts")
DefaultAssay(ORN) <- "raw_RNA"
TF <- TF[c(2,3)]
DimPlot(ORN)
pdf("/data/R04/liyh526/project/Honeybee/11_Fig/FeaturePlot_LOC411133_LOC552079_C8.PDF",width=6,height=6)
FeaturePlot(ORN, TF[1], ncol = 1, max.cutoff = 3)
FeaturePlot(ORN, TF[2], ncol = 1, max.cutoff = 30)
dev.off()
library("viridis")

colors_for_exp_pattern<- c("grey","#D51B07")
Vln <- VlnPlot(ORN,"LOC411133",ncol = 1,pt.size = 0,assay = "raw_RNA", slot = "counts") +
  scale_color_manual(values =colors_for_exp_pattern)
Vln_df <- Vln$data

pdf("/data/R04/liyh526/project/Honeybee/11_Fig/FigS4D_MP_ViolinPlot_LOC552079_C8.PDF",width=4,height=6)
VlnPlot(ORN, TF[2], cols = c("grey","#D51B07"), ncol = 1,pt.size = 0,assay = "raw_RNA", slot = "counts") + NoLegend()
dev.off()

pdf("/data/R04/liyh526/project/Honeybee/11_Fig/FigS4E_MP_ViolinPlot_LOC411133_C8.PDF",width=4,height=6)
VlnPlot(ORN, TF[1], cols = c("grey","#18798B"), ncol = 1,pt.size = 0,assay = "raw_RNA", slot = "counts") + NoLegend()
dev.off()

pdf("/data/R04/liyh526/project/Honeybee/11_Fig/FigS4D_MP_boxplot_LOC552079_C8.PDF",width=4,height=6)
ggboxplot(Vln_df, x="ident", y="LOC552079",color = "ident",width=0.6, notch = F)+
  stat_compare_means(label.x = 1.35, label.y = max(Vln_df$LOC411133)+1)+theme(legend.position="none") +
  xlab("") +
  ylab("TF Expression") +
  scale_color_manual(values = c("grey","#18798B")) +
  ggtitle("LOC552079") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("/data/R04/liyh526/project/Honeybee/11_Fig/FigS4E_MP_boxplot_LOC411133_C8.PDF",width=4,height=6)
ggboxplot(Vln_df, x="ident", y="LOC411133",color = "ident",width=0.6, notch = F)+
  stat_compare_means(label.x = 1.35, label.y = max(Vln_df$LOC411133)+1)+theme(legend.position="none") +
  xlab("") +
  ylab("TF Expression") +
  scale_color_manual(values = c("grey","#18798B")) +
  ggtitle("LOC411133") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# 绘制Bubble Plot
p <- ggplot(data, aes(x = average_bindscore, y = tau, size = exp, fill = avg_log2FC)) +
  geom_point(shape = 21, stroke = 0, alpha = 2) +
  geom_text_repel(data = filtered_data, aes(label = symbol), color = "black", box.padding = 0.5, point.padding = 0, size = 3) +
  scale_size_identity() +
  scale_size(range = c(5, 15), name="Exp") +
  scale_fill_viridis(discrete=FALSE, 
                     option="F", 
                     direction = -1,
                     end = 0.9,
                     begin = 0.3) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  xlim(-10,20) +
  ylim(0.25,1) +
  labs(fill = "Log2FC") +
  xlab("BindScore") +
  ylab("Tau")

# Add a threshold line
p <- p + geom_vline(xintercept = threshold_x, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = threshold_y, linetype = "dashed", color = "gray")

pdf("/data/R04/liyh526/project/Honeybee/11_Fig/FigS4C_exp_tau_bindscore_OR128_OR129_bubbleplot.pdf",width=12,height=10)
print(p)
dev.off()

DefaultAssay(ORN) <- "peaks_ORN_subcluster"
# motid plot
peak_plot <- PeakPlot(
  object = ORN,
  region = "Group4-6013575-6019462"
)
peak_df <- peak_plot$data
peak.gr <- GRanges(peak_df$seqnames,IRanges(start=peak_df$start,end=peak_df$end),strand = "*")
# read in gr of the four transcription factors
setwd("/data/R04/liyh526/project/Honeybee/16_TF")
dir()
LOC100577980.gr <- readPeakFile("LOC100577980_motifs.bed", as = "GRanges")
LOC411133.gr <- readPeakFile("LOC411133_motifs.bed", as = "GRanges")
LOC552079.gr <- readPeakFile("LOC552079_motifs.bed", as = "GRanges")
LOC100577980_pos_in_peak.gr <- LOC100577980.gr[unique(queryHits(findOverlaps(LOC100577980.gr,peak.gr))),]
LOC411133_pos_in_peak.gr <- LOC411133.gr[unique(queryHits(findOverlaps(LOC411133.gr,peak.gr))),]
LOC552079_pos_in_peak.gr <- LOC552079.gr[unique(queryHits(findOverlaps(LOC552079.gr,peak.gr))),]
LOC100577980_pos_plot <- ggplot(data=as.data.frame(LOC100577980_pos_in_peak.gr)) +
  geom_segment(aes(x=(start + end) / 2, y=ifelse(V6 == "+", 0.1, -0.1) ,xend = (start + end) / 2, yend = 0),size = 3,color="#059749")+
  theme_classic() +
  ylab(label = "LOC100577980") +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),plot.margin=unit(c(0,0,0,0), "cm"),panel.background=element_rect(fill='transparent', color='black',linetype="solid")) +
  xlim(c(6012975, 6019962)) +
  geom_rect(aes(xmin=6015801,xmax=6016000,ymin=-Inf,ymax=Inf),
            fill="gray",alpha=0.2) +
  geom_rect(aes(xmin=6019083,xmax=6019282,ymin=-Inf,ymax=Inf),
            fill="gray",alpha=0.2)
LOC411133_pos_plot <- ggplot(data=as.data.frame(LOC411133_pos_in_peak.gr)) +
  geom_segment(aes(x=(start + end) / 2, y=ifelse(V6 == "+", 0.1, -0.1) ,xend = (start + end) / 2, yend = 0),size = 4,color="#18798B")+
  theme_classic() +
  ylab(label = "LOC411133") +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),plot.margin=unit(c(0,0,0,0), "cm"),panel.background=element_rect(fill='transparent', color='black',linetype="solid")) +
  xlim(c(6012975, 6019962)) +
  geom_rect(aes(xmin=6015801,xmax=6016000,ymin=-Inf,ymax=Inf),
            fill="gray",alpha=0.2) +
  geom_rect(aes(xmin=6019083,xmax=6019282,ymin=-Inf,ymax=Inf),
            fill="gray",alpha=0.2)
LOC552079_pos_plot <- ggplot(data=as.data.frame(LOC552079_pos_in_peak.gr)) +
  geom_segment(aes(x=(start + end) / 2, y=ifelse(V6 == "+", 0.1, -0.1) ,xend = (start + end) / 2, yend = 0),size = 4,color="#D51B07")+
  theme_classic() +
  ylab(label = "LOC552079") +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),plot.margin=unit(c(0,0,0,0), "cm"),panel.background=element_rect(fill='transparent', color='black',linetype="solid")) +
  xlim(c(6012975, 6019962)) +
  geom_rect(aes(xmin=6015801,xmax=6016000,ymin=-Inf,ymax=Inf),
            fill="gray",alpha=0.2) +
  geom_rect(aes(xmin=6019083,xmax=6019282,ymin=-Inf,ymax=Inf),
            fill="gray",alpha=0.2)

gene_plot <- AnnotationPlot(object = ORN, 
                            region = "Group4-6013575-6019462", 
                            extend.upstream = 600,
                            extend.downstream = 600)

pdf("/data/R04/liyh526/project/Honeybee/11_Fig/Fig4_TFBS_C8_Or128_Or129_motifplot.pdf",height=8,width=16)
LOC411133_pos_plot / LOC100577980_pos_plot / LOC552079_pos_plot / gene_plot + plot_layout(heights=c(0.15,0.15,0.15,0.15))
dev.off()

