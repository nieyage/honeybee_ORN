# intergenic regions are drawn together and plotted with All Or
setwd("/data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/GC_output")
GC_table <- read.table("honeybee_Or_genebody_All.txt")
colnames(GC_table) <- c("sequence","GC")
GC_table$label <- "H genebody"
GC_table$AT <- 1-GC_table$GC

GC_table2 <- read.table("honeybee_All_Or_upstream_intergenic.txt")
colnames(GC_table2) <- c("sequence","GC")
GC_table2$label <- "H intergenic"
GC_table2$AT <- 1-GC_table2$GC

full_GC_table <- rbind(GC_table, GC_table2)

GC_table3 <- read.table("honeybee_All_Or_downstream_intergenic.txt")
colnames(GC_table3) <- c("sequence","GC")
GC_table3$label <- "H intergenic"
GC_table3$AT <- 1-GC_table3$GC

full_GC_table <- rbind(full_GC_table, GC_table3)

setwd("/data/R04/liyh526/project/Honeybee/15_Other_species_AT/AaegL5.0/GC_output")
GC_table10 <- read.table("AaegL_Or_genebody.txt")
colnames(GC_table10) <- c("sequence","GC")
GC_table10$label <- "A genebody"
GC_table10$AT <- 1-GC_table10$GC

full_GC_table <- rbind(full_GC_table, GC_table10)

GC_table11 <- read.table("AaegL_Or_upstream_intergenic.txt")
colnames(GC_table11) <- c("sequence","GC")
GC_table11$label <- "A intergenic"
GC_table11$AT <- 1-GC_table11$GC

full_GC_table <- rbind(full_GC_table, GC_table11)

GC_table12 <- read.table("AaegL_Or_downstream_intergenic.txt")
colnames(GC_table12) <- c("sequence","GC")
GC_table12$label <- "A intergenic"
GC_table12$AT <- 1-GC_table12$GC

full_GC_table <- rbind(full_GC_table, GC_table12)

setwd("/data/R04/liyh526/project/Honeybee/15_Other_species_AT/fly/GC_output")
GC_table7 <- read.table("Fly_Or_genebody.txt")
colnames(GC_table7) <- c("sequence","GC")
GC_table7$label <- "F genebody"
GC_table7$AT <- 1-GC_table7$GC

full_GC_table <- rbind(full_GC_table, GC_table7)

GC_table8 <- read.table("Fly_Or_upstream_intergenic.txt")
colnames(GC_table8) <- c("sequence","GC")
GC_table8$label <- "F intergenic"
GC_table8$AT <- 1-GC_table8$GC

full_GC_table <- rbind(full_GC_table, GC_table8)

GC_table9 <- read.table("Fly_Or_downstream_intergenic.txt")
colnames(GC_table9) <- c("sequence","GC")
GC_table9$label <- "F intergenic"
GC_table9$AT <- 1-GC_table9$GC

full_GC_table <- rbind(full_GC_table, GC_table9)

my_comparisons=list(c("H genebody","A genebody"),
                    c("H genebody","F genebody"),
                    c("H intergenic","A intergenic"),
                    c("H intergenic","F intergenic"))

pdf("/data/R04/liyh526/project/Honeybee/11_Fig/Three_species_AT_content_body_intergenic_T.test_among.pdf.pdf",width=10,height=6)
ggviolin(full_GC_table, x="label", y="AT", fill = "label",
         palette = c("#C3004D", "#C3004D", 
                     "#324B77", "#324B77", 
                     "#2F4858", "#2F4858"),
         add = "boxplot",
         add.params = list(fill="white")) +
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5)) +
  NoLegend() +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", 
                     label.y = c(0.96, 1, 1.04, 1.08))
dev.off()

# Group2 intergenic
setwd("/data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/GC_output")
GC_table <- read.table("honeybee_Or_genebody_Group2.txt")
colnames(GC_table) <- c("sequence","GC")
GC_table$label <- "H genebody"
GC_table$AT <- 1-GC_table$GC

GC_table2 <- read.table("honeybee_Group2_Or_upstream_intergenic.txt")
colnames(GC_table2) <- c("sequence","GC")
GC_table2$label <- "H intergenic"
GC_table2$AT <- 1-GC_table2$GC

full_GC_table <- rbind(GC_table, GC_table2)

GC_table3 <- read.table("honeybee_Group2_Or_downstream_intergenic.txt")
colnames(GC_table3) <- c("sequence","GC")
GC_table3$label <- "H intergenic"
GC_table3$AT <- 1-GC_table3$GC

full_GC_table <- rbind(full_GC_table, GC_table3)

my_comparisons=list(c("H genebody","H intergenic"))

pdf("/data/R04/liyh526/project/Honeybee/11_Fig/Honeybee_Group2_AT_content_body_intergenic_T.test.pdf.pdf",width=6,height=8)
ggviolin(full_GC_table, x="label", y="AT", fill = "label",
         palette = c("#C3004D", "#C3004D", "#C3004D"),
         add = "boxplot",
         add.params = list(fill="white")) +
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5)) +
  NoLegend() +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", 
                     label.y = c(0.92))
dev.off()

# The map was regrouped, and the upper and lower streams were separated, with three in each group
setwd("/data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/GC_output")
GC_table <- read.table("honeybee_Or_genebody_All.txt")
colnames(GC_table) <- c("sequence","GC")
GC_table$label <- "H genebody"
GC_table$AT <- 1-GC_table$GC

GC_table2 <- read.table("honeybee_All_Or_upstream_intergenic.txt")
colnames(GC_table2) <- c("sequence","GC")
GC_table2$label <- "H upstream intergenic"
GC_table2$AT <- 1-GC_table2$GC

full_GC_table <- rbind(GC_table, GC_table2)

GC_table3 <- read.table("honeybee_All_Or_downstream_intergenic.txt")
colnames(GC_table3) <- c("sequence","GC")
GC_table3$label <- "H downstream intergenic"
GC_table3$AT <- 1-GC_table3$GC

full_GC_table <- rbind(full_GC_table, GC_table3)

setwd("/data/R04/liyh526/project/Honeybee/15_Other_species_AT/AaegL5.0/GC_output")
GC_table10 <- read.table("AaegL_Or_genebody.txt")
colnames(GC_table10) <- c("sequence","GC")
GC_table10$label <- "A genebody"
GC_table10$AT <- 1-GC_table10$GC

full_GC_table <- rbind(full_GC_table, GC_table10)

GC_table11 <- read.table("AaegL_Or_upstream_intergenic.txt")
colnames(GC_table11) <- c("sequence","GC")
GC_table11$label <- "A upstream intergenic"
GC_table11$AT <- 1-GC_table11$GC

full_GC_table <- rbind(full_GC_table, GC_table11)

GC_table12 <- read.table("AaegL_Or_downstream_intergenic.txt")
colnames(GC_table12) <- c("sequence","GC")
GC_table12$label <- "A downstream intergenic"
GC_table12$AT <- 1-GC_table12$GC

full_GC_table <- rbind(full_GC_table, GC_table12)

setwd("/data/R04/liyh526/project/Honeybee/15_Other_species_AT/fly/GC_output")
GC_table7 <- read.table("Fly_Or_genebody.txt")
colnames(GC_table7) <- c("sequence","GC")
GC_table7$label <- "F genebody"
GC_table7$AT <- 1-GC_table7$GC

full_GC_table <- rbind(full_GC_table, GC_table7)

GC_table8 <- read.table("Fly_Or_upstream_intergenic.txt")
colnames(GC_table8) <- c("sequence","GC")
GC_table8$label <- "F upstream intergenic"
GC_table8$AT <- 1-GC_table8$GC

full_GC_table <- rbind(full_GC_table, GC_table8)

GC_table9 <- read.table("Fly_Or_downstream_intergenic.txt")
colnames(GC_table9) <- c("sequence","GC")
GC_table9$label <- "F downstream intergenic"
GC_table9$AT <- 1-GC_table9$GC

full_GC_table <- rbind(full_GC_table, GC_table9)

# upstream
dim(full_GC_table)
upstream_table <- full_GC_table[-which(full_GC_table$label == "H genebody"), ]
upstream_table <- upstream_table[-which(upstream_table$label == "F genebody"), ]
upstream_table <- upstream_table[-which(upstream_table$label == "A genebody"), ]
dim(upstream_table)
upstream_table <- upstream_table[-which(upstream_table$label == "H downstream intergenic"), ]
upstream_table <- upstream_table[-which(upstream_table$label == "F downstream intergenic"), ]
upstream_table <- upstream_table[-which(upstream_table$label == "A downstream intergenic"), ]
dim(upstream_table)

my_comparisons=list(c("H upstream intergenic","A upstream intergenic"),
                    c("H upstream intergenic","F upstream intergenic"),
                    c("A upstream intergenic","F upstream intergenic"))

pdf("/data/R04/liyh526/project/Honeybee/11_Fig/Three_species_AT_content_upstream_T.test.pdf",width=8,height=6)
ggviolin(upstream_table, x="label", y="AT", fill = "label",
         palette = c("#C3004D","#324B77", "#2F4858"),
         add = "boxplot",
         add.params = list(fill="white")) +
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5)) +
  NoLegend() +
  xlab("") +
  stat_compare_means(comparisons = my_comparisons,method = "t.test", 
                     label.y = c(0.97, 1.01, 1.06))
dev.off()

# downstream
dim(full_GC_table)
downstream_table <- full_GC_table[-which(full_GC_table$label == "H genebody"), ]
downstream_table <- downstream_table[-which(downstream_table$label == "F genebody"), ]
downstream_table <- downstream_table[-which(downstream_table$label == "A genebody"), ]
dim(downstream_table)
downstream_table <- downstream_table[-which(downstream_table$label == "H upstream intergenic"), ]
downstream_table <- downstream_table[-which(downstream_table$label == "F upstream intergenic"), ]
downstream_table <- downstream_table[-which(downstream_table$label == "A upstream intergenic"), ]
dim(downstream_table)

my_comparisons=list(c("H downstream intergenic","A downstream intergenic"),
                    c("H downstream intergenic","F downstream intergenic"),
                    c("A downstream intergenic","F downstream intergenic"))

pdf("/data/R04/liyh526/project/Honeybee/11_Fig/Three_species_AT_content_downstream_T.test.pdf",width=8,height=6)
ggviolin(downstream_table, x="label", y="AT", fill = "label",
         palette = c("#C3004D","#324B77", "#2F4858"),
         add = "boxplot",
         add.params = list(fill="white")) +
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5)) +
  NoLegend() +
  xlab("") +
  stat_compare_means(comparisons = my_comparisons,method = "t.test", 
                     label.y = c(0.99, 1.03, 1.08))
dev.off()


## Tandem repeat counting
# Read all the files I need
setwd("/data/R04/liyh526/project/Honeybee/15_Other_species_AT/repeat_region")
dir()
## Fly
simpleRepeat_UCSC_fly <- read.table("RepeatMasker_UCSC_fly.txt", header = F)

colnames(simpleRepeat_UCSC_fly) <- c("bin","swScore","milliDiv","milliDel","milliIns",
                                     "genoName","genoStart","genoEnd","genoLeft",
                                     "strand","repName","repClass","repFamily",
                                     "repStart","repEnd","repLeft","id")
# Remove unwanted columns
simpleRepeat_UCSC_fly <- simpleRepeat_UCSC_fly[,c(-1,-2,-3,-4,-5,-9,-13,-14,-15,-16,-17)]
# add Key
simpleRepeat_UCSC_fly$Key <- paste0("Fly",1:nrow(simpleRepeat_UCSC_fly))
# Save the bed file
simpleRepeat_UCSC_fly$genoName <- gsub("chr","",simpleRepeat_UCSC_fly$genoName)
# Remove unwanted areas
idx1 <- which(simpleRepeat_UCSC_fly$genoName == "2L")
length(idx1)
idx2 <- which(simpleRepeat_UCSC_fly$genoName == "2R")
length(idx2)
idx3 <- which(simpleRepeat_UCSC_fly$genoName == "3L")
length(idx3)
idx4 <- which(simpleRepeat_UCSC_fly$genoName == "3R")
length(idx4)
idx5 <- which(simpleRepeat_UCSC_fly$genoName == "4")
length(idx5)
idx <- c(idx1,idx2,idx3,idx4,idx5)
length(idx)

simpleRepeat_UCSC_fly <- simpleRepeat_UCSC_fly[idx,]
dim(simpleRepeat_UCSC_fly)

write.table(simpleRepeat_UCSC_fly[,c(1,2,3,7)], "simpleRepeat_UCSC_fly.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## Honeybee
simpleRepeat_UCSC_honeybee <- read.table("RepeatMasker_UCSC_honeybee.txt", header = F)

colnames(simpleRepeat_UCSC_honeybee) <- c('chrom','chromStart','chromEnd','name')
# add Key
simpleRepeat_UCSC_honeybee$Key <- paste0("Honeybee",1:nrow(simpleRepeat_UCSC_honeybee))
# Save the bed file
table(simpleRepeat_UCSC_honeybee$chrom)
simpleRepeat_UCSC_honeybee$chrom <- gsub("chrLG1","Group1",simpleRepeat_UCSC_honeybee$chrom)
simpleRepeat_UCSC_honeybee$chrom <- gsub("chrLG2","Group2",simpleRepeat_UCSC_honeybee$chrom)
simpleRepeat_UCSC_honeybee$chrom <- gsub("chrLG3","Group3",simpleRepeat_UCSC_honeybee$chrom)
simpleRepeat_UCSC_honeybee$chrom <- gsub("chrLG4","Group4",simpleRepeat_UCSC_honeybee$chrom)
simpleRepeat_UCSC_honeybee$chrom <- gsub("chrLG5","Group5",simpleRepeat_UCSC_honeybee$chrom)
simpleRepeat_UCSC_honeybee$chrom <- gsub("chrLG6","Group6",simpleRepeat_UCSC_honeybee$chrom)
simpleRepeat_UCSC_honeybee$chrom <- gsub("chrLG7","Group7",simpleRepeat_UCSC_honeybee$chrom)
simpleRepeat_UCSC_honeybee$chrom <- gsub("chrLG8","Group8",simpleRepeat_UCSC_honeybee$chrom)
simpleRepeat_UCSC_honeybee$chrom <- gsub("chrLG9","Group9",simpleRepeat_UCSC_honeybee$chrom)
simpleRepeat_UCSC_honeybee$chrom <- gsub("chrLG10","Group10",simpleRepeat_UCSC_honeybee$chrom)
simpleRepeat_UCSC_honeybee$chrom <- gsub("chrLG11","Group11",simpleRepeat_UCSC_honeybee$chrom)
simpleRepeat_UCSC_honeybee$chrom <- gsub("chrLG12","Group12",simpleRepeat_UCSC_honeybee$chrom)
simpleRepeat_UCSC_honeybee$chrom <- gsub("chrLG13","Group13",simpleRepeat_UCSC_honeybee$chrom)
simpleRepeat_UCSC_honeybee$chrom <- gsub("chrLG14","Group14",simpleRepeat_UCSC_honeybee$chrom)
simpleRepeat_UCSC_honeybee$chrom <- gsub("chrLG15","Group15",simpleRepeat_UCSC_honeybee$chrom)
simpleRepeat_UCSC_honeybee$chrom <- gsub("chrLG16","Group16",simpleRepeat_UCSC_honeybee$chrom)
table(simpleRepeat_UCSC_honeybee$chrom)
idx1 <- which(simpleRepeat_UCSC_honeybee$chrom == "Group1")
length(idx1)
idx2 <- which(simpleRepeat_UCSC_honeybee$chrom == "Group2")
length(idx2)
idx3 <- which(simpleRepeat_UCSC_honeybee$chrom == "Group3")
length(idx3)

idx4 <- which(simpleRepeat_UCSC_honeybee$chrom == "Group4")
length(idx4)
idx5 <- which(simpleRepeat_UCSC_honeybee$chrom == "Group5")
length(idx5)
idx6 <- which(simpleRepeat_UCSC_honeybee$chrom == "Group6")
length(idx6)

idx7 <- which(simpleRepeat_UCSC_honeybee$chrom == "Group7")
length(idx7)
idx8 <- which(simpleRepeat_UCSC_honeybee$chrom == "Group8")
length(idx8)
idx9 <- which(simpleRepeat_UCSC_honeybee$chrom == "Group9")
length(idx9)

idx10 <- which(simpleRepeat_UCSC_honeybee$chrom == "Group10")
length(idx10)
idx11 <- which(simpleRepeat_UCSC_honeybee$chrom == "Group11")
length(idx11)
idx12 <- which(simpleRepeat_UCSC_honeybee$chrom == "Group12")
length(idx12)

idx13 <- which(simpleRepeat_UCSC_honeybee$chrom == "Group13")
length(idx13)
idx14 <- which(simpleRepeat_UCSC_honeybee$chrom == "Group14")
length(idx14)
idx15 <- which(simpleRepeat_UCSC_honeybee$chrom == "Group15")
length(idx15)
idx16 <- which(simpleRepeat_UCSC_honeybee$chrom == "Group16")
length(idx16)

idx <- c(idx1,idx2,idx3,idx4,idx5,idx6,
         idx7,idx8,idx9,idx10,idx11,idx12,
         idx13,idx14,idx15,idx16)
length(idx)
simpleRepeat_UCSC_honeybee <- simpleRepeat_UCSC_honeybee[idx,]
dim(simpleRepeat_UCSC_honeybee)
write.table(simpleRepeat_UCSC_honeybee[,c(1,2,3,5)], "simpleRepeat_UCSC_honeybee.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## mosquito
simpleRepeat_UCSC_mosquito <- read.table("RepeatMasker_UCSC_mosquito.txt", header = F)
colnames(simpleRepeat_UCSC_mosquito) <- c('chrom','chromStart','chromEnd','name')
# add Key
simpleRepeat_UCSC_mosquito$Key <- paste0("mosquito",1:nrow(simpleRepeat_UCSC_mosquito))
# Save the bed file
table(simpleRepeat_UCSC_mosquito$chrom)
simpleRepeat_UCSC_mosquito$chrom <- gsub("chr1","NC_035107.1",simpleRepeat_UCSC_mosquito$chrom)
simpleRepeat_UCSC_mosquito$chrom <- gsub("chr2","NC_035108.1",simpleRepeat_UCSC_mosquito$chrom)
simpleRepeat_UCSC_mosquito$chrom <- gsub("chr3","NC_035109.1",simpleRepeat_UCSC_mosquito$chrom)
# Remove unwanted areas
idx1 <- which(simpleRepeat_UCSC_mosquito$chrom == "NC_035107.1")
length(idx1)
idx2 <- which(simpleRepeat_UCSC_mosquito$chrom == "NC_035108.1")
length(idx2)
idx3 <- which(simpleRepeat_UCSC_mosquito$chrom == "NC_035109.1")
length(idx3)
idx <- c(idx1,idx2,idx3)
length(idx)
simpleRepeat_UCSC_mosquito <- simpleRepeat_UCSC_mosquito[idx,]
dim(simpleRepeat_UCSC_mosquito)
write.table(simpleRepeat_UCSC_mosquito[,c(1,2,3,5)], "simpleRepeat_UCSC_mosquito.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

##### linux #####
# Calculate AT ratio
# honeybee
bedtools getfasta -fi /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/genome.fa \
-bed /data/R04/liyh526/project/Honeybee/15_Other_species_AT/repeat_region/simpleRepeat_UCSC_honeybee.bed \
-fo /data/R04/liyh526/project/Honeybee/15_Other_species_AT/repeat_region/fasta/simpleRepeat_UCSC_honeybee.fasta

cd /data/R04/liyh526/project/Honeybee/15_Other_species_AT/repeat_region/fasta
geecee -sequence simpleRepeat_UCSC_honeybee.fasta -outfile ../GC_output/simpleRepeat_UCSC_honeybee.txt
# mosquito
bedtools getfasta -fi /data/R04/liyh526/project/Honeybee/15_Other_species_AT/AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.fna \
-bed /data/R04/liyh526/project/Honeybee/15_Other_species_AT/repeat_region/simpleRepeat_UCSC_mosquito.bed \
-fo /data/R04/liyh526/project/Honeybee/15_Other_species_AT/repeat_region/fasta/simpleRepeat_UCSC_mosquito.fasta

cd /data/R04/liyh526/project/Honeybee/15_Other_species_AT/repeat_region/fasta
geecee -sequence simpleRepeat_UCSC_mosquito.fasta -outfile ../GC_output/simpleRepeat_UCSC_mosquito.txt
# fly
bedtools getfasta -fi /data/R04/liyh526/project/Honeybee/15_Other_species_AT/fly/remain.fasta \
-bed /data/R04/liyh526/project/Honeybee/15_Other_species_AT/repeat_region/simpleRepeat_UCSC_fly.bed \
-fo /data/R04/liyh526/project/Honeybee/15_Other_species_AT/repeat_region/fasta/simpleRepeat_UCSC_fly.fasta

cd /data/R04/liyh526/project/Honeybee/15_Other_species_AT/repeat_region/fasta
geecee -sequence simpleRepeat_UCSC_fly.fasta -outfile ../GC_output/simpleRepeat_UCSC_fly.txt

# We read in the GC content that we got
setwd("/data/R04/liyh526/project/Honeybee/15_Other_species_AT/repeat_region/GC_output")
dir()
fly_GC <- read.table("simpleRepeat_UCSC_fly.txt")
simpleRepeat_UCSC_fly$GC <- fly_GC$V2 
honeybee_GC <- read.table("simpleRepeat_UCSC_honeybee.txt")
simpleRepeat_UCSC_honeybee$GC <- honeybee_GC$V2
mosquito_GC <- read.table("simpleRepeat_UCSC_mosquito.txt")
simpleRepeat_UCSC_mosquito$GC <- mosquito_GC$V2

## Read the intergenic regions
setwd("/data/R04/liyh526/project/Honeybee/15_Other_species_AT/AaegL5.0")
AaegL_Or_downstream_intergenic <- read.table("AaegL_Or_downstream_intergenic.bed", header = F)
AaegL_Or_upstream_intergenic <- read.table("AaegL_Or_upstream_intergenic.bed", header = F)
setwd("/data/R04/liyh526/project/Honeybee/15_Other_species_AT/fly")
Fly_Or_downstream_intergenic <- read.table("Fly_Or_downstream_intergenic.bed", header = F)
Fly_Or_upstream_intergenic <- read.table("Fly_Or_upstream_intergenic.bed", header = F)
setwd("/data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee")
honeybee_All_Or_downstream_intergenic <- read.table("honeybee_All_Or_downstream_intergenic.bed", header = F)
honeybee_All_Or_upstream_intergenic <- read.table("honeybee_All_Or_upstream_intergenic.bed", header = F)

# AT Tandem repeat counting
# Fly down
# Initialize the counter
count <- 0

# Line by line comparison
for (i in 1:nrow(simpleRepeat_UCSC_fly)) {
  for (j in 1:nrow(Fly_Or_downstream_intergenic)) {
    if (simpleRepeat_UCSC_fly[i, "genoName"] == Fly_Or_downstream_intergenic[j, "V1"] & 
        simpleRepeat_UCSC_fly[i, "genoStart"] >= Fly_Or_downstream_intergenic[j, "V2"] & 
        simpleRepeat_UCSC_fly[i, "genoEnd"] <= Fly_Or_downstream_intergenic[j, "V3"] &
        simpleRepeat_UCSC_fly[i, "GC"] <= 0.2) {
      count <- count + 1
    }
  }
}

print(count)
Fly_downstream_ATrepeat_count <- count

# Fly up
# Initialize the counter
count <- 0

# Line by line comparison
for (i in 1:nrow(simpleRepeat_UCSC_fly)) {
  for (j in 1:nrow(Fly_Or_upstream_intergenic)) {
    if (simpleRepeat_UCSC_fly[i, "genoName"] == Fly_Or_upstream_intergenic[j, "V1"] & 
        simpleRepeat_UCSC_fly[i, "genoStart"] >= Fly_Or_upstream_intergenic[j, "V2"] & 
        simpleRepeat_UCSC_fly[i, "genoEnd"] <= Fly_Or_upstream_intergenic[j, "V3"] &
        simpleRepeat_UCSC_fly[i, "GC"] <= 0.2) {
      count <- count + 1
    }
  }
}

print(count)
Fly_upstream_ATrepeat_count <- count


# Honeybee down
# Initialize the counter
count <- 0

# Line by line comparison
for (i in 1:nrow(simpleRepeat_UCSC_honeybee)) {
  for (j in 1:nrow(honeybee_All_Or_downstream_intergenic)) {
    if (simpleRepeat_UCSC_honeybee[i, "chrom"] == honeybee_All_Or_downstream_intergenic[j, "V1"] & 
        simpleRepeat_UCSC_honeybee[i, "chromStart"] >= honeybee_All_Or_downstream_intergenic[j, "V2"] & 
        simpleRepeat_UCSC_honeybee[i, "chromEnd"] <= honeybee_All_Or_downstream_intergenic[j, "V3"] &
        simpleRepeat_UCSC_honeybee[i, "GC"] <= 0.2) {
      count <- count + 1
    }
  }
}

print(count)
Honeybee_downstream_ATrepeat_count <- count

# Honeybee up
# Initialize the counter
count <- 0

# Line by line comparison
for (i in 1:nrow(simpleRepeat_UCSC_honeybee)) {
  for (j in 1:nrow(honeybee_All_Or_upstream_intergenic)) {
    if (simpleRepeat_UCSC_honeybee[i, "chrom"] == honeybee_All_Or_upstream_intergenic[j, "V1"] & 
        simpleRepeat_UCSC_honeybee[i, "chromStart"] >= honeybee_All_Or_upstream_intergenic[j, "V2"] & 
        simpleRepeat_UCSC_honeybee[i, "chromEnd"] <= honeybee_All_Or_upstream_intergenic[j, "V3"] &
        simpleRepeat_UCSC_honeybee[i, "GC"] <= 0.2) {
      count <- count + 1
    }
  }
}

print(count)
Honeybee_upstream_ATrepeat_count <- count


# Mosquito down
# Initialize the counter
count <- 0

# Line by line comparison
for (i in 1:nrow(simpleRepeat_UCSC_mosquito)) {
  for (j in 1:nrow(AaegL_Or_downstream_intergenic)) {
    if (simpleRepeat_UCSC_mosquito[i, "chrom"] == AaegL_Or_downstream_intergenic[j, "V1"] & 
        simpleRepeat_UCSC_mosquito[i, "chromStart"] >= AaegL_Or_downstream_intergenic[j, "V2"] & 
        simpleRepeat_UCSC_mosquito[i, "chromEnd"] <= AaegL_Or_downstream_intergenic[j, "V3"] &
        simpleRepeat_UCSC_mosquito[i, "GC"] <= 0.2) {
      count <- count + 1
    }
  }
}

print(count)
Mosquito_downstream_ATrepeat_count <- count

# Mosquito up
# Initialize the counter
count <- 0

# Line by line comparison
for (i in 1:nrow(simpleRepeat_UCSC_mosquito)) {
  for (j in 1:nrow(AaegL_Or_upstream_intergenic)) {
    if (simpleRepeat_UCSC_mosquito[i, "chrom"] == AaegL_Or_upstream_intergenic[j, "V1"] & 
        simpleRepeat_UCSC_mosquito[i, "chromStart"] >= AaegL_Or_upstream_intergenic[j, "V2"] & 
        simpleRepeat_UCSC_mosquito[i, "chromEnd"] <= AaegL_Or_upstream_intergenic[j, "V3"] &
        simpleRepeat_UCSC_mosquito[i, "GC"] <= 0.2) {
      count <- count + 1
    }
  }
}

print(count)
Mosquito_upstream_ATrepeat_count <- count

results <- data.frame(
  Type = c("H downstream", "H upstream",
           "A downstream", "A upstream",
           "F downstream", "F upstream"),
  Count = c(Honeybee_downstream_ATrepeat_count, Honeybee_upstream_ATrepeat_count,
            Mosquito_downstream_ATrepeat_count, Mosquito_upstream_ATrepeat_count,
            Fly_downstream_ATrepeat_count, Fly_upstream_ATrepeat_count)
)

results$Type <- factor(results$Type,levels = c("H downstream", "H upstream",
                                               "A downstream", "A upstream",
                                               "F downstream", "F upstream"))

results_down <- results[c(1,3,5),]
results_up <- results[c(2,4,6),]

pdf("/data/R04/liyh526/project/Honeybee/11_Fig/Three_species_AT_tandem_RepeatMasker.pdf",width=6,height=6)
ggplot(results, aes(x = Type, y = Count, fill = Type)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("H downstream" = "#C3004D", 
                               "H upstream" = "#C3004D",
                               "A downstream" = "#324B77", 
                               "A upstream" = "#324B77",
                               "F downstream" = "#2F4858",
                               "F upstream" = "#2F4858")) +
  labs(title = "Intergenic AT tandem Repeat", y = "Count") + 
  theme_classic() +
  NoLegend() +
  xlab("")
dev.off()




### Intergenic_region_AT_percent
# AT percent
# Fly
# Upstream and downstream separation
bedtools getfasta -fi /data/R04/liyh526/project/Honeybee/15_Other_species_AT/fly/remain.fasta \
                  -bed /data/R04/liyh526/project/Honeybee/15_Other_species_AT/fly/Fly_Or_upstream_intergenic.bed \
                  -fo /data/R04/liyh526/project/Honeybee/15_Other_species_AT/fly/fasta/Fly_Or_upstream_intergenic.fasta

cd /data/R04/liyh526/project/Honeybee/15_Other_species_AT/fly/fasta
geecee -sequence Fly_Or_upstream_intergenic.fasta -outfile ../GC_output/Fly_Or_upstream_intergenic.txt

bedtools getfasta -fi /data/R04/liyh526/project/Honeybee/15_Other_species_AT/fly/remain.fasta \
                  -bed /data/R04/liyh526/project/Honeybee/15_Other_species_AT/fly/Fly_Or_downstream_intergenic.bed \
                  -fo /data/R04/liyh526/project/Honeybee/15_Other_species_AT/fly/fasta/Fly_Or_downstream_intergenic.fasta

cd /data/R04/liyh526/project/Honeybee/15_Other_species_AT/fly/fasta
geecee -sequence Fly_Or_downstream_intergenic.fasta -outfile ../GC_output/Fly_Or_downstream_intergenic.txt

bedtools getfasta -fi /data/R04/liyh526/project/Honeybee/15_Other_species_AT/fly/remain.fasta \
                  -bed /data/R04/liyh526/project/Honeybee/15_Other_species_AT/fly/Fly_Or_genebody.bed \
                  -fo /data/R04/liyh526/project/Honeybee/15_Other_species_AT/fly/fasta/Fly_Or_genebody.fasta

cd /data/R04/liyh526/project/Honeybee/15_Other_species_AT/fly/fasta
geecee -sequence Fly_Or_genebody.fasta -outfile ../GC_output/Fly_Or_genebody.txt

# mosquito
# Upstream and downstream separation
bedtools getfasta -fi /data/R04/liyh526/project/Honeybee/15_Other_species_AT/AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.fna \
                  -bed /data/R04/liyh526/project/Honeybee/15_Other_species_AT/AaegL5.0/AaegL_Or_upstream_intergenic.bed \
                  -fo /data/R04/liyh526/project/Honeybee/15_Other_species_AT/AaegL5.0/fasta/AaegL_Or_upstream_intergenic.fasta

cd /data/R04/liyh526/project/Honeybee/15_Other_species_AT/AaegL5.0/fasta
geecee -sequence AaegL_Or_upstream_intergenic.fasta -outfile ../GC_output/AaegL_Or_upstream_intergenic.txt

bedtools getfasta -fi /data/R04/liyh526/project/Honeybee/15_Other_species_AT/AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.fna \
                  -bed /data/R04/liyh526/project/Honeybee/15_Other_species_AT/AaegL5.0/AaegL_Or_downstream_intergenic.bed \
                  -fo /data/R04/liyh526/project/Honeybee/15_Other_species_AT/AaegL5.0/fasta/AaegL_Or_downstream_intergenic.fasta

cd /data/R04/liyh526/project/Honeybee/15_Other_species_AT/AaegL5.0/fasta
geecee -sequence AaegL_Or_downstream_intergenic.fasta -outfile ../GC_output/AaegL_Or_downstream_intergenic.txt

bedtools getfasta -fi /data/R04/liyh526/project/Honeybee/15_Other_species_AT/AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.fna \
                  -bed /data/R04/liyh526/project/Honeybee/15_Other_species_AT/AaegL5.0/AaegL_Or_genebody.bed \
                  -fo /data/R04/liyh526/project/Honeybee/15_Other_species_AT/AaegL5.0/fasta/AaegL_Or_genebody.fasta

cd /data/R04/liyh526/project/Honeybee/15_Other_species_AT/AaegL5.0/fasta
geecee -sequence AaegL_Or_genebody.fasta -outfile ../GC_output/AaegL_Or_genebody.txt

# honeybee
# Upstream and downstream separation
# Group2
bedtools getfasta -fi /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/genome.fa \
                  -bed /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/honeybee_Group2_Or_upstream_intergenic.bed \
                  -fo /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/fasta/honeybee_Group2_Or_upstream_intergenic.fasta

cd /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/fasta
geecee -sequence honeybee_Group2_Or_upstream_intergenic.fasta -outfile ../GC_output/honeybee_Group2_Or_upstream_intergenic.txt

bedtools getfasta -fi /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/genome.fa \
                  -bed /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/honeybee_Group2_Or_downstream_intergenic.bed \
                  -fo /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/fasta/honeybee_Group2_Or_downstream_intergenic.fasta

cd /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/fasta
geecee -sequence honeybee_Group2_Or_downstream_intergenic.fasta -outfile ../GC_output/honeybee_Group2_Or_downstream_intergenic.txt

bedtools getfasta -fi /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/genome.fa \
                  -bed /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/honeybee_Or_genebody_Group2.bed \
                  -fo /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/fasta/honeybee_Or_genebody_Group2.fasta

cd /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/fasta
geecee -sequence honeybee_Or_genebody_Group2.fasta -outfile ../GC_output/honeybee_Or_genebody_Group2.txt

# All
bedtools getfasta -fi /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/genome.fa \
                  -bed /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/honeybee_All_Or_upstream_intergenic.bed \
                  -fo /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/fasta/honeybee_All_Or_upstream_intergenic.fasta

cd /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/fasta
geecee -sequence honeybee_All_Or_upstream_intergenic.fasta -outfile ../GC_output/honeybee_All_Or_upstream_intergenic.txt

bedtools getfasta -fi /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/genome.fa \
                  -bed /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/honeybee_All_Or_downstream_intergenic.bed \
                  -fo /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/fasta/honeybee_All_Or_downstream_intergenic.fasta

cd /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/fasta
geecee -sequence honeybee_All_Or_downstream_intergenic.fasta -outfile ../GC_output/honeybee_All_Or_downstream_intergenic.txt

bedtools getfasta -fi /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/genome.fa \
                  -bed /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/honeybee_Or_genebody_All.bed \
                  -fo /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/fasta/honeybee_Or_genebody_All.fasta

cd /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/fasta
geecee -sequence honeybee_Or_genebody_All.fasta -outfile ../GC_output/honeybee_Or_genebody_All.txt


# TSS AT
bedtools getfasta -fi /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/genome.fa \
                  -bed /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/honeybee_All_Or_TSS_UP_DOWN_1k.bed \
                  -fo /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/fasta/honeybee_All_Or_TSS_UP_DOWN_1k.fasta

cd /data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/fasta
geecee -sequence honeybee_All_Or_TSS_UP_DOWN_1k.fasta -outfile ../GC_output/honeybee_All_Or_TSS_UP_DOWN_1k.txt
