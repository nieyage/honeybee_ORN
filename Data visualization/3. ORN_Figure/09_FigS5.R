library(Signac)
library(Seurat)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
library(ggplot2)
library(dplyr)
set.seed(1234)
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(c("LOC412949","LOC100577496","LOC102653640","LOC727346","LOC100578352","LOC552552","LOC726019","LOC551704","LOC410623","LOC100576097","LOC409777"))
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
all_receptor_gene <- unique(c(Orco,OR_gene,IR_gene,GR_gene))
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_latest.rds")
DefaultAssay(ORN)<-"raw_RNA"

tree_anno<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS5/OR_gene_tree_anno2.csv")
TSS_AT<- read.table("/data/R04/liyh526/project/Honeybee/15_Other_species_AT/honeybee/GC_output/honeybee_All_Or_TSS_ATcontent_UP_DOWN_1k.bed")
OR_target<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS5/OR_iORdb_pubChem.csv")

# merge the 3 dataframe 
tree_anno$Pherobase<- OR_target$Pherobase[match(tree_anno$OR_gene,OR_target$OR_gene)]
tree_anno$Name<- OR_target$Name[match(tree_anno$OR_gene,OR_target$OR_gene)]
tree_anno$Pubchem<- OR_target$Pubchem[match(tree_anno$OR_gene,OR_target$OR_gene)]
tree_anno$origin<- OR_target$origin[match(tree_anno$OR_gene,OR_target$OR_gene)]

TSS_AT$V4<- gsub("_TSS","",TSS_AT$V4)
tree_anno$TSS_AT_percent<- TSS_AT$V5[match(tree_anno$OR_gene,TSS_AT$V4)]

# calculate the TES length

gtf <- rtracklayer::import('/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz')
gtf$gene_biotype<-"protein_coding"
gene.coords <- gtf[gtf$gene_biotype == 'protein_coding']
gene.coords <- gene.coords[gene.coords$gene_name %in% all_receptor_gene]
gene.coords <- gene.coords[gene.coords$type == "transcript"]
gene.coords <- as.data.frame(gene.coords)
strand_p <- gene.coords[gene.coords$strand == "+",]
strand_n <- gene.coords[gene.coords$strand == "-",]

# for strand "+" OR gene 
strand_p <- strand_p[!duplicated(strand_p$gene_name),]
strand_p <- strand_p[order(strand_p$seqnames,strand_p$start,decreasing=F),]
OR_TES_location<- data.frame()
for (i in 1:nrow(strand_p)){
	OR1 <- strand_p[i,]
	if(i!=nrow(strand_p)){
	OR2 <- strand_p[i+1,]
	if(OR1$seqnames==OR2$seqnames){
		OR1_TES_start <- OR1$end+1;
		intersections <- abs(OR1$end-OR2$start);
		OR1_TES_end <- OR2$start-1;
		OR1_TES_location <- data.frame(OR=OR1$gene_name,seqnames=OR1$seqnames,start=OR1_TES_start,end=OR1_TES_end,strand=OR1$strand)
		OR_TES_location <- rbind(OR_TES_location,OR1_TES_location)
	}
}
}
OR_TES_location_strand_p<- OR_TES_location;

# for strand "-" OR gene 
strand_n <- strand_n[!duplicated(strand_n$gene_name),]
strand_n <- strand_n[order(strand_n$seqnames,strand_n$start,decreasing=T),]

OR_TES_location<- data.frame()
for (i in 1:nrow(strand_n)){
	OR1 <- strand_n[i,]
	if(i!=nrow(strand_n)){
	OR2 <- strand_n[i+1,]
	if(OR1$seqnames==OR2$seqnames){
		OR1_TES_start <- OR1$start-1;
		intersections <- abs(OR1$start-OR2$end);
		OR1_TES_end <- OR2$end+1;
		OR1_TES_location <- data.frame(OR=OR1$gene_name,seqnames=OR1$seqnames,start=OR1_TES_end,end=OR1_TES_start,strand=OR1$strand)
		OR_TES_location <- rbind(OR_TES_location,OR1_TES_location)
	}
}
}
OR_TES_location$strand<- "-"
OR_TES_location<- rbind(OR_TES_location,OR_TES_location_strand_p)
OR_TES_location$seqnames<- as.character(OR_TES_location$seqnames)
OR_TES_location$TES_length<- abs(OR_TES_location$start- OR_TES_location$end) 
tree_anno$TES_length<- OR_TES_location$TES_length[match(tree_anno$OR_gene,OR_TES_location$OR)]
#tree_anno$TES_length[is.na(tree_anno$TES_length)]<-"undefine"
tree_anno$origin[is.na(tree_anno$origin)]<-"undefine"
tree_anno$Pubchem[is.na(tree_anno$Pubchem)]<-"undefine"
tree_anno$Name[is.na(tree_anno$Name)]<-"undefine"
tree_anno$Pherobase[is.na(tree_anno$Pherobase)]<-"undefine"
tree_anno$Pherobase[tree_anno$Pherobase=="-"]<-"undefine"

write.csv(tree_anno,"./00_Figure/FigS5/FigS5_OR_gene_tree_annotation_info.csv")

# step1: plot the OR ggtree 
# sequence tree 
library(Biostrings)
library(muscle)
library(ape);
library(ggtree)
library(tidytree)
OR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/OR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
GR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/GR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
IR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/IR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
supply_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/supply.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
#OR2 is placed in the last column;
all_OR_gene_fasta<- c(OR_fasta,supply_fasta)

all_OR_gene_fasta<- all_OR_gene_fasta[which(names(all_OR_gene_fasta)%in% tree_anno$OR_gene),]
names(all_OR_gene_fasta) <- tree_anno$last_name[match(names(all_OR_gene_fasta),tree_anno$OR_gene)]

aln <- muscle::muscle(all_OR_gene_fasta)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming")
clust <- hclust(sdist,method="complete")#"ward.D"’, ‘"ward.D2"’,‘"single"’, ‘"complete"’, ‘"average"’ (= UPGMA), ‘"mcquitty"’, ‘"median"’ or ‘"centroid"’ (= UPGMC)
tree <- as.phylo(clust)

library(ggtreeExtra)
library(ggplot2)
library(ggtree)
library(treeio)
library(ggnewscale)
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                    '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', 
                    '#E0D4CA', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', 
                    '#6778AE', '#B53E2B', '#DCC1DD', '#CCE0F5', '#625D9E', 
                    '#68A180', '#968175', '#FF9999', '#344CB7', '#FFCC1D', 
                    '#24A19C', '#FF9A76',"#BC80BD", "#CCEBC5", "#FFED6F", 
                    "#E95C59", "#476D87",'#9FA3A8')


tree_anno<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS5/OR_gene_tree_anno2.csv")
Group =    c("Group1"=myUmapcolors[1],    "Group2"=myUmapcolors[2],    "Group3"=myUmapcolors[3],    "Group4"=myUmapcolors[4],    "Group5"=myUmapcolors[5],    "Group6"=myUmapcolors[6],    "Group7"=myUmapcolors[7],    "Group8"=myUmapcolors[8],    "Group9"=myUmapcolors[9],    "Group10"=myUmapcolors[10],
    "Group11"=myUmapcolors[11],    "Group12"=myUmapcolors[12],    "Group13"=myUmapcolors[13],    "Group14"=myUmapcolors[14],    "Group15"=myUmapcolors[15],    "Group16"=myUmapcolors[16],    "GroupUN3"=myUmapcolors[17],
    "GroupUN226"=myUmapcolors[18],    "GroupUN243"=myUmapcolors[19],    "GroupUN248"=myUmapcolors[20])
rownames(tree_anno)<- tree_anno$last_name
data<-  as.data.frame(tree_anno[,c(5,6,9,8)])
data_long<- melt(data, id.vars = c("last_name"), #需保留的不参与聚合的变量列名
                  measure.vars = c('seqnames','subfamily','exp_pattern2'),
                  variable.name = c('POS'),#聚合变量的新列名
                  value.name = 'value')
order<- c("Group1","Group2" ,"Group4","Group5","Group7","Group9","Group10",
"Group11","Group12","Group13","Group14","Group15","Group16","GroupUN243","GroupUN248",
unique(data$subfamily)[2:19],
unique(data$exp_pattern)[c(1,2)],
"undefine"
)
data_long$value<- factor(data_long$value,levels=order)
pdf("./00_Figure/FigS5/FigS5-OR_sequence_protein_similarity-tree_add_anno.pdf",width=15,height=15)
p1<- ggtree(tree, layout="fan",size=0.5) + geom_tiplab(size=3)
p1
p2<- p1+ new_scale_fill() + 
      geom_fruit(
          data=data_long,
          geom=geom_tile,
          mapping=aes(y=last_name, x=POS,fill=value),
          offset=0.2,   # The distance between external layers, default is 0.03 times of x range of tree.
          pwidth=0.2 # width of the external layer, default is 0.2 times of x range of tree.
      ) +
      scale_fill_manual(
          values=c(myUmapcolors),
          guide=guide_legend(keywidth=1, keyheight=1, order=3)
      ) 
p2
##data2<-  tree_anno[,c(5,17)]
#p3<-p2+
#      geom_fruit(
#          data = data2,
#          geom = geom_col,
#          mapping = aes(y=last_name, x=TSS_AT_percent,),
#      ) + 
#      new_scale_fill()
#data3<-  tree_anno[,c(5,18)]
#p4<-p3+
#      geom_fruit(
#          data = data3,
#          geom = geom_col,
#          mapping = aes(y=last_name, x=TES_length),
#      ) + 
#      new_scale_fill()
#
#p4
dev.off()

cluster_info<-as.data.frame(table(dotplot_data$id))
tree_anno<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS5/OR_gene_tree_anno2.csv")

dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower_latest.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));
cluster_info<-as.data.frame(table(dotplot_data$id))
multiOR_cluster<-as.character(cluster_info[cluster_info$Freq>1,1])
singleOR_cluster<-as.character(cluster_info[cluster_info$Freq==1,1])
multiOR<- unique(dotplot_data[which(dotplot_data$id %in% multiOR_cluster),]$features.plot)
singleOR<- unique(dotplot_data[which(dotplot_data$id %in% singleOR_cluster),]$features.plot)


multiOR_obj<- subset(ORN,idents=multiOR_cluster)
singleOR_obj<- subset(ORN,idents=singleOR_cluster)


order<- c("Group1","Group2" ,"Group4","Group5","Group7","Group9","Group10",
"Group11","Group12","Group13","Group14","Group15","Group16","GroupUN243","GroupUN248",
unique(data$subfamily)[2:19],
unique(data$exp_pattern)[c(1,2,4,5)],
unique(data$Pherobase)[2:9],
unique(data$origin)[2:4],"undefine"
)

# seqnames
tree_anno<- tree_anno[tree_anno$exp_pattern2!="undefine",]
tree_anno$Pherobase<- OR_target$Pherobase[match(tree_anno$OR_gene,OR_target$OR_gene)]
tree_anno$origin<- OR_target$origin[match(tree_anno$OR_gene,OR_target$OR_gene)]
tree_anno$origin[is.na(tree_anno$origin)]<-"undefine"
tree_anno$Pherobase[is.na(tree_anno$Pherobase)]<-"undefine"
tree_anno$Pherobase[tree_anno$Pherobase=="-"]<-"undefine"




pdf("./00_Figure/FigS5/FigS5-coexp-gene_proportion.pdf",width=10,height=4)
data<- as.data.frame(table(tree_anno$seqnames,tree_anno$exp_pattern2))
colnames(data)<-c("seqnames","exp_pattern","Freq")
data$seqnames<- factor(data$seqnames,levels=order)
p<-ggplot(data = data, aes_string(x = "seqnames", y = "Freq", 
        fill = "exp_pattern")) +  xlab(" ") + ylab(" ") + 
        scale_fill_manual(values = c("#B781CC","#7A9BCC","#CC6E7B")) + 
        geom_bar( stat = "identity", width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));
p
data<- as.data.frame(table(tree_anno$subfamily,tree_anno$exp_pattern2))
colnames(data)<-c("subfamily","exp_pattern","Freq")
#data$subfamily<- factor(data$subfamily,levels=order)
p<-ggplot(data = data, aes_string(x = "subfamily", y = "Freq", 
        fill = "exp_pattern")) +  xlab(" ") + ylab(" ") + 
        scale_fill_manual(values = c("#B781CC","#7A9BCC","#CC6E7B")) + 
        geom_bar( stat = "identity", width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));
p
data<- as.data.frame(table(tree_anno$Pherobase,tree_anno$exp_pattern2))
colnames(data)<-c("Pherobase","exp_pattern","Freq")
p<-ggplot(data = data, aes_string(x = "Pherobase", y = "Freq", 
        fill = "exp_pattern")) +  xlab(" ") + ylab(" ") + 
        scale_fill_manual(values = c("#B781CC","#7A9BCC","#CC6E7B")) + 
        geom_bar( stat = "identity", width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));
p
data<- as.data.frame(table(tree_anno$origin,tree_anno$exp_pattern2))
colnames(data)<-c("origin","exp_pattern","Freq")
p<-ggplot(data = data, aes_string(x = "origin", y = "Freq", 
        fill = "exp_pattern")) +  xlab(" ") + ylab(" ") + 
        scale_fill_manual(values = c("#B781CC","#7A9BCC","#CC6E7B")) + 
        geom_bar( stat = "identity", width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));
p
dev.off();

# the OR count 
coexp_OR<- tree_anno[tree_anno$exp_pattern2=="co-exp",]$OR_gene
single_OR<- tree_anno[tree_anno$exp_pattern2=="single_OR",]$OR_gene

all_OR<- c(coexp_OR,single_OR)
OR_data<- t(as.matrix(ORN@assays$RNA[all_OR,]))
obj_data<-colSums(as.data.frame(t(as.matrix(ORN@assays$RNA[all_OR,]))))
OR_data[OR_data>1]=1

dotplot_data<- data.frame(OR=names(obj_data),total=obj_data,total_cell=colSums(OR_data))
dotplot_data$exp_pattern=tree_anno[match(dotplot_data$OR,tree_anno$OR_gene),]$exp_pattern2
dotplot_data$avg_exp<- dotplot_data$total/dotplot_data$total_cell

dotplot_data$OR_symbol=tree_anno[match(dotplot_data$OR,tree_anno$OR_gene),]$last_name
dotplot_data<- dotplot_data[-1,]
# plot the dot plot (order OR)

pdf("./00_Figure/FigS5-OR_exp_order_dotplot.pdf",width=12,height=6)
# the total count in cells
dotplot_data<- dotplot_data[order(dotplot_data$total,decreasing=TRUE),]
dotplot_data$OR_symbol<- factor(dotplot_data$OR_symbol,levels=unique(dotplot_data$OR_symbol))

ggplot(dotplot_data, aes(x = OR_symbol, y = total, color = exp_pattern)) +
  geom_point(alpha = 0.7) +
  labs(x = "OR", y = "Total UMI", title = "Total UMI of OR gene") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

# the total cell in cells
dotplot_data<- dotplot_data[order(dotplot_data$total_cell,decreasing=TRUE),]
dotplot_data$OR_symbol<- factor(dotplot_data$OR_symbol,levels=unique(dotplot_data$OR_symbol))

ggplot(dotplot_data, aes(x = OR_symbol, y = total_cell, color = exp_pattern)) +
  geom_point(alpha = 0.7) +
  labs(x = "OR", y = "Total total_cell", title = "Total total_cell of OR gene") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

# the avg  in cells
dotplot_data<- dotplot_data[order(dotplot_data$avg_exp,decreasing=TRUE),]
dotplot_data$OR_symbol<- factor(dotplot_data$OR_symbol,levels=unique(dotplot_data$OR_symbol))

ggplot(dotplot_data, aes(x = OR_symbol, y = avg_exp, color = exp_pattern)) +
  geom_point(alpha = 0.7) +
  labs(x = "OR", y = "avg exp", title = "The avg exp of OR gene") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()



pdf("./00_Figure/FigS5-OR_exp_order_ECDF.pdf",width=12,height=6)
ggplot(dotplot_data, aes(x = avg_exp,color = exp_pattern)) +
  stat_ecdf()
  labs(x = "OR", y = "avg exp", title = "The ECDF of avg exp") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

dev.off()













