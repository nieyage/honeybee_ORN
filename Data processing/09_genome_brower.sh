
class WaterMelon:
    def __init__(self,gtf,fasta):
        self.gtf = gtf
        self.fasta = fasta
        
    def chromoLen(self):
        chrLen = {}
        for rec in SeqIO.parse(self.fasta,'fasta'):
            chrLen[rec.id] = len(rec.seq)
            
        return chrLen
    
    def gcContent(self,gc_window):
        self.gc_window = gc_window
        
        gc_content = {'chr_id':[],
                     'bin_start':[],
                     'gc':[]}
        chr_len = self.chromoLen()
        
        #print(chr_len)
        for rec in SeqIO.parse(self.fasta,'fasta'):
            for i in range(0,chr_len[rec.id],self.gc_window):
                #print(rec.id)
                gc_content['chr_id'].append(rec.id)
                gc_content['bin_start'].append(i)
                gc_content['gc'].append(round(GC(rec.seq[i:i+self.gc_window]),2))
        return pd.DataFrame(gc_content)
    
    def geneDensity(self,gene_window):
        self.gene_window = gene_window
        final_df = []
        df = pd.read_table(self.gtf,header=None,comment="#",sep="\t",
                           usecols=[0,2,3,4],
                           names="Chromosome Feature Start End".split())
        #df.columns = "Chromosome Source Feature Start End Score Strand Frame Attribute".split()
        df = df[df.Feature=="gene"]
        chrLen = self.chromoLen()
        for chr_id in chrLen.keys():
            print(chr_id)
            df1 = df[df.Chromosome==chr_id]
            gene_start = [int(a) for a in df1.Start]
            gene_start.insert(0,0)
            gene_start.append(round(chrLen[chr_id]/self.gene_window)*self.gene_window+self.gene_window)
            #print(gene_start)
            bin_start = [int(a.left) for a in pd.cut(gene_start,bins=round(chrLen[chr_id]/self.gene_window)+1).value_counts().index]
            bin_start[0] = 0
            gene_count = list(pd.cut(gene_start,bins=round(chrLen[chr_id]/self.gene_window)+1).value_counts().values)
            #print(len(bin_start))
            #print(len(gene_count))
            #print("OK")
            final_df.append(pd.DataFrame({'chr_id':chr_id,'bin_start':bin_start,'gene_count':gene_count}))
            
        return pd.concat(final_df)
        
  



# get the AT repeat file and upload in UCSC 


# set the range chrLG2:9970000-10160000

awk '$1 == "chr1" && $4 >= 1 && $5 <= 10000 && $3 == "gene"' yourfile.gtf > extracted_genes.gtf

library(ggplot2)

# 读取GTF文件，假设GTF文件名为"yourfile.gtf"
gtf_data <- read.table("/Users/fraya/Documents/project/honeybee/chemoreceptor/genes.gtf", header = FALSE, sep = "\t", quote = "")

# 指定要绘制的区域范围，这里假设绘制chr1的1-10000范围内的基因
chr <- "Group2"
start_pos <- 9970000
end_pos <- 10160000

# 筛选符合区域范围的基因
genes <- subset(gtf_data, V1 == chr & V4 >= start_pos & V5 <= end_pos & V3 == "transcript")
# 创建基因位置图并禁用图例
gg <- ggplot(data = genes, aes(x = (V4 + V5) / 2, y = V5, shape = V7, color = V9)) +
    geom_jitter(height = 1000, width = 100, size = 5) +
    scale_shape_manual(values = c(25,24)) +
    theme_minimal() +
    labs(x = "Genomic Position", y = "") +
    theme(axis.text = element_blank(), axis.title = element_text(size = 12)) +
    guides(shape = FALSE, color = FALSE)  # 禁用形状和颜色图例

# 打印图形
print(gg)


library(Biostrings)

genome <- readDNAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa")
gc_content <- letterFrequencyInSlidingView(genome,"GC" )

LOC1007966102