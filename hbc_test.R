library(DESeq2)
library(tximeta)
library(tximport)

samples <- list.files(path = "./hbc_salmon/data", full.names = T, pattern="salmon$")
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "./hbc_salmon/data/", "") %>% 
  str_replace(".salmon", "")
tx2gene <- read.delim("hbc_salmon/tx2gene_grch38_ens94.txt")

txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("tx_id", "ensgene")], countsFromAbundance="lengthScaledTPM")

sampletype <- factor(c(rep("control",3), rep("MOV10_knockdown", 2), rep("MOV10_overexpression", 3)))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))

data <- txi$counts %>% 
  round() %>% 
  data.frame()

dds.hbc <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)
dds.hbc <- DESeq(dds.hbc)

contrast_oe <- c("sampletype", "MOV10_overexpression", "control")
res_tableOE <- lfcShrink(dds.hbc, coef="sampletype_MOV10_overexpression_vs_control", type="apeglm")

res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

mm<-merge(res_tableOE_tb,result_list_anno[[1]],by="geneID")
ggplot(mm,aes(log2FoldChange.x,log2FoldChange.x,colour=threshold_OE))+
  geom_point(alpha=0.4)+theme_bw()

#######
gg="ENSG00000204256"

sig_genes<-result_list_anno[[1]] %>%  
  dplyr::filter(padj<0.05)

sig_genes_hbc<-res_tableOE_tb %>%  
  dplyr::filter(padj<0.05)

all_genes <- result_list_anno[[1]]$geneID
all_genes_hbc <- res_tableOE_tb$gene

nrow(sig_genes_hbc2)

ego <- enrichGO(gene          = as.character(sig_genes_hbc$gene),
                universe      = as.character(all_genes_hbc),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENSEMBL",
                ont           = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = TRUE)

nrow(ego)
dotplot(ego, showCategory=20)

gg=genes %>% dplyr::filter(gene_name=="BRD2") %>% pull(gene_id)
counts(dds,normalized=T) %>% 
  as.data.frame() %>% 
  rownames_to_column("geneID") %>% 
  dplyr::filter(geneID==gg)
