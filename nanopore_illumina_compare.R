library(Matrix)
library(tidyr)
library(dplyr)

setwd("~/data/RNA_modifications/MB_epi")
list.mtx<-list.dirs()
list.mtx<-gsub("\\./","", list.mtx)
list.mtx<-list.mtx[-1]
mtx.path<-file.path(list.mtx, "matrix.abundance.tpm.mtx")
tx<-read.table("7316-1700_kallisto_nac/transcripts.txt")

abundance_matrix<-list()
sample_names<-c()

for (mtx.file in mtx.path) {
  if (file.exists(mtx.file)) {
    print(paste("File exists: ", mtx.file))
    mtx <- readMM(mtx.file)
    sname<-gsub("_kallisto_nac/matrix.abundance.tpm.mtx", "", mtx.file)
    sample_names<-c(sample_names, sname)
    mat<-data.frame(id=tx[,1], sname=t(as.matrix(mtx)))
    abundance_matrix[[sname]]<-mat
  }
  else {
    warning (paste("File not found:", mtx.file))
  }
}

combined_mat<-Reduce(function (x,y) merge(x,y, by="id"), abundance_matrix)
rownames(combined_mat)<-combined_mat$id
combined_mat<-combined_mat[,-1]
colnames(combined_mat)<-sample_names


pbta<-readRDS("/media/pereralab/5706e723-e35e-44e9-b5c5-24cd5bb2f315/OpenPBTA-analysis/data/pbta-isoform-expression-rsem-tpm.stranded.rds")
pbta.meta<-read.table("/media/pereralab/5706e723-e35e-44e9-b5c5-24cd5bb2f315/OpenPBTA-analysis/data/pbta-mb-pathology-subtypes.tsv", header = T, sep = "\t")
pbta.meta<-pbta.meta[pbta.meta$sample_id %in% colnames(combined_mat),]
pbta<-pbta %>%
  separate(transcript_id, into=c("id", "name"), sep = "_")
pbta<-pbta[,-2]
pbta<-pbta[,c("id",pbta.meta$Kids_First_Biospecimen_ID)]
pbta<-pbta %>%
  group_by(id) %>%
  summarise(across(colnames(pbta)[-1], max))
pbta<-data.frame(pbta)
rownames(pbta)<-pbta$id
pbta<-pbta[,-1]
colnames(pbta)<-pbta.meta$sample_id
pbta<-pbta[,colnames(pbta)%in%colnames(combined_mat)]

uid<-unique(gsub("\\..*", "", rownames(pbta)), gsub("\\..*", "",rownames(combined_mat)))

rownames(pbta)<-gsub("\\..*", "", rownames(pbta))
rownames(combined_mat)<-gsub("\\..*", "", rownames(combined_mat))

library(Hmisc)
library(corrplot)

pbta.in<-pbta[uid,]
com.in<-combined_mat[uid,]
#pbta.in[pbta.in>100]<-1000
#com.in[com.in>1000]<-1000

corres<-rcorr(as.matrix(pbta.in), as.matrix(com.in))
cor_mat<-corres$r
p_mat<-corres$P

corrplot(
  cor_mat,
  method = "number",
  type = "lower",
  p.mat = p_mat,
  pch.cex = 1.2,
  pch.col = "darkred",
  mar = c(0,0,2,0)
)

inmat<-data.frame(illumina=pbta.in[,"7316-1700"], nanopore=com.in[,"7316-1700"])
inmat<-log2(na.omit(inmat[rowSums(inmat>2)==2,]))
library(ggplot2)
ggplot(na.omit(inmat), 
       aes(illumina, nanopore)) +
  geom_point() + 
  geom_smooth(method="lm",color="red")

tab<-read.table("7316-1700.hg38.counts.tsv", header = T, sep = "\t")
tab$transcript_name<-gsub("\\..*", "", tab$transcript_name)

cmerge<-merge(tab[,c("transcript_name", "raw")], pbta.in, by.x="transcript_name", by.y=0)
cmerge<-cmerge[,c("transcript_name", "raw", "7316-1700")]

ggplot(cmerge, 
       aes(raw, `7316-1700`)) +
  geom_point() + 
  geom_smooth(method="lm",color="red")









