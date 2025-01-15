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

nc_files<-list.files(pattern = "hg38.NanoCount.tsv")
nc_list<-lapply(nc_files,function (x) read.table(x, header = T, sep = "\t"))
nc_mat<-Reduce(function(x,y) merge(x,y, by="transcript_name", all=T), nc_list)
nc_mat<-nc_mat[,c(1, grep("tpm", colnames(nc_mat)))]
nc_mat[is.na(nc_mat)]<-0
nc_snames<-gsub(".hg38.NanoCount.tsv", "", nc_files)
colnames(nc_mat)<-c("transcriipt_name", nc_snames)
nc_mat$transcriipt_name<-gsub("\\..*", "", nc_mat$transcriipt_name)
rownames(nc_mat)<-nc_mat$transcriipt_name
nc_mat<-nc_mat[,-1]
counts_df3 <- calculate_gene_counts(nc_mat)
# Combine counts into a table
result_table <- data.frame(
DataFrame = c(rep("Illumina", length(counts_df1)),
rep("lr-kallisto", length(counts_df2)), #,
rep("NanoCount", length(counts_df3))),
Column = c(names(counts_df1), names(counts_df2), names(counts_df3)),
Count = c(counts_df1, counts_df2, counts_df3)
)
# Print the table
print(result_table)
cmerge<-merge(nc_mat[,c("7316-1700")], pbta.in, by=0)
head (cmerge)
cmerge<-merge(nc_mat, pbta.in, by=0)
cmerge
cmerge<-cmerge[,c("Row.names", "7316-1700.x", "7316-1700.y")]
cmerge<-merge(cmerge, com.in, by.x="Row.names", by.y=0)
head (cmerge)
cmerge<-cmerge[,c(1,2,3,6)]
rownames(cmerge)<-cmerge$Row.names
cmerge<-cmerge[,-1]
colnames(cmerge)<-c("NanoCount", "Illumina", "lr-kallisto")
cmerge<-log2(na.omit(cmerge[rowSums(cmerge>0.1)==3,]))
corres<-rcorr(as.matrix(cmerge))
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
ggplot(cmerge,
aes(`lr-kallisto`, Illumina)) +
geom_point() +
geom_smooth(method="lm",color="red")
ggplot(cmerge,
aes(NanoCount, Illumina)) +
geom_point() +
geom_smooth(method="lm",color="red")
cor.test(cmerge$Illumina, cmerge$NanoCount)
cor.test(cmerge$Illumina, cmerge$`lr-kallisto`)









