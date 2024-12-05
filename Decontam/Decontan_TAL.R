install.packages("phyloseq")
library(phyloseq)
library(ggplot2)
library(decontam)
library(readr)
#delete the first row and "#OTU ID"
otu<-as.matrix(read.table("raw_otu.tsv",row.names = 1,header = T,sep=""))

# change to csv
tax<-as.matrix(read.csv("taxonomy.csv",row.names = 1,header = T))
class(otu)
class(tax)
sampledate<- read.csv("meta_deco.csv",row.names = 1,header = T)
OTU<-otu_table(otu,taxa_are_rows = TRUE)
TAX<-tax_table(tax)
sampledate<-sample_data(data.frame(sampledate))

#input into phyloseq
sp<- phyloseq(OTU,TAX,sampledate)

# check the labsize
df<-as.data.frame(sample_data(sp))
df$Librarysize <- sample_sums(sp)
df <- df[order(df$Librarysize), ]

sample_data(sp)$is.neg <- sample_data(sp)$Sample_or_Control == "Control"
sample_variables(sp)
contamdf.prev <- isContaminant(sp, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

head(which(contamdf.prev$contaminant))

tax_contaimant_01<- cbind(otu,tax,contamdf.prev)

tax_contaimant_01[tax_contaimant_01=='TRUE'] <- NA
tax_contaimant_01_decontam=na.omit(tax_contaimant_01)

write.csv(tax_contaimant_01,"tax_contaimant_01.csv")
write.csv(tax_contaimant_01,"tax_contaimant_01_decontam.csv")


#################################################

otu<-as.matrix(read.csv("raw_otu.csv",row.names = 1,header = T))
tax<-as.matrix(read.csv("taxonomy.csv",row.names = 1,header = T))
class(otu)
class(tax)


sampledate<- read.csv("meta_deco.csv",row.names = 1,header = T)
OTU<-otu_table(otu,taxa_are_rows = TRUE)
TAX<-tax_table(tax)
sampledate<-sample_data(data.frame(sampledate))
#input into phyloseq
sp<- phyloseq(OTU,TAX,sampledate)

# check the labsize
df<-as.data.frame(sample_data(sp))
df$Librarysize <- sample_sums(sp)
df <- df[order(df$Librarysize), ]



# to idenftify the contaminants by prevalence

sample_data(sp)$is.neg <- sample_data(sp)$Sample_or_Control == "Control"
contamdf<- isContaminant(sp, method="prevalence", neg="is.neg")
table(contamdf$contaminant)


sample_data(sp)$is.neg <- sample_data(sp)$Sample_or_Control == "Control"
sample_variables(sp)
contamdf.prev <- isContaminant(sp, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
#FALSE  TRUE 
#28308    19 
head(which(contamdf.prev$contaminant))

tax_contaimant_01<- cbind(tax,contamdf.prev)
write.csv(tax_contaimant_01,"tax_contaimant_01.csv")

#####try threshold 0.5

contamdf.prev05 <- isContaminant(sp, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
#FALSE  TRUE 
#28157   170  
tax_contaimant_05<- cbind(tax,contamdf.prev05)
write.csv(tax_contaimant_05,"tax_contaimant_05.csv")


