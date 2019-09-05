#Package requirements
require(data.table)
require(tidyr)
require(dplyr)

#Read in MZ twin list with column names FID, IID
mz<-fread(" ")

#Read in inclusion criteria (e.g. European ancestry) with column names FID, IID
include<-fread(" ", h=F)
mz2<-mz[which(mz$IID%in%include$IID),]

#Simulate example phenotype
mz2$pheno<-rnorm(nrow(mz2),0,1)

#Create MZ twin differences
twins <- mz2 %>%
	arrange(FID) 

dupval <- duplicated(twins$FID) # duplicated values = twins


twins_c1 <- twins[!dupval, ]  %>%
	arrange(FID)

twins_c2 <- twins[dupval, ] %>%
	arrange(FID)
stopifnot(twins_c1$FID == twins_c2$FID) # Checking the couples are in the same order

twins_diff <- twins_c1

# Minus twin 1 phenotypes from twin 2 phenotypes
twins_diff[, -c(1:2)] <- abs(twins_c1[, -c(1:2)] - twins_c2[, -c(1:2)])

##Update sample file##

#Read in .sample file with format "ID_1", "ID_2" etc
sample<-fread(" ")
sample$order<-1:nrow(sample) #Need to preserve sample order

#First allign IDs
twins_diff$FID<-NULL
names(twins_diff)<-c("ID_1", "test")

merge<-merge(sample,twins_diff,by="ID_1", all=T)

#Preserve sample order
out<-merge %>%
	arrange(order)

out$order<-NULL

#Label phenotypes in sample file
out$pheno[out$ID_1==0]<-"P"

#Output

write.table(out, "mztwins.sample", quote=F, row.names=F)
