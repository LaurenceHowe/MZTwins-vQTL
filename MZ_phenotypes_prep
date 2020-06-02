##ASSUMES THE FOLLOWING INPUT FILES##

#1) MZ file containing FID and IID of MZ twins
#2) Phenotype file containing phenotypes as well as age and sex with IID linkable to File 1.
#3) Principal components file containing the first 20 PCs with IID linkable to File 1.

#Filepath names, please edit with file paths and output names
mzname <- paste(" ")
phenname <- paste(" ")
pcname <- paste(" ")

outphenname <- paste(" ")
outcovname <- paste(" ")
malelistname <- paste(" ")

#Package requirements
require(data.table)
require(dplyr)
require(tidyr)

#Read in files
mz <- fread(mzname)
phen <- fread(phenname)
pc <- fread(pcname)

#Merge mz/phenotype files
data <- merge(mz, phen, by = "IID") 
#Use below if FID/IID are both in phenotype file
#data <- merge(mz, phen, by = c("FID", "IID"))  


#Double check that there are no singletons
duplicates <- data[duplicated(data$FID), ]
twins <- data[which(data$FID %in% duplicates$FID), ]

#Generate ordered lists of twins
twin1 <- twins[which(! twins$IID %in% duplicates$IID), ]
twin2 <- twins[which(twins$IID %in% duplicates$IID), ]
twin1o <- twin1[order(twin1$FID), ]
twin2o <- twin2[order(twin2$FID), ]

#Extract FID, Age, Sex and IID of one twin. 
#Note that it's assumed that both twins are genotyped. If only one twin is genotyped will need to tweak script to ensure that IIDs here are from the genotyped twin.
info <- data.table(FID = twin1o$FID, IID = twin1o$IID, Sex = twin1o$Sex, Age = twin1o$Age)

#Generate the within-pair mean
#In the case of "-c(1:4)", columns 5 and onwards are the phenotypes which need to be averaged (columns 1-4 are FID, IID, Age and Sex)
# will require editing.  
ave <- 0.5 * (twin1o[, -c(1:4)] + twin2o[, -c(1:4)])
#Again tweak to indicate only phenotype columns
phenlist <- names(twin1o)[5:10]

#Change the names of these variables: Height -> HeightM
namelist <- NULL
for (i in phenlist)
{
temp <- paste(i, "M", sep="")
namelist <- cbind(namelist, temp)
}

names(ave) <- namelist
				
#Generate the pair difference
#Again tweak columns as for the within-pair
mz_diff <- abs(twin1o[, -c(1:4)] - twin2o[, -c(1:4)])

#Generate the PC residualised phenotypes
datapc <- merge(twins, pc, by = "IID")
#datapc <- merge(twins, pc, by = c("FID", "IID")


resid_df <- NULL
for (i in 1:length(phenlist))
{
#K may require tweaking.. should indicate 1st phenotype column
#Add number of PCs depending on study, only 200 twin pairs in UKB so went with 10.
k <- i + 4
datapc$PHEN <- datapc[, ..k]
model <- lm(PHEN ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = datapc)

test <- datapc[ , phenlist[i], with = FALSE]
names(test) <- c("V1")
test$V1[! is.na(test$V1)]  <- c(resid(model))

resid_df <- cbind(resid_df, test)
}



#Again create a list of names: Height -> Height_resid
namelist2 <- NULL
for (i in phenlist)
{
temp <- paste(i, "_resid", sep="")
namelist2 <- cbind(namelist2, temp)
}

names(resid_df) <- namelist2

#Combine to create output files
phenotypes <- cbind(info, mz_diff, resid_df)
covariates <- cbind(info, ave)

#Also generate list of males for sex-stratified analysis
males <- covariates[which(covariates$Sex == 1), ]
males2 <- males[, c(1:2)]

write.table(phenotypes,outphenname, quote = F, row.names = F)
write.table(covariates,outcovname, quote = F, row.names = F)
write.table(males2, malelistname, quote = F, row.names = F)
