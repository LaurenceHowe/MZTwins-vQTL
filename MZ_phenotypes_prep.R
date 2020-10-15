##ASSUMES THE FOLLOWING INPUT FILES##

#1) MZ file containing FID and IID of MZ twins
#2) Phenotype file containing phenotypes as well as age and sex with IID linkable to File 1.
#3) Principal components file containing the first 20 PCs with IID linkable to File 1.


##############################################################
#Filepath names, please edit with file paths and output names#
##############################################################

#Paths/names of the three input files
mzname <- paste(" ")
phenname <- paste(" ")
pcname <- paste(" ")

#Paths/names of the three output files: e.g. TEDS.Phenotypes.txt
outphenname <- paste(" ")
outcovname <- paste(" ")
malelistname <- paste(" ")

######################
#Package requirements#
######################

require(data.table)
require(dplyr)
require(tidyr)
require(RNOmni)

#############
#Main script#
#############

#Read in files
mz <- fread(mzname)
phen <- fread(phenname)
pc <- fread(pcname)

#Merge mz/phenotype files
data <- merge(mz, phen, by = "IID") 
#Use below if FID/IID are both in phenotype file
#data <- merge(mz, phen, by = c("FID", "IID"))  

#Add code to remove any exclusions (e.g. non-europeans)

#Double check that there are no singletons
duplicates <- data[duplicated(data$FID), ]
twins <- data[which(data$FID %in% duplicates$FID), ]

#Optional: Restrict to individuals which have PC data. 
twins <- twins[which(twins$IID %in% phen$IID & twins$IID %in% pc$IID), ]

#Generate ordered lists of twins
twin1 <- twins[which(! twins$IID %in% duplicates$IID), ]
twin2 <- twins[which(twins$IID %in% duplicates$IID), ]
twin1o <- twin1[order(twin1$FID), ]
twin2o <- twin2[order(twin2$FID), ]

#Extract FID, Age, Sex and IID of one twin. 
#Note that it's assumed that both twins are genotyped. If only one twin is genotyped will need to tweak script to ensure that IIDs here are from the genotyped twin.
sexage <- data.table(FID = twin1o$FID, IID = twin1o$IID, Sex = twin1o$Sex, Age = twin1o$Age)


#####################################################
#Generate the within-pair mean for each MZ twin pair#
#####################################################

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

#################################
#Generate the MZ pair difference#
#################################

#Again tweak columns as for the within-pair
mz_diff <- abs(twin1o[, -c(1:4)] - twin2o[, -c(1:4)])

#Rename phenotypes such that Height -> Height_df

namelist <- NULL
for (i in phenlist)
{
temp <- paste(i, "df", sep="")
namelist <- cbind(namelist, temp)
}

names(mz_diff) <- namelist


#Extract IDs in same order as dataframe of MZ differences
IDs <- twin1o[, c(1,2)]

#Merge MZ differences with IDs
mz_diff2 <- cbind(IDs, mz_diff) 

#########################################################
#Merge MZ differences with covariates for transformation#
#########################################################

#Merge with principal components, might be easiest to merge on FID in case only one twin is genotyped
mzdiffpc <- merge(mz_diff2, pc, by = "IID")
#mzdiffpc <- merge(mz_diff2, pc, by = c("FID", "IID")

#Merge with sex and age
mzfull <- merge(mzdiffpc, sexage, by = "IID") 
#mzfull <- merge(mzdiffpc, sexage, by = c("FID", "IID")

####################################################
#Generate the covariate-residualised mz diff scores#
####################################################

#Phenotype a) _t includes age,sex,10pcs
#Phenotype b) _nopc includes age,sex
#Phenotype c) _nx includes age,10pcs 

resid_df <- NULL

for (i in 1:length(phenlist))
{
#K may require tweaking.. should indicate 1st phenotype column
#Add 10 PCs and age and sex
k <- i + 4
mzfull$PHEN <- mzfull[, ..k]

##Generate residuals for _t, _nopc and _nx

model_t <- lm(PHEN ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Age + Sex, data = mzfull)
test_t <- mzfull[ , phenlist[i], with = FALSE]
names(test_t) <- c("V1")
test_t$V1[! is.na(test_t$V1)]  <- c(resid(model_t))
  
model_nopc <- lm(PHEN ~ Age + Sex, data = mzfull)
test_nopc <- mzfull[ , phenlist[i], with = FALSE]
names(test_nopc) <- c("V1")
test_nopc$V1[! is.na(test_nopc$V1)]  <- c(resid(model_nopc))
  
model_nx <- lm(PHEN ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Age, data = mzfull)
test_nx <- mzfull[ , phenlist[i], with = FALSE]
names(test_nx) <- c("V1")
test_nx$V1[! is.na(test_nx$V1)]  <- c(resid(model_nx))
  
resid_all <- cbind(test_t, test_nopc, test_nx)

  
resid_df <- cbind(resid_df, resid_all)
}

#Sort out phenotype names to reference the models

namelist2 <- NULL
for (i in phenlist)
{
temp1 <- paste(i, "_t", sep="")
temp2 <- paste(i, "_nopc", sep="")
temp3 <- paste(i, "_nx", sep="")

temp4 <- cbind(temp1, temp2, temp3)

namelist2 <- cbind(namelist2, temp4)
}

names(resid_df) <- namelist2




############################
#Rank normal transformation#
############################

#Uses the CRAN package 'RNOmni' for the transformation

output <- NULL

for (i in 1:length(namelist2))
{

resid_df$PHEN <- resid_df[, ..i]
  
#Restrict to complete cases for each phenotype
temp <- resid_df[! is.na(resid_df$PHEN), ]

#Rank normalise
ranknorm_phen <- rankNorm(temp$PHEN)
resid_df$PHEN2 <- NULL
resid_df$PHEN2[! is.na(resid_df$PHEN)] <- ranknorm_phen

#Extract rank normalised phenotype
out <- resid_df$PHEN2

output <- cbind(output, out)
}

#Generate output which is same as resid_df but rank normalised
output2 <- as.data.frame(output)

#Add the names
names(output2) <- namelist2

################################
#Combine to create output files#
################################

#Generate the phenotype and covariate files
phenotypes <- cbind(sexage, output2)
covariates <- cbind(sexage, ave)

#Also generate list of males for sex-stratified analysis
males <- covariates[which(covariates$Sex == 1), ]
males2 <- males[, c(1:2)]

#Write output files
write.table(phenotypes,outphenname, quote = F, row.names = F)
write.table(covariates,outcovname, quote = F, row.names = F)
write.table(males2, malelistname, quote = F, row.names = F)
