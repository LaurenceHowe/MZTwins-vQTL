# MZTwins-vQTL
<br>

This repository contains two scripts to help analysts run the MZ-differences GWAS.
<br>
The GWAS will use PLINK binary format files (.bed, .bim, .fam).

<br>
1) The R script can be used to generate the phenotypes required for the 4 GWAS models. Details on input files in the R script.
<br>
2) The MZ_GWAS script contains PLINK code for each of the 4 GWAS models using the phenotype files derived in 1).
<br>
<br>

Note that these scripts are not fully automated so will require some tweaking, mostly to file paths and column numbers!
<br>
<br>

The analysis plan and details on selecting phenotypes can be found [here](https://uob-my.sharepoint.com/:w:/g/personal/lh14833_bristol_ac_uk/ETjDr5gvZMBPpnV0WSAcTS8BatLUpd4jYM81vaq9l48-Qw).

<br>
<br>

Please complete the study descriptives file found [here](https://uob-my.sharepoint.com/:x:/g/personal/lh14833_bristol_ac_uk/EVmXGuMddR1Kh5H1TOzDZSoBerEDh7n54Cap4hy5a-paLg?e=8hwe1S), to be summitted along with GWAS summary statistics.

<br>
<br>

<b> Email Laurence.Howe@bristol.ac.uk / E.Assary@qmul.ac.uk if you have any queries. </b>


<b> FREQUENTLY ASKED QUESTIONS </b>



    Q- Missingness: Should I impute the missing data for incomplete twin pairs, when repeated measures from multiple waves of data collection exist? 

A- We don’t recommend imputing missing data. The reason for this is that we are interested in phenotypic differences between twin-pairs. These differences should reflects twin’s phenotypes at the time of the data collection/response, and using data for twins from different waves would not be a true reflection of their differences. For example, for depression, some questionnaires ask about symptoms in the past 6 months. Using a twin pair’s data from two waves a year apart, would mean their differences is due to different timeline being examined, rather than actual discordance. Similarly, with height, using data from different time points for each twin from a pair, would mean their differences reflect the effect of time on growth (the twin measured later will be taller), specially in younger twin samples. So it is important that the twin-pair’s phenotype is from the same wave/timepoint.  

    Q-Mean values for repeated measures: Should I use a mean phenotypic value, when repeated measures are available ? 

A-We don’t recommend creating an average value across timepoints. Partly because of above, in case of missingness across some timepoints. Also because it would be difficult to discern what the “mean difference score” would represent. We recommend choosing the timepoint with the largest complete twin-pair sample.

    Q-Additional covariates: should I include additional study specific covariates such as chip, batch etc ?

A- Yes, please add appropriate study specific covariates at the phenotype preparation step, along with sex, age and 10 PCs.

    Q- Age: Does the age variable reflect current age, or age at the time of data collection?

A-The age covariate should reflect the age of the twins at the time of data collection/response, rather than current age

    Q- Relatedness exclusion: Should I exclude related individuals from my data (e.g. MZ cousins)?

A-we advocate removing one pair for pairs with kinship > 0.1.
