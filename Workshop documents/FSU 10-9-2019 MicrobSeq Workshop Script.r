
######################################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
######################################################################################################################################################
##																																					##
##                                                   Florida State University 																		##
##                                                   Analyzing Microbial Sequencing Data in R 													    ##
##                                                   October 9, 2019 																			    ##
##                                                   Brian D. Piccolo, Ph.D.                                                                        ##
##																																					##
######################################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
######################################################################################################################################################

## Requirements:

## Review Chapters 1-3 in Web-book (https://bdpiccolo.github.io/FSU-Microbial-Data-Analysis/).  

## R version >= 3.6.0

## Data files can be found at the following website: https://github.com/bdpiccolo/FSU-Microbial-Data-Analysis/tree/master/Data
## UCDT2DM16SExcel.xlsx
## UCDT2DMmetadata.xlsx

## Text preceded by hashtag, #, do not have to be run in the console since a # tells R not to interpret code line.  But running these lines will do no harm.

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##																													##
##                     ~~~~~~~~~~~~~~~~~~~~~~        Introduction         ~~~~~~~~~~~~~~~~~~~~~~ 					##
##																													##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## This is a modified script that will be demonstrated during the Workshop.  It contains the main concepts from 
## Chapters 4-7 from the web-book (https://bdpiccolo.github.io/FSU-Microbial-Data-Analysis/) associated with this 
## course.  This is to ensure we can cover the main themes from the web-book during the 3-hour workshop.  This script 
## has some notes associated with the coding, but it is assumed that R concepts have been reviewed in Chapters 1-3. 
## No new coding will be introduced in this script; it is directly taken from the web-book.  We will use Excel files 
## to upload data, since most nutrition researchers are familiar with this format.  Coding for files derived from 
## microbial informatics pipelines, e.g., Biome and QZA, can be found in the web-book.

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##																													##
##                 ~~~~~~~~~~~~~~~~~~~~~~        Required R Packages          ~~~~~~~~~~~~~~~~~~~~~~ 				##
##																													##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

# Run the following code in your R console will identify whether the listed packages are already installed on your computer.
# If a package is not found, then it will be downloaded

# CRAN packages
install_CRANpackages <- c("tidyverse","readxl", "zCompositions","rgr","vegan","ape")
if (length(setdiff(install_CRANpackages, rownames(installed.packages()))) > 0) {
	install.packages(setdiff(install_CRANpackages, rownames(installed.packages())))
}

# Bioconductor packages
install_BioCpackages <- c("BiocManager", "phyloseq", "microbiome","DESeq2","metagenomeSeq", "ALDEx2")
if (length(setdiff(install_BioCpackages, rownames(installed.packages()))) > 0) {
	BiocManager::install(setdiff(install_BioCpackages, rownames(installed.packages())), ask = FALSE)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##																													##
##                ~~~~~~~~~~~~~~~~~~~~~~        Set Working Directory         ~~~~~~~~~~~~~~~~~~~~~~ 				##
##																													##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

# Sets the file directory that R will look at.  Should be to the file where you will import and export files, or at least the parent directory.

setwd("C:\\my\\folder\\directory\\")
# or
setwd("C:/my/folder/directory/")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##																													##
##                  ~~~~~~~~~~~~~~~~~~~~~~        Import Files         ~~~~~~~~~~~~~~~~~~~~~~ 						##
##																													##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

# R has traditionally utilized .txt or .csv (tab or comma delimited text files) for data importation.
# As most of the colleagues utilized Excel format, which requires convert Excel to .txt or .csv file.  
# Recent contributions to R have made it much easier to import Excel files.
# I try to minimize as much formatting as possible, so we will import Excel files directly.

#####################
# Load required packages
library(readxl) # required for Excel import
library(tidyverse) # Convenience functions for piping (%>%), reshaping data, and summarizing data.

#####################
# Import Count and Taxonomy data

# Need to specify data type if numbers and characters are mixed. See help page for read_excel().
Excel <- read_xlsx(
	path="./Data/UCDT2DM16SExcel.xlsx",  
	col_types=c(rep("text", 8), rep("numeric", 55))
) %>% as.data.frame() # The 'pipe' transfers the results of a function to the first argument of the next function.

# If you're having trouble loading your data, it might be due to this coding: './Data/
# I've included it to tell R that the data is in a sub-directory.  
# Delete it if your files are in the directory you specified above with setwd().

# Preview at the data
head(x=Excel)

# Look at structure of data
str(object=Excel)

# Look at dimension of data
dim(x=Excel)

# We eventually will be creating a phyloseq object with separate parts of these data.
# phyloseq requires the row-name attribute to match between the files.
# Use the OTU identifier as the row names of this data frame.
rownames(x=Excel) <- Excel$OTU

# Remove OTU column
Excel$OTU <- NULL

# The phyloseq object will require 3 separate objects: A count file, the taxonomy, and the sample meta data.
# Therefore, we need to separate the count and taxonomy data from the current data frame 

# Subset columns related to taxonomy
TaxaData <- Excel[, colnames(Excel) %in% c("Kingdom","Phylum","Class","Order","Family","Genus","Species")]

# Convert to a matrix
TaxaMatrix <- TaxaData %>% as.matrix()

# Look at dimension of data
dim(x=TaxaMatrix)

# Use '!()' to reverse logical statements
CountData <- Excel[, !(colnames(Excel) %in% c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))]

# Convert to a matrix
CountMatrix <- CountData %>% as.matrix()

# Look at dimension of data
dim(x=CountMatrix)

#####################
# Import Sample Metadata and convert to data frame
Metadata <- read_xlsx(path="./Data/UCDT2DMmetadata.xlsx") %>% as.data.frame()

# Set row names with column that matches the labels in your count data
# In this case, it is the column labeled "SampleID".
# This is required to create a phyloseq object.
rownames(x=Metadata) <- Metadata$SampleID

# Look at data
head(x=Metadata)

# Look at structure of data
str(object=Metadata)

# Look at dimensions of data
dim(x=Metadata) # Note that the number of rows match the number of columns in CountMatrix.

#####################
## Ensure the names match between the count and metadata file.

# Extract column names
Count_colnames <- colnames(x=CountData)

# Extract 'SampleID' file.  The label, "SampleID", was the name I provided, but yours may be different.
Metadata_rownames <- rownames(x=Metadata)

# Sort both
Count_colnamesSort <- sort(x=Count_colnames)
Metadata_rownamesSort <- sort(x=Metadata_rownames)

# Use == to ensure that each ordered ID is equal - output should be a vector of TRUE values.  A single FALSE means they are not equal
Count_colnamesSort == Metadata_rownamesSort
# or check with the setequal() function.  Will output either a TRUE or FALSE.  TRUE means they match exactly, FALSE means they do not.
setequal(Count_colnamesSort == Metadata_rownamesSort, TRUE)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##																													##
##                  ~~~~~~~~~~~~~~~~~~~~~~        Make phyloseq object         ~~~~~~~~~~~~~~~~~~~~~~ 				##
##																													##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

# The minimum data required to create a phyloseq object are a Count Table, Taxonomy Table, and Sample Metadata.
# The row and column names of each of these files must match as follows:

# 1. Count Table row names ->  Taxa Table row names.

# 2. Count Table column names -> Sample metadata row names

# We've already determined that the Count Table  column names match with the sample metadata, so lets check the Count table and taxa table row names.

# Extract row names from count data
Count_rownames <- rownames(x=CountData)

# Extract row names from taxa data
Taxa_rownames <- rownames(x=CountData)

# Sort both
Count_rownamesSort <- sort(x=Count_rownames)
Taxa_rownamesSort <- sort(x=Taxa_rownames)

# Make object of logical output using ==
CountTaxaCompLogical <- Count_rownamesSort == Taxa_rownamesSort

# Use TRUE as comparison vector (TRUE is recycled!) in the setequal() function.
setequal(CountTaxaCompLogical, TRUE)

#####################
# Load required packages
library(phyloseq)
library(microbiome)

# Convert OTU data into phyloseq otu_table class
otuTABLE <- otu_table(object=CountMatrix, taxa_are_rows = TRUE)

# Convert taxa data into phyloseq taxonomyTable class
taxTABLE <- tax_table(object=TaxaMatrix)

# Convert Sample Metadata data frame into phyloseq sample_data class
sampleDATA <- sample_data(object=Metadata)

# Create phyloseq object. Order of objects in phyloseq() function does not matter.
phylo_obj <- phyloseq(otuTABLE, taxTABLE, sampleDATA)
phylo_obj

#####################
# Extracting components of phyloseq object.

# Extract OTU table using abundances()
OTUdata <- abundances(phylo_obj)

# Extract Sample Metadata using meta()
SampleData <- meta(phylo_obj)

# microbiome package does not have an extractor function for taxonomy data
# Will use phyloseq extractor, tax_table().  However, tax_table() does not return it as a data frame.  
# We need to coerce to a data frame.  The '@.Data' is the location within tax_table() output that contains the actual data.
TAXAData <- as.data.frame(tax_table(phylo_obj)@.Data)
# or
TAXAData <- tax_table(phylo_obj)@.Data %>% as.data.frame()

#####################
# Subsetting 

# Rats in the UCD-T2DM study were collected in 2014 and in 2016.  We're worried about a potential batch effect. Thus, we wish to remove 
# 2014 samples and only examine 2016 samples.  We can extract the metadata and define a logical statement to keep the rats sampled in 2016.  

# We extracted the Sample Metadata above.
CollectionYear <- SampleData$Collection
CollectionYear

# Make the logical statement using the %in% (match) operator
CollectionYearLogic <- CollectionYear %in% "Y16"
CollectionYearLogic

# Use the prune_samples() function
# First argument is the logical statement, the second is the phyloseq object
phylo_obj2016 <- prune_samples(samples=CollectionYearLogic, x=phylo_obj)
# or
phylo_obj2016 <- prune_samples(samples={meta(phylo_obj)$Collection %in% "Y16"}, x=phylo_obj) 
# The curly brackets help contain the logical statement.  Useful if there is a long string of code.

# The 'pruned' phyloseq object
phylo_obj2016

# While the lean Sprague Dawley rats were used to compare metabolic differences between the diabetic animals, they may not be 
# appropriate control animal for the microbiota study due to genetic influence on the composition of the gut microbiota. 
# We need to remove this group from our 2016 object.

# Extract the metadata from phylo_obj2016
Group2016 <- meta(phylo_obj2016)$GroupAbbrev
Group2016

# Make the logical statement using the != (does not equal) operator
Group2016Logic <- Group2016 != "LSD"
Group2016Logic

# Use the prune_samples() function
# First argument is the logical statement, the second is the phyloseq object
phylo_obj2016_noLSD <- prune_samples(samples=Group2016Logic, x=phylo_obj2016)
# or
phylo_obj2016_noLSD <- prune_samples(samples={meta(phylo_obj2016)$GroupAbbrev != "LSD"}, x=phylo_obj2016)

# The 'pruned' phyloseq object
phylo_obj2016_noLSD

# Beware, once you subset a phyloseq object, there is no going back.  

#####################
# Sample depth

# We can also subset samples based on the characteristics of the OTU table.  An important example of this is when we check sample depths.  

# The sample_sums() function will calculate sample depths from your phyloseq object
phyloSampDepth <- sample_sums(x=phylo_obj)
phyloSampDepth

# What is our total depth?
sum(phyloSampDepth)

# Let's filter depths < 40,000.  ***DO NOT USE THIS AS A CUTOFF GUIDE - FOR PURPOSE OF EXAMPLE ONLY***
# We want to keep samples that have sequence depths greater than 40000
Depth40K <- phyloSampDepth > 40000

# Use the prune_samples() function
phylo40K <- prune_samples(samples=Depth40K, x=phylo_obj)
# or
phylo40K <- prune_samples(samples={sample_sums(phylo_obj) > 40000}, x=phylo_obj)

# The 'pruned' phyloseq object
phylo40K

#####################
# Taxonomy

# We have >50,000 taxa in the UCD-T2DM and we're worried that we sequenced too deeply.  Many of these taxa could be spurious, 
# especially if there are only a few reads across all samples.  Thus, we want to filter any taxa that doesn't sum to 50 reads.  

# The taxa_sums() function will calculate sample depths across taxa 
phyloTaxaDepth <- taxa_sums(x=phylo_obj)

# We want to keep samples that have sequence depths greater than 50 (this is abitrary and selected for the purpose of this example).
TaxaDepth50 <- phyloTaxaDepth > 50

# Use the prune_samples() function
phyloTaxa50 <- prune_taxa(taxa=TaxaDepth50, x=phylo_obj)

# The 'pruned' phyloseq object
phyloTaxa50

#####################
# Aggregating

# Many early investigations relating the microbiota to obesity identified an increased ratio of Firmicutes to Bacteroidetes in obese mice 
# relative to their lean counterparts.  This is at the phylum level, so the most broad phylogenetic level of bacteria.  
# Let's aggregate the phyloseq object to the phylum level.

# Aggregate the phyloseq object to Phylum level. 
# Will take 5-10 seconds to process because of the large amount of taxa
phylo_Phylum <- tax_glom(physeq=phylo_obj, taxrank="Phylum")

# The 'aggregated' phyloseq object
phylo_Phylum

# Use the tax_table() function from the phyloseq package to view aggregated taxa table. 
tax_table(phylo_Phylum)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##																													##
##      ~~~~~~~~~~~~~~~~~~~~~~        Pre-processing microbial sequencing data         ~~~~~~~~~~~~~~~~~~~~~~ 		##
##																													##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

# Pre-processing is arguably the most important step in your data analysis pipeline.  We have already been exposed to some 
# pre-processing examples, e.g.,  filtering low abundant reads and identifying samples with low read counts.  

#####################
# Filtering low abundance taxa

# Identify a minimum detection limit within a sample and then identify a whether a certain percentage of your samples 
# meet that criteria. 

# Use the core() function from the microbiome package 
# The `core` function requires 3 arguments.  The first is the phyloseq object, the second defines the 'detection' limit, 
# and the third defines the 'prevalence' (i.e., proportion) that meets the defined limit.  

# Filter taxa so that taxa have at least 10 reads in 25% of the samples.
phyobj_d10p25 <- core(x=phylo_obj, detection = 10, prevalence = 0.25)

# View object
phyobj_d10p25

#####################
# Proportional Abundance

# Proportional abundance, also sometimes referred to relative abundance, can applied to a phyloseq object.  
# Proportional abundances will sum to 1, or 100 if expressed as percentages.

# Use the microbiome::transform() 
phyobj_d10p25_prop <- microbiome::transform(x=phyobj_d10p25, transform="compositional")

# View a subset of the data
head(x=abundances(x=phyobj_d10p25_prop)[,1:10])
# or 
abundances(x=phyobj_d10p25_prop)[,1:10] %>% head()

# Do the taxa sum to 1?
sample_sums(x=phyobj_d10p25_prop)

# It is also common to filter taxa based on proportional abundances.  Many times, a minimum average or median proportion must be met.  
# Say, 1% on a scale of 0-100% (equivalent to 0.01 on a 0 to 1 scale).  We will use the `prune_taxa` function to do this

# Transform counts to proportional abundances.
phylo_obj_PA <- microbiome::transform(x=phylo_obj, "compositional")

# Determine mean proportional abundances by taxa.
phylo_obj_PA_taxamean <- rowMeans(x=abundances(x=phylo_obj_PA))
# or 
phylo_obj_PA_taxamean <- abundances(x=phylo_obj_PA) %>% rowMeans()

# Use logical statement in prune_taxa()Using prune_taxa() function.
phyobj_PA_p01 <- prune_taxa(taxa={phylo_obj_PA_taxamean > 0.01}, x=phylo_obj_PA)

# View phyloseq object
phyobj_PA_p01

#####################
# Handling Zero

# Can shift all counts with a positive number to remove zeros
phyobj_d10p25_shift <- microbiome::transform(x=phyobj_d10p25, transform="shift", shift=1)

# Compare abbreviated taxa
data.frame(Original=abundances(x=phyobj_d10p25)[1,1:15], Shifted=abundances(x=phyobj_d10p25_shift)[1,1:15])

# Bayesian-multiplicative treatment.  

# In order to perform this zero-imputation method, we will need to use the `cmultRepl` function from the zCompositions package. 
# The first argument of `cmultRepl` takes your data, then we will input "CZM" (count zero multiplicative) in the 'method' argument, 
# and "p-counts" (pseudo-counts) in the 'output' argument. We have to do a bit of data wrangling to get the treated data back 
# into the phyloseq object.  

# Load library
library(zCompositions)

# Extract OTU table
d10p25_OTU <- abundances(x=phyobj_d10p25)

# Taxa needs be rows and samples need to be columns
d10p25_OTU_BMz <- cmultRepl(X=d10p25_OTU, method="CZM", output="p-counts")

# Make duplicate phyloseq object with new name
phyobj_d10p25_BMz <- phyobj_d10p25

# Replace OTU table data with zero imputed data
otu_table(object=phyobj_d10p25_BMz) <- otu_table(object=d10p25_OTU_BMz, taxa_are_rows = TRUE)

# View abbreviated data
head(x=abundances(phyobj_d10p25_BMz)[,1:15])

# Be careful with the phyloseq objects that have had zero treatments and/or transformations.  Some of the functionalities may 
# give incorrect or incosistent results, e.g., `tax_glom`.  It is best to subset or aggregate prior to these treatments.  
# If these treatments are required to make a subset/aggregation decision, then I advise making 2 separate phyloseq objects as we 
# did in the above example.  

#####################
### Centered log ratios

# A centered-log ratio (CLR) transformation has been suggested to render compositional data compatible with standard multivariate 
# techniques.  However, the individual counts become ratios between all parts of the sample, so the interpretion is now altered.

# We will use the `clr` function from the rgr package to carry out this function.  
# The CLR transformation requires a counts or proportions values > 0.  Therefore, we will use the `phyobj_d10p25_BMz` object in the analysis.

# Load library
library(rgr)

# Extract OTU table from the zero imputed phyloseq object.
BMz_OTU <- abundances(phyobj_d10p25_BMz)

# Taxa needs be rows and samples need to be columns
BMz_clr_OTU <- clr(xx=BMz_OTU, ifwarn=FALSE)

# Make duplicate phyloseq object with new name
phyobj_BMz_clr <- phyobj_d10p25_BMz

# Replace OTU table data with zero imputed data
otu_table(object=phyobj_BMz_clr) <- otu_table(object=BMz_clr_OTU, taxa_are_rows = TRUE)

# View abbreviated data
head(x=abundances(phyobj_BMz_clr)[,1:15])

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##																													##
##                   ~~~~~~~~~~~~~~~~~~~~~~        Alpha-Diversity         ~~~~~~~~~~~~~~~~~~~~~~ 					##
##																													##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

#####################
## Alpha diversity measurements in phyloseq package

# estimate_richness()
UCDrats_Adiv <- estimate_richness(physeq=phyobj_d10p25)

# View object
head(x=UCDrats_Adiv)

#####################
## Alpha diversity measurements in microbiome package

# richness measurements
UCDrats_Adiv_rich <- richness(x=phyobj_d10p25)

# View abbreviated object
head(x=UCDrats_Adiv_rich)

# dominance measurements
UCDrats_Adiv_dom <- dominance(x=phyobj_d10p25)

# View abbreviated object
head(x=UCDrats_Adiv_dom)

# evenness measurements
UCDrats_Adiv_even <- evenness(x=phyobj_d10p25)

# View abbreviated object
head(x=UCDrats_Adiv_even)

# rarity measurements
UCDrats_Adiv_rar <- rarity(x=phyobj_d10p25)

# View abbreviated object
head(x=UCDrats_Adiv_rar)

#####################
# Data analysis

# Now that we know how to estimate a wide range of alpha-diversity measurements, we now want to know whether the indices are 
# altered by our experimental treatments, e.g., diet intervention.  

# We will continue to work with the UCD-T2DM rat dataset.  Remember, there are 5 groups, but also 2 collection periods.  
# We will use ANOVA with a main effect of group and block by collection year.

# Let's see if there is a difference in Chao1.  First, we need to combine the diversity measurements with our Sample Metadata.  

# Extract Sample Metadata
UCDrats_MetaData <- meta(x=phyobj_d10p25)

# Make matching columns for both Alpha-diversity and metadata data frames using the row names
# Sample Metadata already has matching column, but will need to create if not in data frame.
UCDrats_Adiv$SampleID <- rownames(x=UCDrats_Adiv)

# Join alpha-diversity and metadata
UCDrats_Adiv_DF <- full_join(x=UCDrats_MetaData, y=UCDrats_Adiv, by="SampleID")

# View data
head(x=UCDrats_Adiv_DF)

# Plot Chao1 by group.  
# Boxplot by Groups
boxplot(formula=Chao1 ~ GroupAbbrev, data=UCDrats_Adiv_DF)

# Boxplot by collection year
boxplot(formula=Chao1 ~ Collection, data=UCDrats_Adiv_DF)

# Looks like the pre-diabetics (PD) and recent-diabetics (RD) rats may have higher diversity compared to the rats with more 
# advanced diabetes (D3M and D6M).  Also appears that rats collected in 2016 may have higher diversity compared to those collected 
# in 2014.  Let's see if these differences are statistically significant.

# Assess ANOVA on Chao1
Chao1ANOVA <- aov(formula=Chao1 ~ GroupAbbrev + Collection, data=UCDrats_Adiv_DF)

# View ANOVA table
summary(object=Chao1ANOVA)

# Both 'GroupAbbev' and 'Collection' effects are significant.  Now, which groups differ from each other.

# Assess Tukey HSD test.
TukeyHSD(x=Chao1ANOVA, which="GroupAbbrev")

# The Chao1 index in RD rats is greater than in both the D3M and D6M rats, while the Chao1 index in PD is higher than D3M and 
# approaches significance relative to D6M.  

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##																													##
##                   ~~~~~~~~~~~~~~~~~~~~~~        Beta-Diversity         ~~~~~~~~~~~~~~~~~~~~~~ 					##
##																													##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

# There are distance-, dissimilarity-, and covariance-, and phylogenetic-based estimations of beta-diversity, all of which
# may uncover different aspects of your data. On top of that, normalization procedures will also influence the results. 
# The first step is to calculate a matrix of your beta-diversity index.

# The `phyloseq::distance` can calculate several beta diveristy indices.  Requires 2 arguments, the phyloseq object and the 
# name of the index.  The indices names can be found with the distanceMethodList function.

# View available beta-diversity indices
distanceMethodList

# The 'Unifrac', 'dpcoa', and 'JSD' require a phylogenetic tree and are computationally more expensive than the others.  
# The indices under 'vegdist' and 'betadiver' are derived from ecology sciences and implemented in the `vegdist` function 
# from the vegan package.  

# The most common diversity indices in microbial based studies tend to be the Unifrac distances, Bray-Curtis Dissimilarities, 
# and the Jaccard index.  We will focus on Bray-Curtis and Jaccard for example purposes. 

#####################
# Beta-diversity using phyloseq package

# Returns a matrix, which won't be shown for the sake of brevity.
phylo_BrayDis <- phyloseq::distance(physeq=phyobj_d10p25, method="bray")

# Jaccard index.
phylo_Jaccard <- phyloseq::distance(physeq=phyobj_d10p25, method="jaccard")

#####################
# Beta-diversity using vegan package

# Load package
library(vegan)

# Extract OTU table.
UCDcounts <- abundances(phyobj_d10p25)

# Samples need to be rows for vegdist() to calculate between sample distances/dissimilarities.
tUCDcounts <- t(UCDcounts)

# Bray index.
vegan_BrayDis <- vegdist(tUCDcounts, "bray")

# Are the two distance matrices equal?  They should be, the `phyloseq::distance` uses `vegdist` underneath it's hood.  

# Use setequal() to test whether all elements are TRUE
setequal(phylo_BrayDis == vegan_BrayDis, TRUE)

#####################
# Ordinations using phyloseq package

# The ordinate() phyloseq function currently supports the following ordinations: "DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"
# It is mostly a wrapper for ordinations in other packages, which means it does not contain the coding itself to perform the ordination.  
# It has to call functions from other packages to conduct the ordination.  Please see the help page for more information: ?ordinate  
# However, it is very convenient because the `ordination` function preps the data based on the required syntax of the called ordination function.  
# Thus, there is no need for additional coding. 

# ordinate() calls the pcoa() function from the ape package.
UCD_bray_PCoA <- ordinate(physeq=phyobj_d10p25, method="PCoA", distance="bray")

# ordinate() calls the metaMDS() function from the vegan package.
UCD_jac_NMDS <- ordinate(physeq=phyobj_d10p25, method="NMDS", distance="jaccard")

# PCoA output is very verbose and skipped here. 
# View NMDS output.   
UCD_jac_NMDS

#####################
# Ordinations using vegan package

### vegan

# We will look closer at the metaMDS() function in the vegan package.  
# There are other ordination functions in vegan (e.g., cca(), rda(), capscale()) that can be used.
# metaMDS() has several arguments that allow extensive customization of the ordination. 

# The first argument, 'comm', can either take your actual data (Counts, proportional abundance, transformed counts, etc.) 
# or a symmetric square distance/dissimilarity matrix.  

# The second argument identifies the distance/dissimilarity index if data is loaded instead of a distance/dissimilarity matrix.  
# I recommend inputting a distance/dissimilarity matrix because you have finer control over the input.  

# The next argument, 'k', sets the number of dimensions. It is generally recommended to use the minimum number of dimensions, 
# unless you are not finding convergence solutions.  The help page recommmends increaseing k by 1 to help find convergence solutions. 

# The next arguments, 'try' and 'trymax' set the minimum and maximum number of iterations. 

# The 'engine' argument defines the MDS variant to run.  Leave this alone unless you have a good reason to alter.  

# Lastly, the 'autotransform' arguments defaults to TRUE and this allows metaMDS() to use detect whether a tranformation is 
# required in the supplied data.  It uses a commonly used standardization in ecological sciences, the Wisconsin Double Standardization, 
# after taking the square root of the data.  The autotranform arguement can be set to `FALSE` if you've already normalized your data 
# (e.g., CLR normalized).

# The phyloseq NMDS was run with the Jaccard index and stated that a convergence was not found, so we will also use the Jaccard index 
# here and alter the dimensions and max trys.  Data is not normalized, so we will keep the 'autotranform' argument at it's default setting. 

# We will use the transposed count data to ensure Samples are rows.
vegan_JaccDis <- vegdist(x=tUCDcounts, method="jaccard")

# Displaying the function on multiple lines to emphasize the multiple arguments.
UCD_jac_veganNMDS <- metaMDS(
	comm = vegan_JaccDis,
	k = 3,
	try=10,
	trymax=50
)

# View output
UCD_jac_veganNMDS

# Using 3 dimensions allows us to reach the best solution with less iterations.  Let's see the output.
# Turns out the 'Stress' statistic is lower in the previous NMDS.  This is likely due to the normalization that metaMDS() applied to 
# the count data prior to determining the Jaccard Indices.  Still, we have good/fair stress level at <0.2 (> 0.3 is considered poor). 

#####################
# Ordinations using ape package

# The `pcoa` funtion is straight forward to run; all you have to input is the symmetric square distance/dissimilarity matrix. 
# We will use Bray-Curtis Dissimilarities on non-normalized counts for this example.  
# Recall from above that the output from `pcoa` is extremely verbose, therefore, we will not print the output here.

# Load library 
library(ape)

# We will use the Bray Dissimilarities object that we calculated earlier.
UCD_bray_apePCoA <- pcoa(D=vegan_BrayDis)

#####################
## Visualizing ordinations with phyloseq

# We use ordinations to summarize high-dimensional data into a set of new components.  These components can be broken down to show the 
# variance associated with samples, often referred to as "scores".  Plotting the "score" values can sometimes show which samples are 
# more similar to each other and dissimilar from others.  More similar scores will have smaller distances between them and cluster 
# closer together, whereas dissimilar scores will have large distances between them and cluster farther away from each other. 

# The `plot_ordination` function from phyloseq is probably the most user-friendly option for plotting.  It integrates the Sample Metadata 
# in the phyloseq object to aid in the visualization of sample differences, so everything is streamlined.  The output of `plot_ordination` 
# uses the popular plotting packge, "ggplot2", as the graphical output, thus additional graphical options (e.g., point sizing, colors, etc)
 # are customizeable using the ggplot2 syntax.  
 
# `plot_ordination` requires 2 arguments at minimum: 1) a phyloseq object that created an 2) ordination from ordinate().  
# A third argument specifies whether you want to plot 'samples', 'taxa', or a 'biplot'.  The latter of these combines samples and 
# taxa information on the same plotting space.  This works well with a small set of samples and taxa, but not with the amoung of taxa 
# in our plot.  That option may work will at the phylum or class level, but not at the more granular taxonomic levels.  

# Let's see if there is any variation by diabetic group.  We can do this by changing the 'color' argument to a column name from the 
# Sample Metadata.  We can obtain these names using the `sample_names` function if we need the exact cases.

# Sample column names
sample_names(physeq=phyobj_d10p25)

# Will use the "GroupAbbrev" column
plot_ordination(physeq=phyobj_d10p25, ordination=UCD_bray_PCoA, color="GroupAbbrev")

# NMDS ordination using Bray-Curtis Dissimilarities.

# We fit NMDS with Jaccard index in the previous example
UCD_bray_NMDS <- ordinate(physeq=phyobj_d10p25, method="NMDS", distance="bray")

# Will use the "GroupAbbrev" column to layer colors on scores.
plot_ordination(physeq=phyobj_d10p25, ordination=UCD_bray_NMDS, color="GroupAbbrev")

# You can also pass arguments from metaMDS() to ordinate() as though they are the same functions.  
# Let's try without the transformation, 3 components, and alter the min/max iterations.
UCD_bray_NMDSconfig <- ordinate(physeq=phyobj_d10p25, method="NMDS", distance="bray", autotransform=FALSE, try=10, trymax=50, k=3)

# Visualize NMDS
plot_ordination(physeq=phyobj_d10p25, ordination=UCD_bray_NMDSconfig, color="GroupAbbrev")

#####################
## Visualizing ordinations with base R 

# There may be workflows that are not compatible with phyloseq, for example, it is recommended to use PCA after using a CLR transformation.

# Impute zeros with Bayesian-multiplicative treatment
UCDcounts_BMz <- cmultRepl(X=UCDcounts, method="CZM", output="p-counts")

# Taxa needs be rows and samples need to be columns
UCD_BMz_clr <- clr(xx=UCDcounts_BMz, ifwarn=FALSE)

# PCA needs samples as columns 
tUCD_BMz_clr <- t(x=UCD_BMz_clr)

# Will use the prcomp() function to perform PCA
# UCD_BMz_clr_pca <- prcomp(tUCD_BMz_clr)
UCD_BMz_clr_pca <- prcomp(x=scale(tUCD_BMz_clr), center=FALSE)

# Extract scores and make new column to join with metadata
UCD_BMz_clr_pca_scores <- UCD_BMz_clr_pca$x %>% as.data.frame()
UCD_BMz_clr_pca_scores$SampleID <- rownames(x=UCD_BMz_clr_pca_scores)

# Join data frames
UCD_BMz_clr_pca_scsDF <- inner_join(x=UCDrats_MetaData, y=UCD_BMz_clr_pca_scores, by="SampleID")

# Plot the data 
plot(
	x = UCD_BMz_clr_pca_scsDF[,"PC1"], 
	y = UCD_BMz_clr_pca_scsDF[,"PC2"], 
	pch = 21, 
	bg = c("red", "brown", "green", "blue", "pink")[factor(UCD_BMz_clr_pca_scsDF$GroupAbbrev)],
	cex = 1.5,
	xlab = "PCA Component 1",
	ylab = "PCA Component 2"
)
legend(x="topleft", legend=c("D3M", "D6M", "LSD", "PD", "RD"), pch=21, pt.bg=c("red", "brown", "green", "blue", "pink"))

#####################
# PERMANOVA

# PERMANOVA is related to an ANOVA in that it partitions variation in response to provided factor(s), but it occurs in multivariate 
# data cloud relating to a dissimilarity or distance matrix.  The statistical inferences (i.e., p-value) is determined via a 
# distribution-free permutational technique. 

# Since PERMANOVA is based on the ANOVA framework, study designs that require ANOVA testing can be adapted for PERMANOVA.  
# For example, 3-way ANOVA design with main effects and all interactions can be utilized with PERMANOVA.  In addition, covariate 
# factors can also be utilized with PERMANOVA.  

# It should be noted again, that PERMANOVA is routinely used to demonstrate whether discriminant clusters in an ordinations are 
# statistically significant, they do not always sync together.  Remember, we can only view up to 3 dimensions of an ordinations like 
# PCoA and PCA.  Sometimes 3 components only represent a small to medium amount of the variance in the data.  The PERMANOVA interprets 
# all of the variation, not just what we see in the first few components.

# We will demonstrate a 2-way PERMANOVA using the UCD-T2DM rat data. Main effects of Group and Collection will be assessed with an 
# interaction term.  

# Calculate dissimilarity matrix
phylo_BrayDis <- phyloseq::distance(physeq=phyobj_d10p25, type="bray")

# Extract sample Metadata
UCDrats_MetaData <- meta(x=phyobj_d10p25)

# Assess PERMANOVA - The exploratory variables are drawn from UCDrats_MetaData
adonis(formula=phylo_BrayDis ~ GroupAbbrev + Collection + GroupAbbrev:Collection, data=UCDrats_MetaData, permutations=500)

# The output provides an ANOVA table.  We see similar statistics that we would see from a regular ANOVA, e.g., degrees of freedom, 
# sum of squares, F-statistic, and a P-value.  Note, the F-statistic is technically a "psuedo F statistic".  These statistics can 
# be intepreted similar to a regular ANOVA.  Hence, the presence of a significant interaction terms suggest we should be assessing 
# Group differences within Collection years.  Let's do it.

# Extract the Collection Year vector
d10p25CollectionYear <- meta(x=phyobj_d10p25)$Collection

# Make the logical statement using the %in% (match) operator
CollectionYear14Logic <- d10p25CollectionYear %in% "Y14"
CollectionYear16Logic <- d10p25CollectionYear %in% "Y16"

# Use the prune_samples() function
# First argument is the logical statement, the second is the phyloseq object
phylo_year2014 <- prune_samples(samples=CollectionYear14Logic, x=phyobj_d10p25)
phylo_year2016 <- prune_samples(samples=CollectionYear16Logic, x=phyobj_d10p25)

# n of 2 in D6M group.  Will drop these samples from analysis.
phylo_year2016 <- prune_samples(samples={meta(x=phylo_year2016)$GroupAbbrev != "D6M"}, x=phylo_year2016)

# Calculate dissimilarity matrix
phylo_BrayDis2014 <- phyloseq::distance(physeq=phylo_year2014, type="bray")
phylo_BrayDis2016 <- phyloseq::distance(physeq=phylo_year2016, type="bray")

# Extract sample Metadata
UCDrats2014_MetaData <- meta(x=phylo_year2014)
UCDrats2016_MetaData <- meta(x=phylo_year2016)

# Assess PERMANOVA year 2014
adonis(formula=phylo_BrayDis2014 ~ GroupAbbrev, data=UCDrats2014_MetaData, permutations=500)
adonis(formula=phylo_BrayDis2016 ~ GroupAbbrev, data=UCDrats2016_MetaData, permutations=500)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##																													##
##               ~~~~~~~~~~~~~~~~~~~~~~       Differential Abundance        ~~~~~~~~~~~~~~~~~~~~~~ 					##
##																													##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##



# Differential Abundance

# The next step in assessing differences in microbial sequencing data is to determine which taxa differ between experimental groups.  

# Similar to beta-diversity, results can be variable depending on the normalization and statistical approach.  

# Differential abundance analysis occurs entirely outside of the phyloseq environment.  

# Many differential abundance methods are directly related to RNA-seq analysis and build an object similar to the phyloseq object, 
# but these new objects are typically unrelated to microbial data.  

#####################
## Rank analysis

# Early studies of microbial sequencing data used non-parametric rank-order statistics for group comparisons, e.g., Mann-Whitney U Test
# and Kruskal-Wallis 1-way analysis of variance by rank. This approach was adopted due to the fact that count data and proportional 
# abundances are rarely normally distributed.  However, some have critisized this approach due to the fact that the underlying data is 
# compositional and should be handled by compositional data analysis approaches.  Furthermore, the high abundance of zeros in low abundant 
# taxa inflate the number of ties in the rank and limit the power to detect differences in low abundant taxa.  

# Use the Mann-Whitney U test to determine which taxa differ between LSD and pre-diabetic (PD) UCD-T2DM rats under various normalizations.  
# Since this is 16S rRNA data and we have poor resolution at Species level, we will evaluate group differences at Genus level.

# Aggregate to Genus level
phyobj_d10p25_G <- tax_glom(physeq=phyobj_d10p25, taxrank="Genus")

# Subset phyloseq object
phyloG_LSDPD <- prune_samples(samples={meta(phyobj_d10p25_G)$GroupAbbrev %in% c("LSD", "PD")}, x=phyobj_d10p25_G)

# Extract Sample Meta Data and make new ID column with rownames
LSDPD_Gmeta <- meta(x=phyloG_LSDPD) %>% as.data.frame()
LSDPD_Gmeta$SampleID <- rownames(LSDPD_Gmeta)

# Exract OTU table 
LSDPD_Gcounts <- abundances(x=phyloG_LSDPD)

# Need to join with meta data, so samples need to be rows.
# transpose data and make new ID column with rownames
tLSDPD_Gcounts <- t(x=LSDPD_Gcounts) %>% as.data.frame()
tLSDPD_Gcounts$SampleID <- rownames(x=tLSDPD_Gcounts)

# Join count and meta data
LSDPD_Gcountmeta <- full_join(x=LSDPD_Gmeta, y=tLSDPD_Gcounts, by="SampleID")

# Reshape data to tall format 
# The '-' prefix denotes columns that should not have data reshaped into the new Taxa column
LSDPD_G_tall <- LSDPD_Gcountmeta %>% gather(key=Taxa, value=counts, -SampleID, -Rat, -GroupAbbrev, -Collection)

# Identify the unique taxa now listed in the "Taxa" column
LSDPD_G_uniqueOTU <- unique(x=LSDPD_G_tall$Taxa)

# Use sapply() to loop the Mann-Whitney U test across all unique taxa
LSDPD_MW_P <- sapply(X=LSDPD_G_uniqueOTU, function(x) {
	# Subset data within function
	DAT <- LSDPD_G_tall[LSDPD_G_tall$Taxa %in% x,]
	
	# Run Mann Whitney U on subsetted data
	MWU <- wilcox.test(formula=counts ~ GroupAbbrev, data=DAT)
	
	# Extract P.value and return
	MWU$p.value
})

# Make new data frame containing taxa identifier, p value, and FDR.
LSDPD_MW_P_df <- data.frame(
	Taxa=names(LSDPD_MW_P), # Taxa unique identifiers
	P=round(LSDPD_MW_P, digits=5), # Raw P-value
	FDR=round(p.adjust(LSDPD_MW_P, "fdr"), digits=5) # FDR corrected Pvalue, rounded to 5 digits
)

# Exract Taxa table, no extractor function, so need to manually extract.
LSDPD_Gtaxa <- tax_table(object=phyloG_LSDPD)@.Data %>% as.data.frame()

# Make Taxa column of row names for unique identifier, to join with p value data frame
LSDPD_Gtaxa$Taxa <- rownames(x=LSDPD_Gtaxa)

# Join P value and taxa data frames
LSDPD_MW_results <- full_join(x=LSDPD_Gtaxa, y=LSDPD_MW_P_df, by="Taxa")

# Subset significant results
LSDPD_MW_Sigresults <- LSDPD_MW_results[LSDPD_MW_results$FDR < 0.05,]

# Dispaly with trimmed columns
LSDPD_MW_Sigresults[,c("Phylum", "Family", "Genus", "FDR")]

# Almost 20 genera are altered between the LSD and PD rats.  However, we do not have any information regarding the dispersion of the data.  
# Let's calculate the median and IQR since we used a non-parametric method. 

# Note the linking pipes, %>%; the function sequence will terminate after the last pipe.
LSDPD_G_MEDIAN <- LSDPD_G_tall %>%

	# Define columns to get groups.  
	# In our case, we want medians by Taxa and GroupAbbrev
	group_by(Taxa, GroupAbbrev) %>%
	
	# Summarize defines the functions that we want to group by
	# Use the column name in which you want the function to run on
	summarise(MEDIAN = median(x=counts), IQR = IQR(x=counts)) 

# View object 
LSDPD_G_MEDIAN	

# We could stop here, but it would be better to organize the Mann-Whitney results with the median/IQR data.  
# With a bit of data wrangling, we can do it.

# Will use tidyverse workflow to do this.
LSDPD_G_MEDIAN_wide <- LSDPD_G_MEDIAN %>% 

	# Use gather() to combine MEDIAN and IQR columns into a single column
	gather(key=STATS, value=value, -Taxa, -GroupAbbrev) %>%
	
	# Use mutate() to make new columns that combines GroupAbbrev and STATS
	mutate(GroupStat = paste(GroupAbbrev, STATS, sep="_")) %>%
	
	# Use dplyr::select() to remove GroupAbbrev and STATS columns
	dplyr::select(-GroupAbbrev, -STATS) %>% 
	
	# Use spread() to distribute value by GroupStat
	spread(key=GroupStat, value=value) 

# Join with LSDPD_MW_results
LSDPD_MW_Final <- full_join(x=LSDPD_MW_results, y=LSDPD_G_MEDIAN_wide, by="Taxa")

# Display significant taxa with trimmed columns
LSDPD_MW_Final[LSDPD_MW_Final$FDR < 0.05, c("Genus", "FDR","LSD_MEDIAN","LSD_IQR","PD_MEDIAN","PD_IQR")]
# or, using tidyverse and pipes
LSDPD_MW_Final %>%
	filter(FDR < 0.05) %>%
	dplyr::select(Genus, FDR, LSD_MEDIAN, LSD_IQR, PD_MEDIAN, PD_IQR)


# Lets try this workflow one more time, but with CLR normalized data.  
# This workflow with CLR transformed data has been referred to as ANCOM (ANalysis of COmposition of Microbes).

# Impute zeros with Bayesian-multiplicative treatment
LSDPD_Gcounts_BMz <- cmultRepl(X=LSDPD_Gcounts, method="CZM", output="p-counts")

# Taxa needs be rows and samples need to be columns
LSDPD_G_BMz_clr <- clr(xx=LSDPD_Gcounts_BMz, ifwarn=FALSE)

# Need to join with meta data, so samples need to be rows.
# transpose data and make new ID column with rownames
tLSDPD_Gclr <- t(x=LSDPD_G_BMz_clr) %>% as.data.frame()
tLSDPD_Gclr$SampleID <- rownames(x=tLSDPD_Gclr)

# Join count and meta data
LSDPD_Gclrmeta <- full_join(x=LSDPD_Gmeta, y=tLSDPD_Gclr, by="SampleID")

# Reshape data to tall format 
# The '-' prefix denotes columns that should not have data reshaped into the new Taxa column
LSDPD_Gclr_tall <- LSDPD_Gclrmeta %>% gather(key=Taxa, value=clr, -SampleID, -Rat, -GroupAbbrev, -Collection)

# Identify the unique taxa now listed in the "Taxa" column
LSDPD_G_uniqueOTU <- unique(x=LSDPD_Gclr_tall$Taxa)

# Use sapply() to loop the Mann-Whitney U test across all unique taxa
LSDPD_clrMW_P <- sapply(X=LSDPD_G_uniqueOTU, function(x) {
	# Subset data within function
	DAT <- LSDPD_Gclr_tall[LSDPD_Gclr_tall$Taxa %in% x,]
	
	# Run Mann Whitney U on subsetted data
	MWU <- wilcox.test(formula=clr ~ GroupAbbrev, data=DAT)
	
	# Extract P.value and return
	MWU$p.value
})

# Make new data frame containing taxa identifier, p value, and FDR.
LSDPD_clrMW_P_DF <- data.frame(
	Taxa=names(LSDPD_clrMW_P), # Taxa unique identifiers
	P=round(LSDPD_clrMW_P, digits=5), # Raw pvalues
	FDR=round(p.adjust(LSDPD_clrMW_P, "fdr"), digits=5) # FDR corrected Pvalue, rounded to 5 digits
)

# Join P value and taxa data frames
# Taxa table extracted from phyloseq object in coding above
LSDPD_clrMW_results <- full_join(x=LSDPD_Gtaxa, y=LSDPD_clrMW_P_DF, by="Taxa")

# Can calculate this from LSDPD_Gclr_tall object using tidyverse
# Note the use of the pipe, %>%, the function sequence will terminate after the last pipe.
LSDPD_Gclr_MEDIANwide <- LSDPD_Gclr_tall %>%

	# Define columns to get groups.  
	# In our case, we want medians by Taxa and GroupAbbrev
	group_by(Taxa, GroupAbbrev) %>%
	
	# Summarize defines the functions that we want to group by
	# Use the column name in which you want the function to run on
	summarise(MEDIAN = median(x=clr), IQR = IQR(x=clr)) %>%

	# We can continue this pipe sequence
	# Use gather() to combine MEDIAN and IQR columns into a single column
	gather(key=STATS, value=value, -Taxa, -GroupAbbrev) %>%
	
	# Use mutate() to make new columns that combines GroupAbbrev and STATS
	mutate(GroupStat = paste(GroupAbbrev, STATS, sep="_")) %>%
	
	# Use dplyr::select() to remove GroupAbbrev and STATS columns
	dplyr::select(-GroupAbbrev, -STATS) %>% 
	
	# Use spread() to distribute value by GroupStat
	spread(key=GroupStat, value=value) 

# Join with LSDPD_MW_results
LSDPD_clrMW_Final <- full_join(x=LSDPD_clrMW_results, y=LSDPD_Gclr_MEDIANwide, by="Taxa")

# Display significant taxa with trimmed columns
LSDPD_clrMW_Final[LSDPD_clrMW_Final$FDR < 0.05, c("Genus", "FDR","LSD_MEDIAN","LSD_IQR","PD_MEDIAN","PD_IQR")]
# or, using tidyverse and pipes
LSDPD_clrMW_Final %>%
	filter(FDR < 0.05) %>%
	dplyr::select(Genus, FDR, LSD_MEDIAN, LSD_IQR, PD_MEDIAN, PD_IQR)

# How many taxa were found to be significant in both Mann Whitney workflows?

# Extract Taxa that were significant in Mann Whitney U test using raw counts
LSDPD_MW_rawcounts_sig <- LSDPD_MW_Final[LSDPD_MW_Final$FDR < 0.05, "Taxa"]

# Extract Taxa that were significant in Mann Whitney U test using CLR transformed data
LSDPD_MW_CLR_sig <- LSDPD_clrMW_Final[LSDPD_clrMW_Final$FDR < 0.05, "Taxa"]

# Match taxa IDs across both vectors and sum the number of TRUE
sum(x=LSDPD_MW_rawcounts_sig %in% LSDPD_MW_CLR_sig)

# Should be 19.  How many total taxa were significant in both analysis?
length(x=LSDPD_MW_rawcounts_sig)
length(x=LSDPD_MW_CLR_sig)

# Both should have lengths of 19, so perfect agreement.

#####################
## DESeq2

# DESeq2 is a Bioconductor based R package that was developed for RNA-seq count data analysis that using a Negative Binomial Wald Test 
# to compare dichtonomous outcomes.  

# The phyloseq_to_deseq2() phyloseq function itransfer a phyloseq object directly into a DESeq2 object.  

# As DESeq2 requires dichtonomous outcomes, we will use the LSD vs PD data from above.

# load library
library(DESeq2)

# phyloseq authors made a function that does all of the heavy lifting to get a phyloseq object into a DESeqDataSet object.
# Must be in DESeqDataSet for DESeq2 to work
# The phyloseq_to_deseq2() requires you to identify the column in the Sample Metadata that is to be compared.
LSDPD_G_deseq2 <- phyloseq_to_deseq2(physeq=phyloG_LSDPD, design = ~ GroupAbbrev)

# Run the Wald test as implemented in DESeq2
LSDPD_G_deseq2Wald <- DESeq(object=LSDPD_G_deseq2, test="Wald", fitType="parametric")

# The results() function creates a DESeq2 table of the results that we can coerce into a data frame
LSDPD_G_deseq2Waldres <- results(object=LSDPD_G_deseq2Wald, cooksCutoff = FALSE) %>% as.data.frame()

# Make a Taxa column to match the column we made in the taxa table data frame
LSDPD_G_deseq2Waldres$Taxa <- rownames(x=LSDPD_G_deseq2Waldres)

# Join DESeq2 results and taxa data frames
# Taxa table extracted from phyloseq object in coding above
LSDPD_DESeq2_results <- full_join(x=LSDPD_Gtaxa, y=LSDPD_G_deseq2Waldres, by="Taxa")

# View abbreviated table
LSDPD_DESeq2_results[LSDPD_DESeq2_results$padj < 0.05,c("Genus","baseMean","log2FoldChange","padj")]

# How many taxa were significant in the DESeq2 analysis?

# Extract Taxa that were significant in DESeq2 results
LSDPD_MW_DESeq2_sig <- LSDPD_DESeq2_results[LSDPD_DESeq2_results$padj < 0.05,"Taxa"]

# Total significant taxa
length(x=LSDPD_MW_DESeq2_sig)

# How many match with raw count MWU approach?
sum(x=LSDPD_MW_DESeq2_sig %in% LSDPD_MW_rawcounts_sig)

# Only 14.  So not exact agreement with earlier test.

#####################
## ALDEx2

# ALDEx stands for Anova Like Differential  EXpression analysis.  It uses a Dirichlet-multinomial model to infer abundances from 
# counts data and then calculates sample variation using standard statistical tools.  Unlike most of the methods we have looked at, 
# ALDEx2 takes into consideration the compositional nature of the microbial sequencing count data.  

# My understanding of ALDEx2 is that it handles the zero counts in a unique way.  Essentially it estimates many possible imputations 
# for zeros and then tests each of these permutations for significance, and then takes the average value for each significance test 
# as the final determination.  So, very similar to the CLR-Mann Whitney U example, but instead of doing this a single time, it runs 
# the test many times with different combinations of the zero-imputations. 

# As the CLR transformation is supposed to tranform the data into an Euclidean space, more traditional statistical approaches, 
# e.g., _t_-test, Mann Whitney U test, Kruskal Wallis test, and general linear models can be conducted with this approach.  Therefore, 
# covariates can be included as well.

# Load library
library(ALDEx2)

# Calculate Monte Carlo samples of the Dirichlet distribution for each sample.
# Each calculation uses CLR transformation
LSDPD_GaldexMCcalc <- aldex.clr(reads=LSDPD_Gcounts, conds=LSDPD_Gmeta$GroupAbbrev, mc.samples=130, denom="all", verbose=FALSE)

# Assess t-test among all Monte Carlo samples and then make Taxa column to merge with Taxonomy table
LSDPD_Galdex_ttest <- aldex.ttest(clr=LSDPD_GaldexMCcalc, paired.test=FALSE)
LSDPD_Galdex_ttest$Taxa <- rownames(x=LSDPD_Galdex_ttest)

# Merge with taxonomy table
LSDPD_Galdex_Final <- full_join(x=LSDPD_Gtaxa, y=LSDPD_Galdex_ttest, by="Taxa")

# Assess effect sizes for all Monte Carlo samples and then make Taxa column to merge with Taxonomy table
LSDPD_Galdex_esize <- aldex.effect(clr=LSDPD_GaldexMCcalc, verbose=FALSE)
LSDPD_Galdex_esize$Taxa <- rownames(x=LSDPD_Galdex_esize)

# Merge with LSDPD_Galdex_Final object
LSDPD_Galdex_Final <- full_join(x=LSDPD_Galdex_Final, y=LSDPD_Galdex_esize, by="Taxa")

# View abbreviated table. "BH" refers to Benjamini and Hochberg FDR correction.
LSDPD_Galdex_Final[LSDPD_Galdex_Final$we.eBH < 0.05,c("Genus", "we.eBH", "effect")]

# Extract Taxa that were significant in ALDEx2 results
LSDPD_MW_ALDEx2_sig <- LSDPD_Galdex_Final[LSDPD_Galdex_Final$we.eBH < 0.05,"Taxa"]

# Total significant taxa
length(x=LSDPD_MW_ALDEx2_sig)

# How many match with raw count MWU approach?
sum(x=LSDPD_MW_ALDEx2_sig %in% LSDPD_MW_rawcounts_sig)

# All of the significant taxa from the ALDEx were significant in the raw count analysis with MWU. 

# The authors of ALDEx2 suggest that differences with an absolute effect size > 1 have stronger evidence for true differences. 
LSDPD_Galdex_Final[LSDPD_Galdex_Final$we.eBH < 0.05 & abs(LSDPD_Galdex_Final$effect) > 1,c("Genus", "we.eBH", "effect")]
 
# Based on this criteria, only 7 genera are significant.  




