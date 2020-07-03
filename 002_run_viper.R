# libraries ---------------------------------------------------------------
library(tidyverse)
library(viper)
library(bcellViper)

# load in the files -------------------------------------------------------
# the expression set is from 001_create_eset.R
# chek the object
exampleSet
# check the pheno data associated to the object
pData(exampleSet)
# check the expressio values of the sample
exampleSet@assayData$exprs

# load the regulatory network file (in this case we use a precompiled version provided from the lab)
# we should build it usign endothelila cells array (see ARACNe documentation)
# load("../data_senescence_ECFC/bcell_u133p2_cleaner_mas5_tfregulon.rda")
# load("../data_senescence_ECFC/bcell_u95av2_cleaner_mas5_tfregulon.rda")
# load("../data_senescence_ECFC/gbm_u133a_cleaner_mas5_tfregulon.rda")
# load("../data_senescence_ECFC/brca_tcga_rnaseq851_signalomeregulon.rda")
# 
# names(regul)
# sum(is.na(names(regul)))
# names(regulon)
# sum(is.na(names(regulon)))
# regul
load(system.file("data","bcellViper.rda",package = "bcellViper"))
regulon
str_subset(names(regulon),pattern = "IRF7|STAT|NFKB|CEBPB")
# adjfile <- system.file("aracne", "bcellaracne.adj", package = "bcellViper")
# # check the structure of the adjfile
# read_tsv(adjfile,col_names = F)
# 
# # create the regulon object from the adjfile and the expression file 
# regul <- aracne2regulon(adjfile, exampleSet, verbose = FALSE)
# # see the content of regul object
# print(regul)
# names(regul)
# -------------------------------------------------------------------------
# Generating the gene expression signatures

signature_sen <- rowTtest(exampleSet, "subtype", c("LP", "X","ET"), "EP")
# signature_sen <- rowTtest(exampleSet, "subtype", c("X"), "EP")
# signature_sen <- rowTtest(exampleSet, "subtype", c("ET"), "EP")
# signature_sen <- rowTtest(exampleSet, "subtype", c("LP"), "EP")
# estimate z-score values for the GES:
signature_sen <- (qnorm(signature_sen$p.value/2, lower.tail = FALSE) * sign(signature_sen$statistic))[, 1]

# Build NULL model by sample permutations
nullmodel <- ttestNull(exampleSet, "subtype", c("LP", "X","ET"), "EP", per = 1000,repos = TRUE, verbose = FALSE)
# nullmodel <- ttestNull(exampleSet, "subtype", c("X"), "EP", per = 1000,repos = TRUE, verbose = FALSE)
# nullmodel <- ttestNull(exampleSet, "subtype", c("ET"), "EP", per = 1000,repos = TRUE, verbose = FALSE)
# nullmodel <- ttestNull(exampleSet, "subtype", c("LP"), "EP", per = 1000,repos = TRUE, verbose = FALSE)

# run msVIPER
mrs_sen <- msviper(signature_sen, regulon, nullmodel, verbose = FALSE)
summary(mrs_sen,mrs = 100)

# plot the result
plot(mrs_sen, cex = .7,mrs = 30)

mrs_sen$regulon$IRF7

# Leading-edge analysis
mrs_sen <- ledge(mrs_sen)
summary(mrs_sen)