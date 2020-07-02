library(viper)
library(bcellViper)
library(tidyverse)

# -------------------------------------------------------------------------
# Generating the regulon object
# As described under ‘Overview of VIPER’ (section 1), msVIPER and VIPER require a gene expression signature and an appropriate cell context-specific regulatory network. This regulatory network is provided in the format of a class regulon object. Regulon objects can be generated from networks reverse engineered with the ARACNe algorithm [1]. This is performed by the function aracne2regulon, which takes two arguments as input: the ARACNe output .adj file, and the expression data-set used by ARACNe to reverse engineer the network. As an example, the package bcellViper provides a subset of the ARACNe output file containing the network for 20 TF regulators (bcellaracne.adj file). For convenience, the full network is also provided, as a regulon class object, together with the gene expression data used to reverse engineer it contained in an EpressionSet object. The B-cell expression data contains 211 samples representing several normal and tumor human B-cell phenotypes profiled on Affymetrix H-GU95Av2 (Gene Expression Omnibus series GSE2350)[1]. The provided dataset was generated from custom probe-clusters obtained by the the cleaner algorithm[3] and MAS5[2] normalization. The following lines are an example for the use of aracne2regulon function to generate the regulon object from the ARACNe output data and the expression data:

data(bcellViper, package="bcellViper")
dset
pData(dset)
dset@assayData$exprs
table(pData(dset)$description)
adjfile <- system.file("aracne", "bcellaracne.adj", package = "bcellViper")
read_tsv(adjfile,col_names = F)
regul <- aracne2regulon(adjfile, dset, verbose = FALSE)
print(regul)


# -------------------------------------------------------------------------
# 6 Master Regulator Analysis performed by msVIPER
# To illustrate this section, we analyze part of the expression data from [1], consitent on 5 na¨ıve human B-cell, 5 memory B-cell, 5 centroblast and 5 centrocyte B-cell samples profiled on Affymetrix H-GU95Av2 gene arrays. The complete dataset is available from Gene Expression Omnibus (GSE2350), and here for convenience, we have included the ‘cleaner’[3] processed and MAS5[2] normalized samples in the bcellViper package.

# 6.1 Generating the gene expression signatures
# Lets assume that we are interested in identifying transcriptional regulators associated with the Germinal Center (GC) reaction. GC are the peripheral lymphoid organs where antigen-driven somatic hypermutation of the genes encoding the immunoglobulin variable region occurs, and are the main source of memory B cells and plasma cells that produce high-affinity antibodies. Centroblast and centrocyte are the main Bcell phenotypes present in the GCs, they are derived from antigen-stimulated peripheral blood B-cells, and represent the most proliferative cellular physiologic phenotypes of the adult human body. Thus, we can obtain a gene expression signature for the GC formation by comparing GC (centroblasts and centrocytes) against na¨ıve B-cells. The ‘ExpressionSet’ object available from the bcellViper data package contains 5 centroblast samples (CB), 5 centrocyte samples (CC) and 5 na¨ıve peripheral blood B-cell samples (N). The viper package includes the function rowTtest that efficiently performs Student’s t-test for each row of a dataset. The rowTtest function conveniently takes an ‘ExpressionSet’ object as argument and produces a list object containing the Student’s t-statistic (statistic) and the test’s p-value (p.value), that by default is estimated by a 2-tail test.

signature <- rowTtest(dset, "description", c("CB", "CC"), "N")

# It can also take two matrixes as arguments, the first one containing the ‘test’ samples and the second the ‘reference’ samples. While we could define the Gene Expression Signature (GES) by using the t-statistic, to be consistent with the z-score based null model for msVIPER (see section 6.2), we will estimate z-score values for the GES:

signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * sign(signature$statistic))[, 1]

# 6.2 NULL model by sample permutations
# A uniform distribution of the targets on the GES is not a good prior for msVIPER. Given the high degree of co-regulation in transcriptional networks, the assumption of statistical independence of gene expression is unrealistic an can potentially lead to p-value underestimates. To account for the correlation structure between genes, we define a null model for msVIPER by using a set of signatures obtained after permuting the samples at random. The function ttestNull performs such process by shuffling the samples among the ‘test’ and ‘reference’ sets, according to the re-sampling mode and number of permutations indicated by the parameters repos and per, respectively.
nullmodel <- ttestNull(dset, "description", c("CB", "CC"), "N", per = 1000,repos = TRUE, verbose = FALSE)

# As output, the ttestNull function produces a numerical matrix of z-scores, with genes/probes in rows and permutation iterations in columns, than can be used as null model for the msVIPER analysis.

# 6.3 msVIPER
# The last element required by msVIPER that we are still missing is an apropriate cell context-specific regulatory network. We have included a B-cell regulatory network in the bcellViper package, and additional networks described in [12] for human B-cell, glioma and breast carcinoma can be obtained from figshare (Table 1).
regulon

# The msVIPER analysis is performed by the msVIPER function. It requires a GES, regulon object and null model as arguments, and produces an object of class ‘msVIPER’, containing the GES, regulon and estimated enrichment, including the Normalized Enrichment Score (NES) and p-value, as output.

mrs <- msviper(signature, regulon, nullmodel, verbose = FALSE)

# The reults can be summarized by the generic function summary, which takes the msviper object and either the number of top regulators to report or a specific set of regulators to list. The default for this parameter is the top 10 master regulators (MRs).
summary(mrs)

# A graphics representation of the results (msVIPER plot) can be obtained by the generic function plot (shown in Fig. 1). It takes the msviper object and either, the number of top differentially active regulators, or the names of the regulators to include in the plot as arguments. The default behavior is to plot the top 10 most differentially active MRs.
plot(mrs, cex = .7)

# 6.3.1 Leading-edge analysis
# msVIPER infers the relative activity of a regulatory gene based on the enrichment of its most closelyregulated targets on a given GES, but does not identify which are the target genes enriched in the GES. Subramanian et al. [14] proposed a method called leading-edge analysis to identify the genes driving the enrichment of a gene-set on a GES based on Gene Set Enrichment Analysis (GSEA). We implemented the leading-edge analysis in the ledge function of the viper package. The function only has a ‘msviper’ class object as argument and generates an updated ‘msviper’ object that now includes a ‘ledge’ slot.
mrs <- ledge(mrs)
summary(mrs)

# 7 Beyond msVIPER --------------------------------------------------------
# 7.1 Bootstrap msVIPER
# The effect of outlier samples on the gene expression signature can be reduced by the use of resampling techniques. msVIPER is capable of performing the analysis with bootstrap if a matrix of bootstraped signatures, instead of a vector, is given as signature argument. We implemened the function bootstrapTtest in the viper package to generate this kind of bootstraped GES matrixes from the ‘test’ and ‘reference’ datasets. The function produces 100 bootstrap interactions by default.
signature <- bootstrapTtest(dset, "description", c("CB", "CC"), "N", verbose = FALSE)
mrs <- msviper(signature, regulon, nullmodel, verbose = FALSE)

# By default, msviper integrates the regulator activity results across all bootstraped iteration using the average, but this can be easily modified to use the median or mode values by the bootstrapmsviper function:
mrs <- bootstrapmsviper(mrs, "mode")

# Bootstraped msviper results can be displayed in the same way as non-bootstraped results (Fig. 2):
plot(mrs, cex = .7)
