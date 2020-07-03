# libraries ---------------------------------------------------------------
library(tidyverse)
library(viper)
library(bcellViper)

# load the data and wrangling ---------------------------------------------
# these are the data already normalized (ExpressionSet object)
df <- read_tsv("../data_senescence_ECFC/mRNA_Expression_tidy.txt")

# load the lut to convert the probenames into entrez genenames
# load the lut and conver the probenames into genenames
lut <- read_tsv(file = "../data_senescence_ECFC/agilent_full_fix.txt")

# clean the dataset
df2 <- inner_join(df,lut,by=c("NAME"="ProbeName"))%>%
  # select only the column of interest
  select(GeneSymbol,h_192_EP:h_B2_IF) %>%
  # select(EntrezID,h_192_EP:h_B2_IF)%>%
  # select only the non identifyers
  filter(!is.na(GeneSymbol))%>%
  # filter(!is.na(EntrezID))%>%
  # collapse the same genes by average them
  group_by(GeneSymbol)%>%
  # group_by(EntrezID)%>%
  summarise_at(colnames(.)[-1], mean, na.rm = TRUE) %>%
  rename(rowname = GeneSymbol)
  # rename(rowname = EntrezID)


head(df2,n = 20)
dim(df2)
dim(lut)
# build eset --------------------------------------------------------------

# -------------------------------------------------------------------------
# 01 create an Expression set object
#
# exp should be a matrix of normalized expression
exp <- df2 %>%
  # add the probename to the rows
  column_to_rownames() %>%
  # as to be a matrix to be accepted
  as.matrix()

# -------------------------------------------------------------------------
# 02 pheno data 
# is a complex object build from a df with the organization of the samples (the classes of the sample). notice that the rownames has to be the same of the colnames of the expression set
samples <- colnames(exp)

# build the data frame from the samples
df_pheno <- data.frame(rowname = samples) %>%
  mutate(subtype = str_replace(rowname,pattern = "h_.*_",replacement = ""),
         clone = str_replace_all(rowname,pattern = "h_|_EP|_X|_LP|_ET|_IF",replacement = "")) %>%
  column_to_rownames()

# create the pheno object
pheno <- new("AnnotatedDataFrame",data=df_pheno)
head(pheno)

# -------------------------------------------------------------------------
# 03 experiment data 
# they are some metadata abou the run can be set as generic
expData <- new("MIAME",
               name="test",
               lab="test",
               contact="test",
               title="test",
               abstract="test",
               url="test",
               other=list(
                 notes="test"
               ))

class(expData)

# -------------------------------------------------------------------------
# 04 the annotation object
# it holds the info about the probe set, is just a character string can be omitted

# -------------------------------------------------------------------------
# 05 build a new expression set from its component
# annotate everything with the expression set and use org.Hs.eg.db as annotation to do that I would have needed to use the entrez genename as rows
# with annotation info
# exampleSet <- ExpressionSet(assayData = exp,
#                             phenoData = pheno,
#                             experimentData = expData,
#                             annotation = "org.Hs.eg.db")

# no annotation info
exampleSet <- ExpressionSet(assayData = exp,
                            phenoData = pheno,
                            experimentData = expData)

exampleSet
# confirm the object is an eset
class(exampleSet)
