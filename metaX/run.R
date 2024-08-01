#document: https://github.com/wenbostar/metaX/blob/master/vignettes/metaX.pdf
library("metaX")

# Input to the tool is a list of metabolites converted into HMDB Compound IDs
HMDB_metabolites <- c("HMDB00060","HMDB00056","HMDB00064")

# Perform pathway analysis
res <- pathwayAnalysis(id=HMDB_metabolites,
                    id.type="hmdb",
                    outfile="./output/pathway.csv")