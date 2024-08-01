#document:https://rdrr.io/bioc/FELLA/f/vignettes/quickstart.Rmd
# FELLA needs:

# The KEGG graph and other complementary data files. This is stored in a unique FELLA.DATA S4 object.
# A list of affected metabolites (KEGG compounds). This is stored in a unique FELLA.USER S4 object, 
# along with user analyses.

# ################################################################
# ############I. Loading the KEGG data
# ################################################################
library(FELLA)
# we make use of sample data that contains small subgraph of FELLA's KEGG graph (mid 2017 KEGG release)
data("FELLA.sample")
# class(FELLA.sample)
# show(FELLA.sample)

# ################################################################
# ############II. Loading the metabolomics summary data
# ################################################################
#The second block of necessary data is a list of affected metabolites, which shoud be specified as KEGG compound IDs. 
#Provided is a list of hypothetical affected metabolites belonging to the graph, to which some decoys that do not map to 
#the graph are added.
data("input.sample")
input.full <- c(input.sample, paste0("intruder", 1:10))
# show(input.full)


myAnalysis <- defineCompounds(
    compounds = input.full, 
    data = FELLA.sample)

# ################################################################
# ############III. Enriching the data
# ################################################################
#Once the FELLA.DATA and the FELLA.USER with the affected metabolites are ready, the data can be easily enriched.
myAnalysis <- enrich(
    compounds = input.full, 
    method = "diffusion",  #Enrichment methods: "hypergeom", "diffusion", or "pagerank"
    approx = "normality",  #Statistical approximations: "normality" or "simulation"
    data = FELLA.sample)
show(myAnalysis)

# ################################################################
# ############IV. Exporting the results
# ################################################################
myTable <- generateResultsTable(
    object = myAnalysis, 
    method = "diffusion", 
    threshold = 0.1, 
    data = FELLA.sample)

knitr::kable(head(myTable, 20))