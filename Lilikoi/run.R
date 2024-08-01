#document: https://cran.r-project.org/web/packages/lilikoi/vignettes/Vignette.html 

library(lilikoi)

###################################################
############I. LOAD EXAMPLE DATA
###################################################
dt <- lilikoi.Loaddata(file=system.file("extdata", "plasma_breast_cancer.csv", package = "lilikoi"))
Metadata <- dt$Metadata
dataSet <- dt$dataSet


###################################################
############II. ID Mapping
###################################################
#Transform the metabolite names to the HMDB ids using Lilikoi MetaTOpathway function
convertResults=lilikoi.MetaTOpathway(q.type='name')
Metabolite_pathway_table = convertResults$table

# ###################################################
# ############III. Functional analysis
# ###################################################
# # Calculate the Pathway Dysregulation score
PDSmatrix=lilikoi.PDSfun(Metabolite_pathway_table)

# #Select the most signficant pathway related to phenotype.
selected_Pathways_Weka= lilikoi.featuresSelection(PDSmatrix,
                                                    threshold= 0.50, #select the top pathways
                                                    method="gain") # "info" (information gain) or "gain" (gain ratio)
print(selected_Pathways_Weka)