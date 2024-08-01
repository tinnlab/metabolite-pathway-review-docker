#Document: chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.bioconductor.org/packages//2.12/bioc/vignettes/PAPi/inst/doc/PAPiPackage.pdf

###################################################
### I. load necessary R libraries and load example data
###################################################
library(PAPi)
library(svDialogs)
data(papiData)
print(papiData)

###################################################
### II. Download the KEGG pathway database (OPTIONAL)
###################################################
# buildDatabase(save = F, saveAs = "KEGG")

# NOTE: Downloading KEGG is very slow, maybe taking up to 10 hours
# We will choose to run PAPi online

###################################################
### III. Run PAPi
###################################################
papiResults <- papi(papiData, save = FALSE, offline = FALSE)
print(head(papiResults))

###################################################
### IV. Apply ANOVA or t-test to the results produced by PAPi
###################################################
head(papiResults)
ApplyingHtest <- papiHtest(
	papiResults,
	save = FALSE,
	StatTest = "T"
)
head(ApplyingHtest)

