import sspa
import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

#document: https://pypi.org/project/sspa/0.2.1/


###################################################
############I. LOAD EXAMPLE DATA
###################################################
covid_data = sspa.load_example_data(omicstype="metabolomics", processed=False)
# Statistics of Covid data
# Group
# Healthy Donor     133
# COVID19           133
# Name: count, dtype: int64

###################################################
############II. PROCESS DATA
###################################################
# Keep only metabolites (exclude metadata columns)
covid_values = covid_data.iloc[:, :-2]

# Remove features with too many na values
data_filt = covid_values.loc[:, covid_values.isin([' ', np.nan, 0]).mean() < 0.5]

# Impute using the median
imputed_mat = data_filt.fillna(data_filt.median())
# Log transform the data
log2_mat = np.log2(imputed_mat)
# Standardise the data
processed_data = pd.DataFrame(StandardScaler().fit_transform(log2_mat), columns=imputed_mat.columns, index=imputed_mat.index)

###################################################
############III. LOAD PATHWAY DATA FROM KEGG
###################################################
# download pathway data from KEGG via API 
# NOTE: Optional; do this process just to keep the pathway information from KEGG up-to-date. I downloaded it and stored it for future use.
base_dir = '/code'
# kegg_human_latest = sspa.process_kegg("hsa", download_latest=True, filepath= base_dir )

# Download can take hours depending on the users' internet
kegg_human_ptw_db = sspa.process_gmt( os.path.join(base_dir, "KEGG_hsa_pathways_compounds_R110.gmt") )

###################################################
############IV. ID Mapping
###################################################
# Create a table to convert metabolite names into KEGG Compound IDs  
# NOTE: This process is very slow, so I did it and stored the resulting table for future use. 
# If users want to run this tool with their data, they need to do this step!

# compound_names = processed_data.columns.tolist()
# conversion_table = sspa.identifier_conversion(input_type="name", compound_list=compound_names)

conversion_table = pd.read_csv( os.path.join(base_dir, "met_id_mapping.csv") )

# Convert metabolite names to KEGG Compound IDs (e.g., C00315, C02918,...)
processed_data_mapped = sspa.map_identifiers(conversion_table, output_id_type="KEGG", matrix=processed_data)
processed_data_mapped

###################################################
############V. Pathway Analysis
###################################################
# Perform pahtway analysis using kPCA
kpca_scores = sspa.sspa_KPCA(kegg_human_ptw_db, min_entity=3, random_state=1).fit_transform(processed_data_mapped)
# ** Use t-test to get the pvalue for each pathway **

# Perform pahtway analysis using ssClustPA
ssclustpa_res = sspa.sspa_ssClustPA(kegg_human_ptw_db).fit_transform(processed_data_mapped)
# ** Use t-test to get the pvalue for each pathway **