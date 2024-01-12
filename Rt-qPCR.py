# Imports
import pandas as pd
import numpy as np
import statistics
import matplotlib.pyplot as plt
import scipy.stats as stats

# Read original file including path
original_file = pd.read_excel(r"")

# Create main dataframe
df = pd.DataFrame(original_file, columns= ['Group# - TreatmentType', 'Target', 'Cq', 'Avg/group'])

# Delete Rows with undetermined CT
df = df[~df['CT'].isin(['Undetermined'])]

### Create dataframes ###
# Gene of Interest treated
GOI = pd.DataFrame(df, columns= ['Sample Name', 'CT', 'Ct Mean'])

# Housekeeping gene treated
HKG = pd.DataFrame(df, columns= ['Sample Name', 'CT', 'Ct Mean'])

# Gene of interest and housekeeping gene not treated
GOIW = pd.DataFrame(df, columns= ['Sample Name', 'CT', 'Ct Mean'])
HKGW = pd.DataFrame(df, columns= ['Sample Name', 'CT', 'Ct Mean'])

# Organize values for each dataframe from excel results
HKG = df[df['Sample Name'].isin(['preh c', '10m h c', '15m h c'])]
GOI = df[df['Sample Name'].isin(['preH G', '10m h g', '15m h g'])]
HKGW = df[df['Sample Name'].isin(['prew c', '10m w c', '15m w c'])]
GOIW = df[df['Sample Name'].isin(['prew g', '10m w g', '15m w g'])]

# Create control as an integer
GOI["GOI_int"] = GOI["CT"].astype(np.float64)
HKG["HKG_int"] = HKG["CT"].astype(np.float64)
GOIW["GOIW_int"] = GOIW["CT"].astype(np.float64)
HKGW["HKGW_int"] = HKGW["CT"].astype(np.float64)

# Create delta CT(CQ) column
GOI["Delta Ct"] = GOI.iloc[0, 2] - GOI["GOI_int"]
HKG["Delta Ct"] = HKG.iloc[0, 2] - HKG["HKG_int"]
GOIW["Delta Ct"] = GOIW.iloc[0,2] - GOIW["GOIW_int"]
HKGW["Delta Ct"] = HKGW.iloc[0,2] - HKGW["HKGW_int"]

# Create variable for primer efficiencies
PE_GOI = 1.873
PE_HKG = 1.864

### Gene expression ratio ###
# PE_GOI^Delta Ct GOI / PE_HKG^Delta Ct HKG
# Treated
Step_1 = np.power((PE_GOI), GOI["Delta Ct"].reset_index(drop=True))
Step_2 = np.power((PE_HKG), HKG["Delta Ct"].reset_index(drop=True))
GER = (Step_1/Step_2)

# Not treated
Step_1b = np.power((PE_GOI), GOIW["Delta Ct"].reset_index(drop=True))
Step_2b = np.power((PE_HKG), HKGW["Delta Ct"].reset_index(drop=True))
GERb = (Step_1b/Step_2b)

# Concatenate dataFrames
HKG.reset_index(drop=True, inplace=True)
GOI.reset_index(drop=True, inplace=True)
GOI_HKG = pd.concat([GOI, HKG], axis=1)

GOIW.reset_index(drop=True, inplace=True)
HKGW.reset_index(drop=True, inplace=True)
GOIW_HKGW = pd.concat([GOIW, HKGW], axis=1)

# Create Gene Expression Ratio Column
GOI_HKG["Gene Expression Ratio"] = Step_1/Step_2

GOIW_HKGW["Gene Expression Ratio"] = Step_1b/Step_2b

# Rename the columns
GOI_HKG.columns = ["Sample Name", "CT", "Ct Mean", "GOI_int", "Delta Ct", "Sample Name 2", "CT 2", "Ct Mean 2", "HKG_int", "Delta Ct 2", "Gene Expression Ratio"]

GOIW_HKGW.columns = ["Sample Name", "CT", "Ct Mean", "GOIW_int", "Delta Ct", "Sample Name 2", "CT 2", "Ct Mean 2", "HKGW_int", "Delta Ct 2", "Gene Expression Ratio"]

### Verify results with statistics ###
# Create Standard Deviation Column
Standard_Deviation = GOI_HKG.groupby("Sample Name", sort=False)["Gene Expression Ratio"].std()
Standard_DeviationW = GOIW_HKGW.groupby("Sample Name", sort=False)["Gene Expression Ratio"].std()

# Create the Standard error column
for i in range(0, 6):
    GOI_HKG.at[i, 'Standard Deviation'] = Standard_Deviation[GOI_HKG.at[i, 'Sample Name']]
    GOI_HKG.at[i, 'Standard Error'] = (GOI_HKG.at[i, 'Standard Deviation']/np.sqrt(2))
    GOIW_HKGW.at[i, 'Standard Deviation'] = Standard_DeviationW[GOIW_HKGW.at[i, 'Sample Name']]
    GOIW_HKGW.at[i, 'Standard Error'] = (GOIW_HKGW.at[i, 'Standard Deviation']/np.sqrt(2))

# T-test
# treated vs pre-treated
group = df.groupby('Sample Name')['Gene Expression Ratio']
control_expression = group.get_group('control')
treated_expression = group.get_group('treated')

t_statistic, p_value = stats.ttest_ind(control_expression, treated_expression)

### Print statements and return results to excel ###
print("T-Statistic: ", t_statistic)
print("P-Value: ", p_value)
GOIW_HKGW.to_excel("Results.xlsx")