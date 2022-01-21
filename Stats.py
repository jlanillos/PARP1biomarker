# Obtain RNAseq data from the IMmotion151:
#conda activate py36; pyega3 fetch EGAF00004859518
import pandas as pd
from lifelines import KaplanMeierFitter #pip install lifelines
import matplotlib.pyplot as plt
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test
import numpy as np
from lifelines import CoxPHFitter

c = pd.read_excel('mmc2(4).xlsx') # Table S1 from Motzer et al, Cancer Cell 38, 803–817, December 14, 2020
df = pd.read_csv('EGAF00004859518/IMmotion151.expression.data.TPM.anon.20201106.csv',sep=',') # Transcript counts matrix from EGAF00004859518
parp = df.loc[df['symbol'] == 'PARP1'] # Find PARP1 expression
parp.reset_index(inplace=True)
d = dict(zip(list(parp.T.index), list(parp.T[0].values)))

c['PARP1'] = c['RNASEQ_SAMPLE_ID'].map(d) # Assign PARP1 expression values corresponding to each tumor/case
# Calculate mean and median Z-scores for PARP1
c['zPARP1'] = (c['PARP1'] - c['PARP1'].mean()) / (c['PARP1'].std())
c['zmedianPARP1'] = (c['PARP1'] - c['PARP1'].median()) / (c['PARP1'].std())

# cf will be the definite dataframe that contains both clinical and expression data
cf = c.loc[~c['zPARP1'].isnull()]
cf['zgroup'] = ''
cf['zmediangroup'] = ''
cf.loc[cf['zPARP1'] < 0, 'zgroup'] = 0 # low PARP1, group 0, according to mean-Zscore
cf.loc[cf['zPARP1'] > 0, 'zgroup'] = 1 # high PARP1, group 1
cf.loc[cf['zmedianPARP1'] < 0, 'zmediangroup'] = 0 # low PARP1, group 0, according to mean-Zscore
cf.loc[cf['zmedianPARP1'] > 0, 'zmediangroup'] = 1 # high PARP1, group 1, according to mean-Zscore

# Keep only samples with mutational data (there are samples that do not contain their mutational status, they were not included in NGS gene panel)
cf = cf.loc[~cf['FMI_SAMPLE_ID'].isnull()]

# Binarize PBRM1 mutational status and prepare PARP1+PBRM1 info for each sample
cf['PBRM1_mut'] = np.where(cf['FMI_PBRM1'].isin(['loss','rearrangement', 'short-variant']), 1, 0)
cf['PBRM1_mut_term'] = np.where(cf['FMI_PBRM1'].isin(['loss','rearrangement', 'short-variant']), 'PBRM1-mut', 'PBRM1-wt')
cf['PARP1_zgroup_term'] = np.where(cf['zgroup']  == 1, 'PARP1-High', 'PARP1-Low')
cf['PARP1_PBRM1'] = cf['zgroup'].astype(str) + '_' + cf['PBRM1_mut'].astype(str)
cf['PARP1_PBRM1_term'] = cf['PARP1_zgroup_term'].astype(str) + '/' + cf['PBRM1_mut_term'].astype(str)
# MSKCC prognostic group, assign the term "favorable" to "low" values, assign "Int_Poor" otherwise
cf['MSKCCPrognosticGroup'] = np.where(cf['MSKCC_RISK_SCORE'] == 'Low', 'Favorable', 'Int_Poor')
# Keep samples annotated with ccRCC histology. Filter out the rest of them:
cf = cf.loc[cf['HISTOLOGY_SARCOMATOID'].isin(['ccRCC_Sarc', 'ccRCC_nonSarc'])].copy()

#############################################################################
# Save the initial dataframe into a safe copy (cf_aux)
cf_aux = cf.copy()


# Create the figure
# stats dataframe will serve to store output results from statistical analysis (Cox-regressions)
stats = pd.DataFrame(columns=['ARM', 'comp', 'HR', '95CI', 'p'])

################################################################
################# Atezolizumab+Bevacizumab #####################
################################################################
cf = cf_aux.copy()
arm = 'Atezolizumab+Bevacizumab'
cf = cf.loc[cf['ARM'] == arm]
#### Mean Z-score (PARP1 HIGH/ LOW) - Censored and Not Censored
kmf = KaplanMeierFitter()
T = cf['PFS_MONTHS'].loc[cf['zgroup'] == 1]
E = cf['PFS_CENSOR'].loc[cf['zgroup'] == 1]
kmf.fit(T, event_observed=E, label="PARP1-high")
kmfaux = kmf

kmf = KaplanMeierFitter()
T1 = cf['PFS_MONTHS'].loc[cf['zgroup'] == 0]
E1 = cf['PFS_CENSOR'].loc[cf['zgroup'] == 0]
kmf.fit(T1, event_observed=E1, label="PARP1-low")

results=logrank_test(T1,T,event_observed_A=E1, event_observed_B=E)
offset = 0.08
if results.p_value < 0.05:
    print('SUTEN: PARP1-low vs PARP1-high: ' + str(results.p_value))
    #ax[0,0].text(0, 0+ offset, 'PARP1-low vs PARP1-high: ' + str("{0:.3f}".format(results.p_value)), fontsize=fs)
    offset = offset + 0.05


# CoxPHFitter - PARP1 Univariate
cph = CoxPHFitter()
df_r = cf[['zgroup','PFS_MONTHS','PFS_CENSOR']]
df_dummy = pd.get_dummies(df_r, drop_first=True)
colcov = 'zgroup_1'
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR')
#cph.print_summary()
comp = 'PARP1-univ'
hr = np.round(cph.summary.T[colcov]['exp(coef)'], decimals=2)
ci = str(np.round(cph.summary.T[colcov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[colcov]['exp(coef) upper 95%'], decimals=2))
pval = cph.summary.T[colcov]['p']
stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)

offset = 0.08
pval = np.round(cph.summary.T[colcov]['p'], decimals = 4)
fs_text = 10


# CoxPHFitter - MSKCC Univariate
cph = CoxPHFitter()
df_r = cf[['PFS_MONTHS','PFS_CENSOR','MSKCCPrognosticGroup']]
df_dummy = pd.get_dummies(df_r, drop_first=True)
colcov = 'MSKCCPrognosticGroup_Int_Poor'
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR')
comp = 'MSKCC-univ'
hr = np.round(cph.summary.T[colcov]['exp(coef)'], decimals=2)
ci = str(np.round(cph.summary.T[colcov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[colcov]['exp(coef) upper 95%'], decimals=2))
pval = cph.summary.T[colcov]['p']
stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)


# CoxPHFitter - PARP1, MSKCC Multivariate
cph = CoxPHFitter()
df_r = cf[['PFS_MONTHS','PFS_CENSOR','zgroup', 'MSKCCPrognosticGroup']] #cf[['PFS_MONTHS','PFS_CENSOR','AGE','SEX','RACE','MSKCC_RISK_SCORE','IMDC_RISK_SCORE']]
df_dummy = pd.get_dummies(df_r, drop_first=True)
colcovs = ['MSKCCPrognosticGroup_Int_Poor', 'zgroup_1']
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR', formula=' + '.join(colcovs))
comp = 'MSKCC_PARP1-multiv'
for cov in colcovs:
    comp = cov + '-multiv'
    hr = np.round(cph.summary.T[cov]['exp(coef)'], decimals=2)
    ci = str(np.round(cph.summary.T[cov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[cov]['exp(coef) upper 95%'], decimals=2))
    pval = cph.summary.T[cov]['p']
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)

# CoxPHFitter - PARP1 high+PBRM1wt versus others (univariate)
cph = CoxPHFitter()
df_r = cf.copy()
df_r['PARP1highVSothers'] = ''
df_r.loc[df_r['PARP1_PBRM1_term'] == 'PARP1-High/PBRM1-wt', 'PARP1highVSothers'] = 'PARP1-High/PBRM1-wt'
df_r.loc[df_r['PARP1_PBRM1_term'] != 'PARP1-High/PBRM1-wt', 'PARP1highVSothers'] = 'other'
df_r = df_r[['PFS_MONTHS','PFS_CENSOR','PARP1highVSothers']]# 'MSKCCPrognosticGroup']]

df_dummy = pd.get_dummies(df_r, drop_first=True)
colcovs = ['PARP1highVSothers_other'] #,'MSKCCPrognosticGroup_Int_Poor']
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR')
comp = 'PARP1highPBRM1wtVSothers_univ' #comp = 'MSKCC_PARP1-multiv'
for cov in colcovs:
    comp = cov + '-univ'
    hr = np.round(cph.summary.T[cov]['exp(coef)'], decimals=2)
    ci = str(np.round(cph.summary.T[cov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[cov]['exp(coef) upper 95%'], decimals=2))
    pval = cph.summary.T[cov]['p']
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)

# CoxPHFitter - PARP1 high+PBRM1wt versus others (multivariate)
cph = CoxPHFitter()
df_r = cf.copy()
df_r['PARP1highVSothers'] = ''
df_r.loc[df_r['PARP1_PBRM1_term'] == 'PARP1-High/PBRM1-wt', 'PARP1highVSothers'] = 'PARP1-High/PBRM1-wt'
df_r.loc[df_r['PARP1_PBRM1_term'] != 'PARP1-High/PBRM1-wt', 'PARP1highVSothers'] = 'other'
df_r = df_r[['PFS_MONTHS','PFS_CENSOR','PARP1highVSothers', 'MSKCCPrognosticGroup']]

df_dummy = pd.get_dummies(df_r, drop_first=True)
colcovs = ['PARP1highVSothers_other', 'MSKCCPrognosticGroup_Int_Poor'] #,'MSKCCPrognosticGroup_Int_Poor']
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR', formula=' + '.join(colcovs))
comp = 'PARP1highPBRM1wtVSothers_multiv' #comp = 'MSKCC_PARP1-multiv'
for cov in colcovs:
    comp = cov + '-multiv'
    hr = np.round(cph.summary.T[cov]['exp(coef)'], decimals=2)
    ci = str(np.round(cph.summary.T[cov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[cov]['exp(coef) upper 95%'], decimals=2))
    pval = cph.summary.T[cov]['p']
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)

#### Mean Z-score (PARP1 HIGH/ LOW +  PBRM1 Mut/ WT) - CENSORED only
T = cf['PFS_MONTHS'].loc[(cf['zgroup'] == 1) & (cf['PBRM1_mut'] == 1)]
E = cf['PFS_CENSOR'].loc[(cf['zgroup'] == 1) & (cf['PBRM1_mut'] == 1)]
T1 = cf['PFS_MONTHS'].loc[(cf['zgroup'] == 1) & (cf['PBRM1_mut'] == 0)]
E1 = cf['PFS_CENSOR'].loc[(cf['zgroup'] == 1) & (cf['PBRM1_mut'] == 0)]

T2 = cf['PFS_MONTHS'].loc[(cf['zgroup'] == 0) & (cf['PBRM1_mut'] == 1)]
E2 = cf['PFS_CENSOR'].loc[(cf['zgroup'] == 0) & (cf['PBRM1_mut'] == 1)]
T3 = cf['PFS_MONTHS'].loc[(cf['zgroup'] == 0) & (cf['PBRM1_mut'] == 0)]
E3 = cf['PFS_CENSOR'].loc[(cf['zgroup'] == 0) & (cf['PBRM1_mut'] == 0)]

kmf = KaplanMeierFitter()
kmf.fit(T, event_observed=E, label="PARP1-high / PBRM1mut")

kmf = KaplanMeierFitter()
kmf.fit(T1, event_observed=E1, label="PARP1-high / PBRM1wt")

kmf = KaplanMeierFitter()
kmf.fit(T2, event_observed=E2, label="PARP1-low / PBRM1mut")

kmf = KaplanMeierFitter()
kmf.fit(T3, event_observed=E3, label="PARP1-low / PBRM1wt")

# T2 is PARP1-low and PBRM1-mut
logrank_stats = pd.DataFrame(columns=['ARM', 'comp', 'tscore', 'p'])

pvals = dict()
results=logrank_test(T2,T,event_observed_A=E2, event_observed_B=E)
comp = 'low-mut vs high-mut'
logrank_stats = logrank_stats.append(pd.Series([arm, comp, results.test_statistic, results.p_value], index = ['ARM', 'comp', 'tscore', 'p']), ignore_index=True)

if results.p_value < 0.05:
    print('PARP1-low&PBRMmut vs PARP1-high&PBRM1mut: ' + str(results.p_value))
    pvals['PARP1-low&PBRMmut vs PARP1-high&PBRM1mut: '] = str("{0:.3f}".format(results.p_value))

results=logrank_test(T2,T1,event_observed_A=E2, event_observed_B=E1)
comp = 'low-mut vs high-wt'
logrank_stats = logrank_stats.append(pd.Series([arm, comp, results.test_statistic, results.p_value], index = ['ARM', 'comp', 'tscore', 'p']), ignore_index=True)

if results.p_value < 0.05:
    print('PARP1-low&PBRMmut vs PARP1-high&PBRM1wt: ' + str(results.p_value))
    pvals['PARP1-low&PBRMmut vs PARP1-high&PBRM1wt: '] = str("{0:.3f}".format(results.p_value))


results=logrank_test(T2,T3,event_observed_A=E2, event_observed_B=E3)
comp = 'low-mut vs low-wt'
logrank_stats = logrank_stats.append(pd.Series([arm, comp, results.test_statistic, results.p_value], index = ['ARM', 'comp', 'tscore', 'p']), ignore_index=True)
if results.p_value < 0.05:
    print('PARP1-low&PBRMmut vs PARP1-low&PBRM1wt: ' + str(results.p_value))
    pvals['PARP1-low&PBRMmut vs PARP1-low&PBRM1wt: '] = str("{0:.3f}".format(results.p_value))

results=logrank_test(T,T1,event_observed_A=E, event_observed_B=E1)
comp = 'high-mut vs high-wt'
logrank_stats = logrank_stats.append(pd.Series([arm, comp, results.test_statistic, results.p_value], index = ['ARM', 'comp', 'tscore', 'p']), ignore_index=True)
if results.p_value < 0.05:
    print('PARP1-high&PBRMmut vs PARP1-high&PBRM1wt: ' + str(results.p_value))
    pvals['PARP1-high&PBRMmut vs PARP1-high&PBRM1wt: '] = str("{0:.3f}".format(results.p_value))

results=logrank_test(T1,pd.concat([T,T2,T3]),event_observed_A=E1, event_observed_B=pd.concat([E,E2,E3]))
comp = 'high-wt vs others'
logrank_stats = logrank_stats.append(pd.Series([arm, comp, results.test_statistic, results.p_value], index = ['ARM', 'comp', 'tscore', 'p']), ignore_index=True)
if results.p_value < 0.05:
    print('PARP1-high&PBRMwt vs others: ' + str(results.p_value))
    pvals['PARP1-high&PBRMwt vs others: '] = str("{0:.3f}".format(results.p_value))

for cov in colcovs:
    comp = cov + '-multiv'
    hr = np.round(cph.summary.T[cov]['exp(coef)'], decimals=2)
    ci = str(np.round(cph.summary.T[cov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[cov]['exp(coef) upper 95%'], decimals=2))
    pval = cph.summary.T[cov]['p']
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)


cph = CoxPHFitter()
df_r = cf[['zgroup','PFS_MONTHS','PFS_CENSOR','AGE','SEX', 'MSKCC_RISK_SCORE', 'PBRM1_mut']]#,'zgroup', #cf[['PFS_MONTHS','PFS_CENSOR','AGE','SEX','RACE','MSKCC_RISK_SCORE','IMDC_RISK_SCORE']]
df_dummy = pd.get_dummies(df_r, drop_first=True)
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR', formula='zgroup_1 + AGE + SEX_M + MSKCC_RISK_SCORE_Low + MSKCC_RISK_SCORE_Intermediate + PBRM1_mut')
cph.print_summary()

# CoxPHFitter - PARP1 high+PBRM1wt versus others (univariate)
cph = CoxPHFitter()
df_r = cf.copy()
df_r['PARP1highVSothers'] = ''
df_r.loc[df_r['PARP1_PBRM1_term'] == 'PARP1-High/PBRM1-wt', 'PARP1highVSothers'] = 'PARP1-High/PBRM1-wt'
df_r.loc[df_r['PARP1_PBRM1_term'] != 'PARP1-High/PBRM1-wt', 'PARP1highVSothers'] = 'other'
df_r = df_r[['PFS_MONTHS','PFS_CENSOR','PARP1highVSothers']]# 'MSKCCPrognosticGroup']]

df_dummy = pd.get_dummies(df_r, drop_first=True)
colcovs = ['PARP1highVSothers_other'] #,'MSKCCPrognosticGroup_Int_Poor']
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR')
comp = 'PARP1highPBRM1wtVSothers_univ' #comp = 'MSKCC_PARP1-multiv'
for cov in colcovs:
    comp = cov + '-univ'
    hr = np.round(cph.summary.T[cov]['exp(coef)'], decimals=2)
    ci = str(np.round(cph.summary.T[cov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[cov]['exp(coef) upper 95%'], decimals=2))
    pval = cph.summary.T[cov]['p']
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)


# CoxPHFitter - PARP1 high+PBRM1wt versus others (multivariate)
cph = CoxPHFitter()
df_r = cf.copy()
df_r['PARP1highVSothers'] = ''
df_r.loc[df_r['PARP1_PBRM1_term'] == 'PARP1-High/PBRM1-wt', 'PARP1highVSothers'] = 'PARP1-High/PBRM1-wt'
df_r.loc[df_r['PARP1_PBRM1_term'] != 'PARP1-High/PBRM1-wt', 'PARP1highVSothers'] = 'other'
df_r = df_r[['PFS_MONTHS','PFS_CENSOR','PARP1highVSothers', 'MSKCCPrognosticGroup']]

df_dummy = pd.get_dummies(df_r, drop_first=True)
colcovs = ['PARP1highVSothers_other', 'MSKCCPrognosticGroup_Int_Poor'] #,'MSKCCPrognosticGroup_Int_Poor']
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR', formula=' + '.join(colcovs))
comp = 'PARP1highPBRM1wtVSothers_multiv' #comp = 'MSKCC_PARP1-multiv'
for cov in colcovs:
    comp = cov + '-multiv'
    hr = np.round(cph.summary.T[cov]['exp(coef)'], decimals=2)
    ci = str(np.round(cph.summary.T[cov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[cov]['exp(coef) upper 95%'], decimals=2))
    pval = cph.summary.T[cov]['p']
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)


# CoxPHFitter - PARP1 high+PBRM1wt versus PARP1high+PBRM1mut
cph = CoxPHFitter()
df_r = cf.loc[cf['PARP1_PBRM1_term'].str.contains('|'.join(['PARP1-High/PBRM1-wt', 'PARP1-High/PBRM1-mut']))]
df_r = df_r[['PFS_MONTHS','PFS_CENSOR','PARP1_PBRM1_term']]# 'MSKCCPrognosticGroup']]

df_dummy = pd.get_dummies(df_r, drop_first=True)
colcovs = ['PARP1_PBRM1_term_PARP1-High/PBRM1-wt'] #,'MSKCCPrognosticGroup_Int_Poor']
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR')
comp = 'PARP1high_univ' #comp = 'MSKCC_PARP1-multiv'
for cov in colcovs:
    comp = cov + '-univ'
    hr = np.round(cph.summary.T[cov]['exp(coef)'], decimals=2)
    ci = str(np.round(cph.summary.T[cov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[cov]['exp(coef) upper 95%'], decimals=2))
    pval = cph.summary.T[cov]['p']
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)


################################################################
######################### Sunitinib ############################
################################################################
cf = cf_aux.copy()
arm = 'Sunitinib'
cf = cf.loc[cf['ARM'] == arm]
#### Mean Z-score (PARP1 HIGH/ LOW) - Censored and Not Censored
kmf = KaplanMeierFitter()
T1 = cf['PFS_MONTHS'].loc[cf['zgroup'] == 1]
E1 = cf['PFS_CENSOR'].loc[cf['zgroup'] == 1]
kmf.fit(T1, event_observed=E1, label="PARP1-high")
kmfaux = kmf

kmf = KaplanMeierFitter()
T1 = cf['PFS_MONTHS'].loc[cf['zgroup'] == 0]
E1 = cf['PFS_CENSOR'].loc[cf['zgroup'] == 0]
kmf.fit(T1, event_observed=E1, label="PARP1-low")

results=logrank_test(T1,T,event_observed_A=E1, event_observed_B=E)
offset = 0.08
if results.p_value < 0.05:
    print('SUTEN: PARP1-low vs PARP1-high: ' + str(results.p_value))
    offset = offset + 0.05

# CoxPHFitter - PARP1 Univariate
cph = CoxPHFitter()
df_r = cf[['zgroup','PFS_MONTHS','PFS_CENSOR']]
df_dummy = pd.get_dummies(df_r, drop_first=True)
colcov = 'zgroup_1'
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR')
#cph.print_summary()
comp = 'PARP1-univ'
hr = np.round(cph.summary.T[colcov]['exp(coef)'], decimals=2)
ci = str(np.round(cph.summary.T[colcov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[colcov]['exp(coef) upper 95%'], decimals=2))
pval = cph.summary.T[colcov]['p']
stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)

pval = cph.summary.T[colcov]['p']

# CoxPHFitter - MSKCC Univariate
cph = CoxPHFitter()
df_r = cf[['PFS_MONTHS','PFS_CENSOR','MSKCCPrognosticGroup']]
df_dummy = pd.get_dummies(df_r, drop_first=True)
colcov = 'MSKCCPrognosticGroup_Int_Poor'
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR')
comp = 'MSKCC-univ'
hr = np.round(cph.summary.T[colcov]['exp(coef)'], decimals=2)
ci = str(np.round(cph.summary.T[colcov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[colcov]['exp(coef) upper 95%'], decimals=2))
pval = cph.summary.T[colcov]['p']
stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)

# CoxPHFitter - PARP1, MSKCC Multivariate
cph = CoxPHFitter()
df_r = cf[['PFS_MONTHS','PFS_CENSOR','zgroup', 'MSKCCPrognosticGroup']] #cf[['PFS_MONTHS','PFS_CENSOR','AGE','SEX','RACE','MSKCC_RISK_SCORE','IMDC_RISK_SCORE']]
df_dummy = pd.get_dummies(df_r, drop_first=True)
colcovs = ['MSKCCPrognosticGroup_Int_Poor', 'zgroup_1']
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR', formula=' + '.join(colcovs))
comp = 'MSKCC_PARP1-multiv'
for cov in colcovs:
    comp = cov + '-multiv'
    hr = np.round(cph.summary.T[cov]['exp(coef)'], decimals=2)
    ci = str(np.round(cph.summary.T[cov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[cov]['exp(coef) upper 95%'], decimals=2))
    pval = cph.summary.T[cov]['p']
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)

# CoxPHFitter - PARP1 high+PBRM1wt versus others (univariate)
cph = CoxPHFitter()
df_r = cf.copy()
df_r['PARP1highVSothers'] = ''
df_r.loc[df_r['PARP1_PBRM1_term'] == 'PARP1-High/PBRM1-wt', 'PARP1highVSothers'] = 'PARP1-High/PBRM1-wt'
df_r.loc[df_r['PARP1_PBRM1_term'] != 'PARP1-High/PBRM1-wt', 'PARP1highVSothers'] = 'other'
df_r = df_r[['PFS_MONTHS','PFS_CENSOR','PARP1highVSothers']]# 'MSKCCPrognosticGroup']]

df_dummy = pd.get_dummies(df_r, drop_first=True)
colcovs = ['PARP1highVSothers_other'] #,'MSKCCPrognosticGroup_Int_Poor']
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR')
comp = 'PARP1highPBRM1wtVSothers_univ' #comp = 'MSKCC_PARP1-multiv'
for cov in colcovs:
    comp = cov + '-univ'
    hr = np.round(cph.summary.T[cov]['exp(coef)'], decimals=2)
    ci = str(np.round(cph.summary.T[cov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[cov]['exp(coef) upper 95%'], decimals=2))
    pval = cph.summary.T[cov]['p']
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)

# CoxPHFitter - PARP1 high+PBRM1wt versus others (multivariate)
cph = CoxPHFitter()
df_r = cf.copy()
df_r['PARP1highVSothers'] = ''
df_r.loc[df_r['PARP1_PBRM1_term'] == 'PARP1-High/PBRM1-wt', 'PARP1highVSothers'] = 'PARP1-High/PBRM1-wt'
df_r.loc[df_r['PARP1_PBRM1_term'] != 'PARP1-High/PBRM1-wt', 'PARP1highVSothers'] = 'other'
df_r = df_r[['PFS_MONTHS','PFS_CENSOR','PARP1highVSothers', 'MSKCCPrognosticGroup']]

df_dummy = pd.get_dummies(df_r, drop_first=True)
colcovs = ['PARP1highVSothers_other', 'MSKCCPrognosticGroup_Int_Poor'] #,'MSKCCPrognosticGroup_Int_Poor']
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR', formula=' + '.join(colcovs))
comp = 'PARP1highPBRM1wtVSothers_multiv' #comp = 'MSKCC_PARP1-multiv'
for cov in colcovs:
    comp = cov + '-multiv'
    hr = np.round(cph.summary.T[cov]['exp(coef)'], decimals=2)
    ci = str(np.round(cph.summary.T[cov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[cov]['exp(coef) upper 95%'], decimals=2))
    pval = cph.summary.T[cov]['p']
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)

#### Mean Z-score (PARP1 HIGH/ LOW +  PBRM1 Mut/ WT) - CENSORED only
T = cf['PFS_MONTHS'].loc[(cf['zgroup'] == 1) & (cf['PBRM1_mut'] == 1)]
E = cf['PFS_CENSOR'].loc[(cf['zgroup'] == 1) & (cf['PBRM1_mut'] == 1)]
T1 = cf['PFS_MONTHS'].loc[(cf['zgroup'] == 1) & (cf['PBRM1_mut'] == 0)]
E1 = cf['PFS_CENSOR'].loc[(cf['zgroup'] == 1) & (cf['PBRM1_mut'] == 0)]

T2 = cf['PFS_MONTHS'].loc[(cf['zgroup'] == 0) & (cf['PBRM1_mut'] == 1)]
E2 = cf['PFS_CENSOR'].loc[(cf['zgroup'] == 0) & (cf['PBRM1_mut'] == 1)]
T3 = cf['PFS_MONTHS'].loc[(cf['zgroup'] == 0) & (cf['PBRM1_mut'] == 0)]
E3 = cf['PFS_CENSOR'].loc[(cf['zgroup'] == 0) & (cf['PBRM1_mut'] == 0)]

kmf = KaplanMeierFitter()
kmf.fit(T, event_observed=E, label="PARP1-high / PBRM1mut")

kmf = KaplanMeierFitter()
kmf.fit(T1, event_observed=E1, label="PARP1-high / PBRM1wt")

kmf = KaplanMeierFitter()
kmf.fit(T2, event_observed=E2, label="PARP1-low / PBRM1mut")

kmf = KaplanMeierFitter()
kmf.fit(T3, event_observed=E3, label="PARP1-low / PBRM1wt")

# T2 is PARP1-low and PBRM1-mut
pvals = dict()
results=logrank_test(T2,T,event_observed_A=E2, event_observed_B=E)
comp = 'low-mut vs high-mut'
logrank_stats = logrank_stats.append(pd.Series([arm, comp, results.test_statistic, results.p_value], index = ['ARM', 'comp', 'tscore', 'p']), ignore_index=True)

if results.p_value < 0.05:
    print('PARP1-low&PBRMmut vs PARP1-high&PBRM1mut: ' + str(results.p_value))
    pvals['PARP1-low&PBRMmut vs PARP1-high&PBRM1mut: '] = str("{0:.3f}".format(results.p_value))

results=logrank_test(T2,T1,event_observed_A=E2, event_observed_B=E1)
comp = 'low-mut vs high-wt'
logrank_stats = logrank_stats.append(pd.Series([arm, comp, results.test_statistic, results.p_value], index = ['ARM', 'comp', 'tscore', 'p']), ignore_index=True)

if results.p_value < 0.05:
    print('PARP1-low&PBRMmut vs PARP1-high&PBRM1wt: ' + str(results.p_value))
    pvals['PARP1-low&PBRMmut vs PARP1-high&PBRM1wt: '] = str("{0:.3f}".format(results.p_value))

results=logrank_test(T2,T3,event_observed_A=E2, event_observed_B=E3)
comp = 'low-mut vs low-wt'
logrank_stats = logrank_stats.append(pd.Series([arm, comp, results.test_statistic, results.p_value], index = ['ARM', 'comp', 'tscore', 'p']), ignore_index=True)
if results.p_value < 0.05:
    print('PARP1-low&PBRMmut vs PARP1-low&PBRM1wt: ' + str(results.p_value))
    pvals['PARP1-low&PBRMmut vs PARP1-low&PBRM1wt: '] = str("{0:.3f}".format(results.p_value))

results=logrank_test(T,T1,event_observed_A=E, event_observed_B=E1)
comp = 'high-mut vs high-wt'
logrank_stats = logrank_stats.append(pd.Series([arm, comp, results.test_statistic, results.p_value], index = ['ARM', 'comp', 'tscore', 'p']), ignore_index=True)
if results.p_value < 0.05:
    print('PARP1-high&PBRMmut vs PARP1-high&PBRM1wt: ' + str(results.p_value))
    pvals['PARP1-high&PBRMmut vs PARP1-high&PBRM1wt: '] = str("{0:.3f}".format(results.p_value))

results=logrank_test(T1,pd.concat([T,T2,T3]),event_observed_A=E1, event_observed_B=pd.concat([E,E2,E3]))
comp = 'high-wt vs others'
logrank_stats = logrank_stats.append(pd.Series([arm, comp, results.test_statistic, results.p_value], index = ['ARM', 'comp', 'tscore', 'p']), ignore_index=True)
if results.p_value < 0.05:
    print('PARP1-high&PBRMwt vs others: ' + str(results.p_value))
    pvals['PARP1-high&PBRMwt vs others: '] = str("{0:.3f}".format(results.p_value))

cph = CoxPHFitter()
cf['PBRM1_mut'] = np.where(cf['FMI_PBRM1'].isin(['loss','rearrangement', 'short-variant']), 1, 0)
cf['PBRM1_mut_PARP1_low'] = np.where((cf['PBRM1_mut'] == 1) & (cf['zgroup'] == 0), 1, 0)
df_r = cf[['PBRM1_mut_PARP1_low','PFS_MONTHS','PFS_CENSOR','AGE','SEX']]#'PBRM1_mut','zgroup', #cf[['PFS_MONTHS','PFS_CENSOR','AGE','SEX','RACE','MSKCC_RISK_SCORE','IMDC_RISK_SCORE']]
df_dummy = pd.get_dummies(df_r, drop_first=True)
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR')#, step_size = 0.5)   ## Fit the data to train the model
cph.print_summary()

cph = CoxPHFitter()
df_r = cf[['zgroup','PFS_MONTHS','PFS_CENSOR','AGE','SEX', 'MSKCC_RISK_SCORE', 'PBRM1_mut']]#,'zgroup', #cf[['PFS_MONTHS','PFS_CENSOR','AGE','SEX','RACE','MSKCC_RISK_SCORE','IMDC_RISK_SCORE']]
df_dummy = pd.get_dummies(df_r, drop_first=True)
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR', formula='zgroup_1 + AGE + SEX_M + MSKCC_RISK_SCORE_Low + MSKCC_RISK_SCORE_Intermediate + PBRM1_mut')
cph.print_summary()


# CoxPHFitter - PARP1 high+PBRM1wt versus others (univariate)
cph = CoxPHFitter()
df_r = cf.copy()
df_r['PARP1highVSothers'] = ''
df_r.loc[df_r['PARP1_PBRM1_term'] == 'PARP1-High/PBRM1-wt', 'PARP1highVSothers'] = 'PARP1-High/PBRM1-wt'
df_r.loc[df_r['PARP1_PBRM1_term'] != 'PARP1-High/PBRM1-wt', 'PARP1highVSothers'] = 'other'
df_r = df_r[['PFS_MONTHS','PFS_CENSOR','PARP1highVSothers']]# 'MSKCCPrognosticGroup']]

df_dummy = pd.get_dummies(df_r, drop_first=True)
colcovs = ['PARP1highVSothers_other'] #,'MSKCCPrognosticGroup_Int_Poor']
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR')
comp = 'PARP1highPBRM1wtVSothers_univ' #comp = 'MSKCC_PARP1-multiv'
for cov in colcovs:
    comp = cov + '-univ'
    hr = np.round(cph.summary.T[cov]['exp(coef)'], decimals=2)
    ci = str(np.round(cph.summary.T[cov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[cov]['exp(coef) upper 95%'], decimals=2))
    pval = cph.summary.T[cov]['p']
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)

# CoxPHFitter - PARP1 high+PBRM1wt versus others (multivariate)
cph = CoxPHFitter()
df_r = cf.copy()
df_r['PARP1highVSothers'] = ''
df_r.loc[df_r['PARP1_PBRM1_term'] == 'PARP1-High/PBRM1-wt', 'PARP1highVSothers'] = 'PARP1-High/PBRM1-wt'
df_r.loc[df_r['PARP1_PBRM1_term'] != 'PARP1-High/PBRM1-wt', 'PARP1highVSothers'] = 'other'
df_r = df_r[['PFS_MONTHS','PFS_CENSOR','PARP1highVSothers', 'MSKCCPrognosticGroup']]

df_dummy = pd.get_dummies(df_r, drop_first=True)
colcovs = ['PARP1highVSothers_other', 'MSKCCPrognosticGroup_Int_Poor'] #,'MSKCCPrognosticGroup_Int_Poor']
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR', formula=' + '.join(colcovs))
comp = 'PARP1highPBRM1wtVSothers_multiv' #comp = 'MSKCC_PARP1-multiv'
for cov in colcovs:
    comp = cov + '-multiv'
    hr = np.round(cph.summary.T[cov]['exp(coef)'], decimals=2)
    ci = str(np.round(cph.summary.T[cov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[cov]['exp(coef) upper 95%'], decimals=2))
    pval =cph.summary.T[cov]['p']
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)

# CoxPHFitter - PARP1 high+PBRM1wt versus PARP1high+PBRM1mut
cph = CoxPHFitter()
df_r = cf.loc[cf['PARP1_PBRM1_term'].str.contains('|'.join(['PARP1-High/PBRM1-wt', 'PARP1-High/PBRM1-mut']))]
df_r = df_r[['PFS_MONTHS','PFS_CENSOR','PARP1_PBRM1_term']]# 'MSKCCPrognosticGroup']]

df_dummy = pd.get_dummies(df_r, drop_first=True)
colcovs = ['PARP1_PBRM1_term_PARP1-High/PBRM1-wt'] #,'MSKCCPrognosticGroup_Int_Poor']
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR')
comp = 'PARP1high_univ' #comp = 'MSKCC_PARP1-multiv'
for cov in colcovs:
    comp = cov + '-univ'
    hr = np.round(cph.summary.T[cov]['exp(coef)'], decimals=2)
    ci = str(np.round(cph.summary.T[cov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[cov]['exp(coef) upper 95%'], decimals=2))
    pval = cph.summary.T[cov]['p']
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)



# CoxPHFitter - PARP1 high+PBRM1wt versus PARP1low+PBRM1wt
cph = CoxPHFitter()
df_r = cf.loc[cf['PARP1_PBRM1_term'].str.contains('|'.join(['PARP1-High/PBRM1-wt', 'PARP1-Low/PBRM1-wt']))]
df_r = df_r[['PFS_MONTHS','PFS_CENSOR','PARP1_PBRM1_term']]# 'MSKCCPrognosticGroup']]

df_dummy = pd.get_dummies(df_r, drop_first=True)
colcovs = ['PARP1_PBRM1_term_PARP1-Low/PBRM1-wt'] #,'MSKCCPrognosticGroup_Int_Poor']
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR')
comp = 'PARP1high_univ' #comp = 'MSKCC_PARP1-multiv'
for cov in colcovs:
    comp = cov + '-univ'
    hr = np.round(cph.summary.T[cov]['exp(coef)'], decimals=2)
    ci = str(np.round(cph.summary.T[cov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[cov]['exp(coef) upper 95%'], decimals=2))
    pval = cph.summary.T[cov]['p']
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)

# CoxPHFitter - PARP1 high+PBRM1wt versus PARP1low+PBRM1mut
cph = CoxPHFitter()
df_r = cf.loc[cf['PARP1_PBRM1_term'].str.contains('|'.join(['PARP1-High/PBRM1-wt', 'PARP1-Low/PBRM1-mut']))]
df_r = df_r[['PFS_MONTHS','PFS_CENSOR','PARP1_PBRM1_term']]# 'MSKCCPrognosticGroup']]

df_dummy = pd.get_dummies(df_r, drop_first=True)
colcovs = ['PARP1_PBRM1_term_PARP1-Low/PBRM1-mut'] #,'MSKCCPrognosticGroup_Int_Poor']
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR')
comp = 'PARP1high_univ' #comp = 'MSKCC_PARP1-multiv'
for cov in colcovs:
    comp = cov + '-univ'
    hr = np.round(cph.summary.T[cov]['exp(coef)'], decimals=2)
    ci = str(np.round(cph.summary.T[cov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[cov]['exp(coef) upper 95%'], decimals=2))
    pval = cph.summary.T[cov]['p']
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)


cf = cf_aux.copy()

stats.to_csv('Cox_analysis.csv',sep='\t',index = None)
logrank_stats.to_csv('logrank_stats_analysis.csv',sep='\t',index = None)


#############################################
####### Motzer's molecular subtypes #########
#############################################
## Analyzis of the Molecular subtypes described by Motzer et al in IMmotion151 and PARP1/PBRM1 status
# dfangio will contain the number of tumors within each molecular subtype (labeled from 1 to 7) for PARP1/PBRM1 groups
# This is useful to calculate Chi2 between PARP1-High/PBRM1-wt versus other groups
angio = cf.groupby(['ARM','PARP1_PBRM1_term','NMF_GROUP']).count()['RNASEQ_SAMPLE_ID'].reset_index()
dfangio = pd.DataFrame({'group':list(set(angio['NMF_GROUP']))})
for arm in list(set(angio['ARM'])):
    for group in list(set(angio['PARP1_PBRM1_term'])):
        dff = angio.loc[(angio['ARM'] == arm) & (angio['PARP1_PBRM1_term'] == group)]
        dfangio[arm + '_' + group] = dfangio['group'].map(dict(zip(list(dff['NMF_GROUP']), list(dff['RNASEQ_SAMPLE_ID']))))

for i in [col for col in list(dfangio.columns) if 'group' not in col]:
    dfangio['pct_' + i] = 100*(dfangio[i] / dfangio[i].sum())

molecular_terms = ['Angio/stroma', 'Angiogenic', 'Complement-Ox', 'T-eff/ Proliferative', 'Proliferative', 'Stromal/Proliferative', 'snoRNA']
group = list(dfangio['group'])
d = dict(zip(group,molecular_terms))
dfangio['molecular_subgroups'] = dfangio['group'].map(d) # Add molecular subgroup terms to final table

dfangio.to_csv('molecular_subtype_tto_parp1ANDpbrm1.csv',sep='\t')

# Calculate Chi2 test:
####### chi2 to compare mol.sign groups PARP1highPBRM1wt vs the rest of groups
from scipy import stats
a = pd.read_csv('molecular_subtype_tto_parp1ANDpbrm1.csv',sep='\t')
a = a[[x for x in list(a.columns) if not 'pct_' in x]] # eliminate columns with percentages
mergeCols = list(set([x.split('_')[1] for x in list(a.columns) if not 'group' in x]))
for cols in mergeCols:
    a[cols] = a['Atezolizumab+Bevacizumab_' + cols] + a['Sunitinib_' + cols]
a = a[['group', 'molecular_subgroups'] + mergeCols]
a.columns= a.columns.str.replace('/','_')

maingroup = 'PARP1-High_PBRM1-wt'
groups2compare = ['PARP1-High_PBRM1-mut','PARP1-Low_PBRM1-mut','PARP1-Low_PBRM1-wt']

l = list()
for x in list(a.sum(axis=0).values):
    if not isinstance(x, str):
        l.append(x)
    else:
        l.append('all')
l[0] = 'NAN'
cols = ['group', 'molecular_subgroups'] + [maingroup] + groups2compare
a = a.append(pd.Series(l, index = cols), ignore_index=True)

molecular_subgroups = [x for x in list(a['molecular_subgroups']) if x != 'all']

for compgroup in groups2compare:
    pvals = list()
    for molgroup in molecular_subgroups:
        contTable = pd.DataFrame([[a[maingroup].loc[a['molecular_subgroups'] == molgroup].values[0],a[compgroup].loc[a['molecular_subgroups'] == molgroup].values[0]],[a[maingroup].loc[a['molecular_subgroups'] == 'all'].values[0]- a[maingroup].loc[a['molecular_subgroups'] == molgroup].values[0],a[compgroup].loc[a['molecular_subgroups'] == 'all'].values[0]- a[compgroup].loc[a['molecular_subgroups'] == molgroup].values[0]]])
        chi2, p, dof, ex = stats.chi2_contingency(contTable)
        pvals.append(p)
    pvals.append('NAN')
    colname = 'p_' + compgroup
    a[colname] = pvals
a.to_csv('angiogroups_pvals.csv',sep='\t',index = None)
