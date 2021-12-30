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

c['PARP1'] = c['RNASEQ_SAMPLE_ID'].map(d)
# Calculate mean Z-scores for PARP1
c['zPARP1'] = (c['PARP1'] - c['PARP1'].mean()) / (c['PARP1'].std())
c['zmedianPARP1'] = (c['PARP1'] - c['PARP1'].median()) / (c['PARP1'].std())

# cf will be a dataframe that contains both clinical and expression data
cf = c.loc[~c['zPARP1'].isnull()]
cf['zgroup'] = ''
cf['zmediangroup'] = ''
cf.loc[cf['zPARP1'] < 0, 'zgroup'] = 0 # low PARP1, group 0
cf.loc[cf['zPARP1'] > 0, 'zgroup'] = 1
cf.loc[cf['zmedianPARP1'] < 0, 'zmediangroup'] = 0
cf.loc[cf['zmedianPARP1'] > 0, 'zmediangroup'] = 1

# Keep only samples sequenced (there are samples that do not contain their mutational status)
cf = cf.loc[~cf['FMI_SAMPLE_ID'].isnull()]

# Binarize PBRM1 mutational status and prepare PARP1+PBRM1 info for each sample
cf['PBRM1_mut'] = np.where(cf['FMI_PBRM1'].isin(['loss','rearrangement', 'short-variant']), 1, 0)
cf['PBRM1_mut_term'] = np.where(cf['FMI_PBRM1'].isin(['loss','rearrangement', 'short-variant']), 'PBRM1-mut', 'PBRM1-wt')
cf['PARP1_zgroup_term'] = np.where(cf['zgroup']  == 1, 'PARP1-High', 'PARP1-Low')
cf['PARP1_PBRM1'] = cf['zgroup'].astype(str) + '_' + cf['PBRM1_mut'].astype(str)
cf['PARP1_PBRM1_term'] = cf['PARP1_zgroup_term'].astype(str) + '/' + cf['PBRM1_mut_term'].astype(str)
# MSKCC prognostic group
cf['MSKCCPrognosticGroup'] = np.where(cf['MSKCC_RISK_SCORE'] == 'Low', 'Favorable', 'Int_Poor')
# Keep samples with ccRCC histology:
cf = cf.loc[cf['HISTOLOGY_SARCOMATOID'].isin(['ccRCC_Sarc', 'ccRCC_nonSarc'])].copy()

#############################################################################

# Save the initial dataframe into a safe copy (cf_aux)
cf_aux = cf.copy()
# Create the figure
fig, ax = plt.subplots(2, 2, figsize=(10, 10))
fslabel = 14
fs = 10
# stats dataframe will serve to store output results from statistical analysis
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
kmf.plot_survival_function(ax=ax[0,0], show_censors = True, ci_show=False, color='k', marker='x', markersize = 1)
kmfaux = kmf

kmf = KaplanMeierFitter()
T1 = cf['PFS_MONTHS'].loc[cf['zgroup'] == 0]
E1 = cf['PFS_CENSOR'].loc[cf['zgroup'] == 0]
kmf.fit(T1, event_observed=E1, label="PARP1-low")
kmf.plot_survival_function(ax=ax[0,0], show_censors = True, ci_show=False, color='gray', marker='x', markersize = 1)
add_at_risk_counts(kmfaux, kmf, ax=ax[0,0])

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
ax[0,0].text(0, 0+ offset, 'HR=' + str(hr) + ' (95% CI=' + ci  + '), ' + 'p=' +  str(pval), fontsize=fs_text, bbox=dict(facecolor='none', edgecolor='black', boxstyle='square,pad=0.7'))

# Set figure label: A
ax[0,0].text(-7, 1.1, 'A', fontsize=fslabel)

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
    pval = np.round(cph.summary.T[cov]['p'], decimals=3)
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
    pval = np.round(cph.summary.T[cov]['p'], decimals=3)
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
    pval = np.round(cph.summary.T[cov]['p'], decimals=3)
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)


#cph.plot()
#plt.show()

ax[0,0].set_xlabel('Time (months)')
ax[0,0].set_ylabel("Progression Free Survival Probability")
ax[0,0].title.set_text('Atezolizumab plus Bevacizumab arm') # 'ATEZO/BEVA: PARP1-HIGH/Low and PBRM1-mut/wt'


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
kmf.plot_survival_function(ax=ax[1,0], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)
#add_at_risk_counts(kmf, ax=ax[0,1])


kmf = KaplanMeierFitter()
kmf.fit(T1, event_observed=E1, label="PARP1-high / PBRM1wt")
kmf.plot_survival_function(ax=ax[1,0], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)
#add_at_risk_counts(kmf, ax=ax[0,1])


kmf = KaplanMeierFitter()
kmf.fit(T2, event_observed=E2, label="PARP1-low / PBRM1mut")
kmf.plot_survival_function(ax=ax[1,0], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)
#add_at_risk_counts(kmf, ax=ax[0,1])

kmf = KaplanMeierFitter()
kmf.fit(T3, event_observed=E3, label="PARP1-low / PBRM1wt")
kmf.plot_survival_function(ax=ax[1,0], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)
#add_at_risk_counts(kmf, ax=ax[0,1])

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
    pval = np.round(cph.summary.T[cov]['p'], decimals=3)
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)



offset = 0.08
for key in pvals.keys():
    #ax[1,0].text(0, 0+ offset, key + pvals[key], fontsize=fs)
    offset = offset + 0.05


cph = CoxPHFitter()
df_r = cf[['zgroup','PFS_MONTHS','PFS_CENSOR','AGE','SEX', 'MSKCC_RISK_SCORE', 'PBRM1_mut']]#,'zgroup', #cf[['PFS_MONTHS','PFS_CENSOR','AGE','SEX','RACE','MSKCC_RISK_SCORE','IMDC_RISK_SCORE']]
df_dummy = pd.get_dummies(df_r, drop_first=True)
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR', formula='zgroup_1 + AGE + SEX_M + MSKCC_RISK_SCORE_Low + MSKCC_RISK_SCORE_Intermediate + PBRM1_mut')
cph.print_summary()
#cph.plot()
#plt.show()

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
    pval = np.round(cph.summary.T[cov]['p'], decimals=3)
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
    pval = np.round(cph.summary.T[cov]['p'], decimals=3)
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
    pval = np.round(cph.summary.T[cov]['p'], decimals=3)
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)


# Set figure label: C
ax[1,0].text(-7, 1.1, 'C', fontsize=fslabel)


ax[1,0].set_xlabel('Time (months)')
ax[1,0].set_ylabel("Progression Free Survival Probability")
ax[1,0].title.set_text('Atezolizumab plus Bevacizumab arm') # 'ATEZO/BEVA: PARP1-HIGH/Low and PBRM1-mut/wt'





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
kmf.plot_survival_function(ax=ax[0,1], show_censors = True, ci_show=False, color='k', marker='x', markersize = 1)
kmfaux = kmf

kmf = KaplanMeierFitter()
T1 = cf['PFS_MONTHS'].loc[cf['zgroup'] == 0]
E1 = cf['PFS_CENSOR'].loc[cf['zgroup'] == 0]
kmf.fit(T1, event_observed=E1, label="PARP1-low")
kmf.plot_survival_function(ax=ax[0,1], show_censors = True, ci_show=False, color='gray', marker='x', markersize = 1)
add_at_risk_counts(kmfaux, kmf, ax=ax[0,1])

# Set figure label: B
ax[0,1].text(-5, 1.1, 'B', fontsize=fslabel)


ax[0,1].set_xlabel('Time (months)')
#ax[0,1].set_ylabel("Progression Free Survival Probability")
ax[0,1].title.set_text('Sunitinib arm')


results=logrank_test(T1,T,event_observed_A=E1, event_observed_B=E)
offset = 0.08
if results.p_value < 0.05:
    print('SUTEN: PARP1-low vs PARP1-high: ' + str(results.p_value))
    #ax[0,1].text(0, 0+ offset, 'PARP1-low vs PARP1-high: ' + str("{0:.3f}".format(results.p_value)), fontsize=fs)
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
ax[0,1].text(0, 0+ offset, 'HR=' + str(hr) + ' (95% CI=' + ci  + '), ' + 'p=' +  str(pval), fontsize=fs_text, bbox=dict(facecolor='none', edgecolor='black', boxstyle='square,pad=0.7'))

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
    pval = np.round(cph.summary.T[cov]['p'], decimals=3)
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
    pval = np.round(cph.summary.T[cov]['p'], decimals=3)
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
    pval = np.round(cph.summary.T[cov]['p'], decimals=3)
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
kmf.plot_survival_function(ax=ax[1,1], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)
#add_at_risk_counts(kmf, ax=ax[1,1])


kmf = KaplanMeierFitter()
kmf.fit(T1, event_observed=E1, label="PARP1-high / PBRM1wt")
kmf.plot_survival_function(ax=ax[1,1], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)
#add_at_risk_counts(kmf, ax=ax[1,1])


kmf = KaplanMeierFitter()
kmf.fit(T2, event_observed=E2, label="PARP1-low / PBRM1mut")
kmf.plot_survival_function(ax=ax[1,1], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)
#add_at_risk_counts(kmf, ax=ax[1,1])

kmf = KaplanMeierFitter()
kmf.fit(T3, event_observed=E3, label="PARP1-low / PBRM1wt")
kmf.plot_survival_function(ax=ax[1,1], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)
#add_at_risk_counts(kmf, ax=ax[1,1])

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



offset = 0.08
for key in pvals.keys():
    #ax[1,1].text(0, 0+ offset, key + pvals[key], fontsize=fs)
    offset = offset + 0.05

# Set figure label: D
ax[1,1].text(-5, 1.1, 'D', fontsize=fslabel)


ax[1,1].set_xlabel('Time (months)')
#ax[1,1].set_ylabel("Progression Free Survival Probability")
ax[1,1].title.set_text('Sunitinib arm')


cph = CoxPHFitter()
cf['PBRM1_mut'] = np.where(cf['FMI_PBRM1'].isin(['loss','rearrangement', 'short-variant']), 1, 0)
cf['PBRM1_mut_PARP1_low'] = np.where((cf['PBRM1_mut'] == 1) & (cf['zgroup'] == 0), 1, 0)
df_r = cf[['PBRM1_mut_PARP1_low','PFS_MONTHS','PFS_CENSOR','AGE','SEX']]#'PBRM1_mut','zgroup', #cf[['PFS_MONTHS','PFS_CENSOR','AGE','SEX','RACE','MSKCC_RISK_SCORE','IMDC_RISK_SCORE']]
df_dummy = pd.get_dummies(df_r, drop_first=True)
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR')#, step_size = 0.5)   ## Fit the data to train the model
cph.print_summary()
#cph.plot()
#plt.show()
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
    pval = np.round(cph.summary.T[cov]['p'], decimals=3)
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
    pval = np.round(cph.summary.T[cov]['p'], decimals=3)
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
    pval = np.round(cph.summary.T[cov]['p'], decimals=3)
    stats = stats.append(pd.Series([arm, comp, hr, ci, pval], index = ['ARM', 'comp', 'HR', '95CI', 'p']), ignore_index=True)




# Set the ticks and ticklabels for all axes, plot and save the figure
ytickvals = [round(x,1) for x in list(np.arange(0,1.1,0.1))]
plt.setp(ax, yticks=ytickvals)
cf = cf_aux.copy()
plt.tight_layout()
plt.savefig('/mnt/64716603-5b56-4f9a-b195-c11560647a3a/data/IMMOTION151/Figure_1.png', dpi = 1000)
plt.show()

stats.to_csv('Cox_analysis.csv',sep='\t',index = None)
logrank_stats.to_csv('logrank_stats_analysis.csv',sep='\t',index = None)


#############################################
####### Motzer's molecular subtypes #########
#############################################
# Analyzis of the Molecular subtypes described by Motzer et al in IMmotion151 and PARP1/PBRM1 status
angio = cf.groupby(['ARM','PARP1_PBRM1_term','NMF_GROUP']).count()['RNASEQ_SAMPLE_ID'].reset_index()
dfangio = pd.DataFrame({'group':list(set(angio['NMF_GROUP']))})
for arm in list(set(angio['ARM'])):
    for group in list(set(angio['PARP1_PBRM1_term'])):
        dff = angio.loc[(angio['ARM'] == arm) & (angio['PARP1_PBRM1_term'] == group)]
        dfangio[arm + '_' + group] = dfangio['group'].map(dict(zip(list(dff['NMF_GROUP']), list(dff['RNASEQ_SAMPLE_ID']))))

for i in [col for col in list(dfangio.columns) if 'group' not in col]:
    dfangio['pct_' + i] = 100*(dfangio[i] / dfangio[i].sum())

#plotangio = dfangio.set_index('group').T.reset_index()
#ax = plotangio.plot(x="index", y=1, kind="bar")
#plotangio.plot(x="index", y=2, kind="bar", ax=ax, color="C2")
