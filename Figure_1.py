# Obtain RNAseq data from the IMmotion151:
#conda activate py36; pyega3 fetch EGAF00004859518
import pandas as pd
from lifelines import KaplanMeierFitter #pip install lifelines
import matplotlib.pyplot as plt
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test
import numpy as np
from lifelines import CoxPHFitter

c = pd.read_excel('Table_S1.xlsx') # Table S1 from Motzer et al, Cancer Cell 38, 803–817, December 14, 2020
df = pd.read_csv('IMmotion151.expression.data.TPM.anon.csv',sep=',') # Transcript counts matrix from EGAF00004859518
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


#############################################################################
####################### Create the figure ###################################
#############################################################################

fig, ax = plt.subplots(2, 2, figsize=(10, 10))
plt.subplots_adjust(hspace=0.9,wspace=0.40,left=0.105, right=0.985, top=0.955, bottom=0.050)

fslabel = 14
fs = 10

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
kmf.plot_survival_function(ax=ax[0,0], show_censors = True, ci_show=False, color='gold', marker='x', markersize = 1)
add_at_risk_counts(kmfaux, kmf, ax=ax[0,0])


# CoxPHFitter - PARP1 Univariate
cph = CoxPHFitter()
df_r = cf[['zgroup','PFS_MONTHS','PFS_CENSOR']]
df_dummy = pd.get_dummies(df_r, drop_first=True)
colcov = 'zgroup_1'
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR')
comp = 'PARP1-univ'
hr = np.round(cph.summary.T[colcov]['exp(coef)'], decimals=2)
ci = str(np.round(cph.summary.T[colcov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[colcov]['exp(coef) upper 95%'], decimals=2))
pval = cph.summary.T[colcov]['p']

offset = 0.08
pval = np.round(cph.summary.T[colcov]['p'], decimals = 4)
fs_text = 10
ax[0,0].text(0, 0+ offset, 'HR=' + str(hr) + ' (95% CI=' + ci  + '), ' + 'p=' +  str(pval), fontsize=fs_text, bbox=dict(facecolor='none', edgecolor='black', boxstyle='square,pad=0.7'))


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

kmf = KaplanMeierFitter()
kmf.fit(T1, event_observed=E1, label="PARP1-high / PBRM1wt")
kmf.plot_survival_function(ax=ax[1,0], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)

kmf = KaplanMeierFitter()
kmf.fit(T2, event_observed=E2, label="PARP1-low / PBRM1mut")
kmf.plot_survival_function(ax=ax[1,0], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)

kmf = KaplanMeierFitter()
kmf.fit(T3, event_observed=E3, label="PARP1-low / PBRM1wt")
kmf.plot_survival_function(ax=ax[1,0], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)


################################################################
######################### Sunitinib ############################
################################################################
cf = cf_aux.copy()
arm = 'Sunitinib'
cf = cf.loc[cf['ARM'] == arm]

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

kmf = KaplanMeierFitter()
kmf.fit(T1, event_observed=E1, label="PARP1-high / PBRM1wt")
kmf.plot_survival_function(ax=ax[1,1], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)

kmf = KaplanMeierFitter()
kmf.fit(T2, event_observed=E2, label="PARP1-low / PBRM1mut")
kmf.plot_survival_function(ax=ax[1,1], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)

kmf = KaplanMeierFitter()
kmf.fit(T3, event_observed=E3, label="PARP1-low / PBRM1wt")
kmf.plot_survival_function(ax=ax[1,1], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)


#### Mean Z-score (PARP1 HIGH/ LOW) - Censored and Not Censored
kmf = KaplanMeierFitter()
T = cf['PFS_MONTHS'].loc[cf['zgroup'] == 1]
E = cf['PFS_CENSOR'].loc[cf['zgroup'] == 1]
kmf.fit(T, event_observed=E, label="PARP1-high")
kmf.plot_survival_function(ax=ax[0,1], show_censors = True, ci_show=False, color='k', marker='x', markersize = 1)
kmfaux = kmf

kmf = KaplanMeierFitter()
T1 = cf['PFS_MONTHS'].loc[cf['zgroup'] == 0]
E1 = cf['PFS_CENSOR'].loc[cf['zgroup'] == 0]
kmf.fit(T1, event_observed=E1, label="PARP1-low")
kmf.plot_survival_function(ax=ax[0,1], show_censors = True, ci_show=False, color='gold', marker='x', markersize = 1)
add_at_risk_counts(kmfaux, kmf, ax=ax[0,1])


# CoxPHFitter - PARP1 Univariate
cph = CoxPHFitter()
df_r = cf[['zgroup','PFS_MONTHS','PFS_CENSOR']]
df_dummy = pd.get_dummies(df_r, drop_first=True)
colcov = 'zgroup_1'
cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR')
comp = 'PARP1-univ'
hr = np.round(cph.summary.T[colcov]['exp(coef)'], decimals=2)
ci = str(np.round(cph.summary.T[colcov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[colcov]['exp(coef) upper 95%'], decimals=2))
pval = cph.summary.T[colcov]['p']

offset = 0.08
pval = np.round(cph.summary.T[colcov]['p'], decimals = 4)
fs_text = 10
ax[0,1].text(0, 0+ offset, 'HR=' + str(hr) + ' (95% CI=' + ci  + '), ' + 'p=' +  str(pval), fontsize=fs_text, bbox=dict(facecolor='none', edgecolor='black', boxstyle='square,pad=0.7'))


# Set figure label: A and axis/title
ax[0,0].text(-7, 1.1, 'A', fontsize=fslabel)
ax[0,0].set_xlabel('Time (months)')
ax[0,0].set_ylabel("Progression Free Survival Probability")
ax[0,0].title.set_text('Atezolizumab plus Bevacizumab arm') # 'ATEZO/BEVA: PARP1-HIGH/Low and PBRM1-mut/wt'

# Set figure label: C and axis/title
ax[1,0].text(-7, 1.05, 'C', fontsize=fslabel)
ax[1,0].set_xlabel('Time (months)')
ax[1,0].set_ylabel("Progression Free Survival Probability")
ax[1,0].title.set_text('Atezolizumab plus Bevacizumab arm') # 'ATEZO/BEVA: PARP1-HIGH/Low and PBRM1-mut/wt'

# Set figure label: D and axis/title
ax[1,1].text(-4, 1.05, 'D', fontsize=fslabel)
ax[1,1].set_xlabel('Time (months)')
ax[1,1].set_ylabel("Progression Free Survival Probability")
ax[1,1].title.set_text('Sunitinib arm')

# Set figure label: B and axis/title
ax[0,1].text(-4, 1.1, 'B', fontsize=fslabel)
ax[0,1].set_xlabel('Time (months)')
ax[0,1].set_ylabel("Progression Free Survival Probability")
ax[0,1].title.set_text('Sunitinib arm')
# Set the ticks and ticklabels for all axes, plot and save the figure
ytickvals = [round(x,1) for x in list(np.arange(0,1.1,0.1))]
plt.setp(ax, yticks=ytickvals)
plt.savefig('/path/to/outdir/Figure_1.png', dpi = 1000)
cf = cf_aux.copy()
plt.close()
#plt.show()
