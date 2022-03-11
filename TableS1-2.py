# This script helps to calculate all the statistical analysis which ended up becoming Supplementary Tables 2 and 3
import pandas as pd
from lifelines.statistics import logrank_test
import numpy as np
from lifelines import CoxPHFitter
from scipy.stats import chi2_contingency

# Load data
clinical_data_file = 'mmc2(4).xlsx' # From Motzer et al, Cancer Cell 2020 (https://doi.org/10.1016/j.ccell.2020.10.011): Table S1. Please, remove first four metadata lines and keep table header
ega_expression_file = 'immotion151.filtered.uqTPMs.csv' # Output file from UQ-normalization.R
# Choose gene from which mutational status will be analyzed together with PARP1 expression
gene = 'PBRM1'
# Choose stratification methods
group = 'mediangroup'


# Data matrix function
def preparematrix(clinical_data_file, ega_expression_file, gene, group):
    c = pd.read_excel(clinical_data_file)
    df = pd.read_csv(ega_expression_file,sep=',')
    df.set_index('symbol', inplace=True)
    colnames = dict(zip(list(df.columns), list(df.columns.str.replace('.','-'))))
    df.rename(columns=colnames, inplace=True)
    df = df.apply(lambda x: 1000*x) # Multiply by 1000: https://www.biostars.org/p/162089/;https://www.biostars.org/p/106127/
    df = df.apply(lambda x: np.log2(x + 1)) # Converto back to log2(TPM)

    c = c.loc[c['HISTOLOGY_SARCOMATOID'].isin(['ccRCC_Sarc', 'ccRCC_nonSarc'])].copy() # Keep samples with ccRCC histology
    parp = df.loc['PARP1'] # Find PARP1 expression
    d = dict(zip(list(parp.index), list(parp.values)))
    c['PARP1'] = c['RNASEQ_SAMPLE_ID'].map(d)

    cf = c.loc[~c['PARP1'].isnull()] # cf will be a dataframe that contains both clinical and expression data

    cf[group] = '' # Stratify patients according PARP1 median value
    cf.loc[cf['PARP1'] <= cf['PARP1'].median(), group] = 0
    cf.loc[cf['PARP1'] > cf['PARP1'].median(), group] = 1

    # Binarize gene mutational status and prepare PARP1+gene info for each sample
    cf[gene + '_mut'] = np.where(cf['FMI_' + gene].isin(['loss','rearrangement', 'short-variant']), 1, 0)
    cf.loc[(cf['FMI_SAMPLE_ID'].isnull()),gene + '_mut']='NAN'
    cf[gene + '_mut_term'] = np.where(cf['FMI_' + gene].isin(['loss','rearrangement', 'short-variant']), gene + '-mut', gene + '-wt')
    cf.loc[(cf['FMI_SAMPLE_ID'].isnull()),gene + '_mut_term']='NAN'
    cf['PARP1_zgroup_term'] = np.where(cf[group] == 1, 'PARP1-High', 'PARP1-Low')
    cf['PARP1_' + gene] = cf[group].astype(str) + '_' + cf[gene + '_mut'].astype(str)
    cf['PARP1_' + gene + '_term'] = cf['PARP1_zgroup_term'].astype(str) + '/' + cf[gene + '_mut_term'].astype(str)
    cf['MSKCCPrognosticGroup'] = np.where(cf['MSKCC_RISK_SCORE'] == 'Low', 'Favorable', 'Int_Poor')
    # Save the initial dataframe into a safe copy (cf_aux)
    cf_aux = cf.copy()
    return cf_aux


# Get data matrix
cf_aux = preparematrix(clinical_data_file, ega_expression_file, gene, group)


def getStats(cf_aux, arm, logrank_stats, savepath):
    if arm == 'all':
        cf = cf_aux.copy()
        arm = 'Entire-cohort'
    else:
        cf = cf_aux.copy()
        cf = cf.loc[cf['ARM'] == arm]
    if arm == 'Atezolizumab+Bevacizumab':
        arm = 'AtezoBeva'
    # Mean Z-score (PARP1 HIGH/ LOW) - Censored and Not Censored
    group = 'mediangroup'
    T = cf['PFS_MONTHS'].loc[cf[group] == 1]
    E = cf['PFS_CENSOR'].loc[cf[group] == 1]
    T1 = cf['PFS_MONTHS'].loc[cf[group] == 0]
    E1 = cf['PFS_CENSOR'].loc[cf[group] == 0]
    results=logrank_test(T1,T,event_observed_A=E1, event_observed_B=E)
    comp = 'low vs high'
    logrank_stats = pd.DataFrame(columns=['ARM', 'comp', 'tscore', 'p'])
    logrank_stats = logrank_stats.append(pd.Series([arm, comp, results.test_statistic, results.p_value], index = ['ARM', 'comp', 'tscore', 'p']), ignore_index=True)

    # Mean Z-score (PARP1 HIGH/ LOW +  gene Mut/ WT) - CENSORED only
    cf = cf.loc[~cf['FMI_SAMPLE_ID'].isnull()] # Keep only samples sequenced (there are samples that do not contain their mutational status)
    T = cf['PFS_MONTHS'].loc[(cf[group] == 1) & (cf[gene + '_mut'] == 1)]
    E = cf['PFS_CENSOR'].loc[(cf[group] == 1) & (cf[gene + '_mut'] == 1)]
    T1 = cf['PFS_MONTHS'].loc[(cf[group] == 1) & (cf[gene + '_mut'] == 0)]
    E1 = cf['PFS_CENSOR'].loc[(cf[group] == 1) & (cf[gene + '_mut'] == 0)]

    T2 = cf['PFS_MONTHS'].loc[(cf[group] == 0) & (cf[gene + '_mut'] == 1)]
    E2 = cf['PFS_CENSOR'].loc[(cf[group] == 0) & (cf[gene + '_mut'] == 1)]
    T3 = cf['PFS_MONTHS'].loc[(cf[group] == 0) & (cf[gene + '_mut'] == 0)]
    E3 = cf['PFS_CENSOR'].loc[(cf[group] == 0) & (cf[gene + '_mut'] == 0)]

    # T2 is PARP1-low and gene-mut
    results=logrank_test(T2,T,event_observed_A=E2, event_observed_B=E)
    comp = 'low-mut vs high-mut'
    logrank_stats = logrank_stats.append(pd.Series([arm, comp, results.test_statistic, results.p_value], index = ['ARM', 'comp', 'tscore', 'p']), ignore_index=True)

    results=logrank_test(T2,T1,event_observed_A=E2, event_observed_B=E1)
    comp = 'low-mut vs high-wt'
    logrank_stats = logrank_stats.append(pd.Series([arm, comp, results.test_statistic, results.p_value], index = ['ARM', 'comp', 'tscore', 'p']), ignore_index=True)

    results=logrank_test(T2,T3,event_observed_A=E2, event_observed_B=E3)
    comp = 'low-mut vs low-wt'
    logrank_stats = logrank_stats.append(pd.Series([arm, comp, results.test_statistic, results.p_value], index = ['ARM', 'comp', 'tscore', 'p']), ignore_index=True)

    results=logrank_test(T,T1,event_observed_A=E, event_observed_B=E1)
    comp = 'high-mut vs high-wt'
    logrank_stats = logrank_stats.append(pd.Series([arm, comp, results.test_statistic, results.p_value], index = ['ARM', 'comp', 'tscore', 'p']), ignore_index=True)

    results=logrank_test(T1,pd.concat([T,T2,T3]),event_observed_A=E1, event_observed_B=pd.concat([E,E2,E3]))
    comp = 'high-wt vs others'
    logrank_stats = logrank_stats.append(pd.Series([arm, comp, results.test_statistic, results.p_value], index = ['ARM', 'comp', 'tscore', 'p']), ignore_index=True)

    # Cox Regression analysis
    covariates = ['AGE','SEX','MSKCCPrognosticGroup', 'LIVER_METASTASES', 'SARCOMATOID', 'PRIMARY_VS_METASTATIC']
    genecovariate = [gene + '_mut']

    # Stratified (mediangroup), univariate
    cph = CoxPHFitter()
    parp1var = 'Stratified'
    covariatetype = 'Univariate'
    group = 'mediangroup'
    df_r = cf[[group, 'PFS_MONTHS','PFS_CENSOR']]
    df_dummy = pd.get_dummies(df_r, drop_first=True)
    cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR', formula= group + '_1')
    cph.summary.to_csv(savepath + '_'.join(['Cox', 'A', arm, parp1var, covariatetype]) + '.csv', sep='\t')

    # Stratified (mediangroup), multivariate simple (MSKCC)
    cph = CoxPHFitter()
    parp1var = 'Stratified'
    MSKCCcovariates = 'MSKCCPrognosticGroup'
    group = 'mediangroup'
    df_r = cf[[group, 'PFS_MONTHS','PFS_CENSOR'] + ['MSKCCPrognosticGroup']]
    df_dummy = pd.get_dummies(df_r, drop_first=True)
    cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR', formula= group + '_1 + MSKCCPrognosticGroup_Int_Poor')
    cph.summary.to_csv(savepath + '_'.join(['Cox', 'B', arm, parp1var, MSKCCcovariates]) + '.csv', sep='\t')

    # Stratified (mediangroup), multivariate
    cph = CoxPHFitter()
    parp1var = 'Stratified'
    covariatetype = 'Multivariate'
    group = 'mediangroup'
    df_r = cf.loc[cf[gene + '_mut'] != 'NAN']
    df_r = df_r[[group, 'PFS_MONTHS','PFS_CENSOR'] + covariates + genecovariate]
    df_dummy = pd.get_dummies(df_r, drop_first=True)
    cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR', formula= group + '_1 + AGE + SEX_M + LIVER_METASTASES_yes + SARCOMATOID_yes + MSKCCPrognosticGroup_Int_Poor + PRIMARY_VS_METASTATIC_PRIMARY + ' + gene + '_mut_1')
    cph.summary.to_csv(savepath + '_'.join(['Cox', 'C', arm, parp1var, covariatetype, genecovariate[0]]) + '.csv', sep='\t')

    # Continuous, univariate
    cph = CoxPHFitter()
    parp1var = 'Continous'
    covariatetype = 'Univariate'
    group = 'PARP1'
    df_r = cf[[group, 'PFS_MONTHS','PFS_CENSOR']]
    df_dummy = pd.get_dummies(df_r, drop_first=True)
    cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR', formula= group)
    cph.summary.to_csv(savepath + '_'.join(['Cox', 'D', arm, parp1var, covariatetype]) + '.csv', sep='\t')

    # Continuous (mediangroup), multivariate simple (MSK)
    cph = CoxPHFitter()
    parp1var = 'Continous'
    MSKCCcovariates = 'MSKCCPrognosticGroup'
    group = 'PARP1'
    df_r = cf[[group, 'PFS_MONTHS','PFS_CENSOR'] + ['MSKCCPrognosticGroup']]
    df_dummy = pd.get_dummies(df_r, drop_first=True)
    cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR', formula= group + ' + MSKCCPrognosticGroup_Int_Poor')
    cph.summary.to_csv(savepath + '_'.join(['Cox', 'E', arm, parp1var, MSKCCcovariates]) + '.csv', sep='\t')

    # Continuous, multivariate
    cph = CoxPHFitter()
    parp1var = 'Continous'
    covariatetype = 'Multivariate'
    group = 'PARP1'
    df_r = cf.loc[cf[gene + '_mut'] != 'NAN']
    df_r = df_r[[group, 'PFS_MONTHS','PFS_CENSOR'] + covariates + genecovariate]
    df_dummy = pd.get_dummies(df_r, drop_first=True)
    cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR', formula= group + '+ AGE + SEX_M + LIVER_METASTASES_yes + SARCOMATOID_yes + MSKCCPrognosticGroup_Int_Poor + PRIMARY_VS_METASTATIC_PRIMARY + ' + gene + '_mut_1')
    cph.summary.to_csv(savepath + '_'.join(['Cox', 'F', arm, parp1var, covariatetype, genecovariate[0]]) + '.csv', sep='\t')

    return logrank_stats


savepath = '/STATS/'
# logrank stats
logrank_stats = pd.DataFrame(columns=['ARM', 'comp', 'tscore', 'p'])
logrank_stats_aux = logrank_stats.copy()
arm = 'Atezolizumab+Bevacizumab' # ATEZOLIZUMAB plus BEVACIZUMAB TREATMENT ARM
logrank_stats = getStats(cf_aux, arm, logrank_stats, savepath)
logrank_stats_aux = logrank_stats_aux.append(logrank_stats)


arm = 'Sunitinib' #SUNITINIB
logrank_stats = getStats(cf_aux, arm, logrank_stats, savepath)
logrank_stats_aux = logrank_stats_aux.append(logrank_stats)


arm = 'all' # ENTIRE COHORT
logrank_stats = getStats(cf_aux, arm, logrank_stats, savepath)
logrank_stats_aux = logrank_stats_aux.append(logrank_stats)

# Save results
logrank_stats = logrank_stats_aux.copy()
logrank_stats.to_csv(gene + '.logrank_stats_analysis.csv',sep='\t',index = None)


# We manually merged all statistical results (Cox-regression) derived from getStats function into Supplementary Table 3
# In the CLI (Bash):
#touch Atezo_ALL.csv; for i in Cox_*Atezo*; do echo $i >> Atezo_ALL.csv; cat $i >> Atezo_ALL.csv; done
#touch Sunitinib_ALL.csv; for i in Cox_*Sunitinib*; do echo $i >> Sunitinib_ALL.csv; cat $i >> Sunitinib_ALL.csv; done
#touch Entire-cohort_ALL.csv; for i in Cox_*Entire-cohort*; do echo $i >> Entire-cohort_ALL.csv; cat $i >> Entire-cohort_ALL.csv; done

########### SUPPLEMENTARY TABLE 2 ###########
####### Motzer's molecular subtypes #########
#############################################
## Analyzis of the Molecular subtypes described by Motzer et al in IMmotion151 and PARP1/GENE status
cf = cf_aux.copy()
cf = cf.loc[~cf['FMI_SAMPLE_ID'].isnull()]
angio = cf.groupby(['ARM','PARP1_' + gene + '_term','NMF_GROUP']).count()['RNASEQ_SAMPLE_ID'].reset_index()

dfangio = pd.DataFrame({'group':list(set(angio['NMF_GROUP']))})
for arm in list(set(angio['ARM'])):
    for group in list(set(angio['PARP1_' + gene + '_term'])):
        dff = angio.loc[(angio['ARM'] == arm) & (angio['PARP1_' + gene + '_term'] == group)]
        dfangio[arm + '_' + group] = dfangio['group'].map(dict(zip(list(dff['NMF_GROUP']), list(dff['RNASEQ_SAMPLE_ID']))))

for i in [col for col in list(dfangio.columns) if 'group' not in col]:
    dfangio['pct_' + i] = 100*(dfangio[i] / dfangio[i].sum())
dfangio.to_csv('molecular_subtype_tto_parp1AND' + gene + '.csv',sep='\t')

cf = cf_aux.copy()
cf = cf.loc[~cf['FMI_SAMPLE_ID'].isnull()]
angio = cf.groupby(['PARP1_' + gene + '_term','NMF_GROUP']).count()['RNASEQ_SAMPLE_ID'].reset_index()


for group in list(set(angio['PARP1_' + gene + '_term'])):
    dff = angio.loc[(angio['PARP1_' + gene + '_term'] == group)]
    dfangio[group] = dfangio['group'].map(dict(zip(list(dff['NMF_GROUP']), list(dff['RNASEQ_SAMPLE_ID']))))

groupscol = ['PARP1-High/PBRM1-wt','PARP1-High/PBRM1-mut','PARP1-Low/PBRM1-wt','PARP1-Low/PBRM1-mut']
dfangio = dfangio[['group'] + groupscol]

maincol = 'PARP1-High/PBRM1-wt'
cols2comp = ['PARP1-High/PBRM1-mut','PARP1-Low/PBRM1-wt','PARP1-Low/PBRM1-mut']

sums = dfangio[groupscol].sum(axis=0)

def computechi(dfangio, group, maincol, i):
    contigency = dfangio[[maincol, i]].loc[dfangio['group'] == group]
    contigency = contigency.append(pd.DataFrame(sums[[maincol, i]]).T)
    c, p, dof, expected = chi2_contingency(contigency)
    val = dfangio[i].loc[dfangio['group'] == group].values[0]
    countingroup = dfangio[[maincol] + cols2comp].loc[dfangio['group'] == group].sum(axis = 1).values[0]
    pct = '(' + str(int(np.round(100*(val/countingroup)))) + ')'
    val = str(val)
    if (p<0.00005):
        p = '***'
    elif (p>=0.00005) and (p<0.005):
        p = '**'
    elif  (p>=0.005) and (p<0.05):
        p = '*'
    else:
        p = ''
    result = val + ' ' + pct + p
    return(result)

for i in cols2comp:
    dfangio['comparison_' + i] = dfangio.apply(lambda x: computechi(dfangio, x['group'], maincol, i), axis = 1)

# Save to a CSV file which we manually modified to the final format as presented in Supplementary Table2
dfangio.to_csv('suppTable2.csv',sep='\t', index = None)
