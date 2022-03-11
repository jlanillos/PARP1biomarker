# This script reproduces Figure 1, which include Kaplan-Meier curves for all treatment arms using PARP1 stratification and PBRM1 status
import pandas as pd
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
from lifelines.plotting import add_at_risk_counts
import numpy as np
from lifelines import CoxPHFitter

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


cf_aux = preparematrix(clinical_data_file, ega_expression_file, gene, group) # Get data matrix


# Figure plot function
def plotFig(cf_aux, fig, ax, column, arm, gene, group, letterkm, letterrisk, nr_decimals):
    addriskyposition = 0.6
    fslabel = 18
    fsylabel = 10
    fsxlabel = 10
    xfs = 9
    ax[1,column].set_visible(False)
    ax[3,column].set_visible(False)
    print(arm)
    if arm == 'Entire cohort':
        cf = cf_aux.copy()
        titletext = 'Entire cohort'
    else:
        cf = cf_aux.loc[cf_aux['ARM'] == arm]
        titletext = arm.replace('+B', ' plus b')
    #### Mean Z-score (PARP1 HIGH/ LOW) - Censored and Not Censored
    kmf = KaplanMeierFitter()
    T = cf['PFS_MONTHS'].loc[cf[group] == 1]
    E = cf['PFS_CENSOR'].loc[cf[group] == 1]
    kmf.fit(T, event_observed=E, label="PARP1-high")
    kmf.plot_survival_function(ax=ax[0,column], show_censors = True, ci_show=False, color='k', marker='x', markersize = 1)
    kmfaux = kmf

    kmf = KaplanMeierFitter()
    T1 = cf['PFS_MONTHS'].loc[cf[group] == 0]
    E1 = cf['PFS_CENSOR'].loc[cf[group] == 0]
    kmf.fit(T1, event_observed=E1, label="PARP1-low")
    kmf.plot_survival_function(ax=ax[0,column], show_censors = True, ci_show=False, color='gold', marker='x', markersize = 1)
    ax[0,column].legend(fontsize=9, handlelength=1)
    add_at_risk_counts(kmfaux, kmf, ax=ax[1,column], rows_to_show=['At risk'], ypos=addriskyposition)#, size=8)
    #ax[0,0].set_xlabel("At risk", fontsize = xfs, ha = "left")

    # CoxPHFitter - PARP1 Univariate
    cph = CoxPHFitter()
    df_r = cf[[group,'PFS_MONTHS','PFS_CENSOR']]
    df_dummy = pd.get_dummies(df_r, drop_first=True)
    colcov = group + '_1'
    cph.fit(df_dummy, 'PFS_MONTHS', event_col='PFS_CENSOR', formula= colcov)
    comp = 'PARP1-univ'
    hr = np.round(cph.summary.T[colcov]['exp(coef)'], decimals=2)
    ci = str(np.round(cph.summary.T[colcov]['exp(coef) lower 95%'], decimals=2)) + '-' + str(np.round(cph.summary.T[colcov]['exp(coef) upper 95%'], decimals=2))
    pval = cph.summary.T[colcov]['p']

    offset = 0.045
    pval = np.round(cph.summary.T[colcov]['p'], decimals = nr_decimals)
    fs_text = 10
    box_style=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=1, edgecolor='white')
    ax[0,column].text(0, 0+ offset, 'HR=' + str(hr) + ' (95% CI=' + ci  + '), ' + 'p=' +  str(pval), fontsize=fs_text, bbox=box_style) #dict(facecolor='none', edgecolor='black', boxstyle='square,pad=0.4')

    #### Mean Z-score (PARP1 HIGH/ LOW +  PBRM1 Mut/ WT) - CENSORED only
    if arm == 'Entire cohort':
        cf = cf_aux.copy()
    else:
        cf = cf_aux.loc[cf_aux['ARM'] == arm]
    # Keep only samples sequenced (there are samples that do not contain their mutational status)
    cf = cf.loc[~cf['FMI_SAMPLE_ID'].isnull()]
    T = cf['PFS_MONTHS'].loc[(cf[group] == 1) & (cf[gene + '_mut'] == 1)]
    E = cf['PFS_CENSOR'].loc[(cf[group] == 1) & (cf[gene + '_mut'] == 1)]
    T1 = cf['PFS_MONTHS'].loc[(cf[group] == 1) & (cf[gene + '_mut'] == 0)]
    E1 = cf['PFS_CENSOR'].loc[(cf[group] == 1) & (cf[gene + '_mut'] == 0)]

    T2 = cf['PFS_MONTHS'].loc[(cf[group] == 0) & (cf[gene + '_mut'] == 1)]
    E2 = cf['PFS_CENSOR'].loc[(cf[group] == 0) & (cf[gene + '_mut'] == 1)]
    T3 = cf['PFS_MONTHS'].loc[(cf[group] == 0) & (cf[gene + '_mut'] == 0)]
    E3 = cf['PFS_CENSOR'].loc[(cf[group] == 0) & (cf[gene + '_mut'] == 0)]

    kmf = KaplanMeierFitter()
    kmf.fit(T, event_observed=E, label="PARP1-high/" + gene + "mut")
    kmf.plot_survival_function(ax=ax[2,column], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)
    kmf1 = kmf

    kmf = KaplanMeierFitter()
    kmf.fit(T1, event_observed=E1, label="PARP1-high/" + gene + "wt")
    kmf.plot_survival_function(ax=ax[2,column], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)
    kmf2 = kmf

    kmf = KaplanMeierFitter()
    kmf.fit(T2, event_observed=E2, label="PARP1-low/" + gene + "mut")
    kmf.plot_survival_function(ax=ax[2,column], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)
    kmf3 = kmf

    kmf = KaplanMeierFitter()
    kmf.fit(T3, event_observed=E3, label="PARP1-low/" + gene + "wt")
    kmf.plot_survival_function(ax=ax[2,column], show_censors = True, ci_show=False, linestyle= 'solid',marker='x', markersize = 1)
    kmf4 = kmf

    ax[2,column].legend(fontsize=9, handlelength=1)

    add_at_risk_counts(kmf1, kmf2, kmf3, kmf4, ax=ax[3,column], rows_to_show=['At risk'], ypos=addriskyposition)#, size=8)
    ax[2,column].set_xlabel("At risk", fontsize = xfs, ha = "left")


    ytickvals = [round(x,1) for x in list(np.arange(0,1.1,0.1))]
    # Set figure label: A and axis/title
    plt.sca(ax[0,0])
    ax[0,column].text(-17, 1.15, letterkm, fontsize=fslabel)
    ax[0,column].set_xlabel('Time (months)', fontsize = fsxlabel)
    ax[0,column].set_ylabel("Progression Free Survival Probability", fontsize=fsylabel)
    ax[0,column].title.set_text(titletext)
    ax[0,column].set_yticks(ytickvals)

    # Set figure label: D and axis/title
    plt.sca(ax[1,0])
    ax[2,column].text(-17, 1.15, letterrisk, fontsize=fslabel)
    ax[2,column].set_xlabel('Time (months)', fontsize = fsxlabel)
    ax[2,column].xaxis.set_label_coords(.36, -.1)
    ax[2,column].set_ylabel("Progression Free Survival Probability", fontsize=fsylabel)
    ax[2,column].title.set_text(titletext)
    ax[2,column].set_yticks(ytickvals)


# Plot Figure
fig, ax = plt.subplots(4, 3, figsize=(16, 10), sharey=False, gridspec_kw={'height_ratios': [9, 1, 9, 1]})
plt.subplots_adjust(hspace=0.5,wspace=0.75,left=0.13, right=0.99, top=0.93, bottom=0.055)

arm = 'Atezolizumab+Bevacizumab' # ATEZOLIZUMAB plus BEVACIZUMAB TREATMENT ARM
# Add subplots
column = 0
letterkm = 'A'
letterrisk = 'D'
nr_decimals = 4
plotFig(cf_aux, fig, ax, column, arm, gene, group, letterkm, letterrisk, nr_decimals)

arm = 'Sunitinib' #SUNITINIB
# Add subplots
column = 1
letterkm = 'B'
letterrisk = 'E'
nr_decimals = 6
plotFig(cf_aux, fig, ax, column, arm, gene, group, letterkm, letterrisk, nr_decimals)

arm = 'Entire cohort' # ENTIRE COHORT
# Add subplots
column = 2
letterkm = 'C'
letterrisk = 'F'
nr_decimals = 8
plotFig(cf_aux, fig, ax, column, arm, gene, group, letterkm, letterrisk, nr_decimals)

plt.savefig('Figure_1NOBORDER_' + gene + '.png', dpi = 500)
plt.close()
