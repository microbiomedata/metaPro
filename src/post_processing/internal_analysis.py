import pandas as pd
import numpy as np
import os
import re
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import minimize

__author__ =  "Gao, Yuqian <Yuqian.Gao@pnnl.gov>"
__maintainer__ = "Anubhav <anubhav@pnnl.gov>"

class InternalAnalysis:

    def __init__(self, parent_folder,):
        self.parent_folder = parent_folder


    def findproteinname(s):
        '''Get Protein type
        :param s:
        :return:
        '''
        p1 = re.compile(r'^Contaminant_')
        p2 = re.compile(r'^XXX_Contaminant_')
        p3 = re.compile(r'^XXX_')
        if p1.search(s) is not None:
            return 'None'
        elif p2.search(s) is not None:
            return 'None'
        elif p3.search(s) is not None:
            return 'Reversed'
        else: # the rest
            return 'Forward'


    def cleansequence(s):
        '''
        clean peptide sequence is the sequence without prefix and postfix but with oxidation
        :param s:
        :return:
        '''
        p = re.compile(r'\.(?P<cleanseq>[A-Z\*@#]+)\.')
        m = p.search(s)
        return m.group('cleanseq')
    
    # Find # of missed cleavages from clean peptide
    def MissedCleavages(s):
        mc = 0
        for i in range(len(s)-1):
            if (s[i]=="K" or s[i]=="R") and (s[i+1]!="P"):
                mc = mc +1
        return mc

    def process_data(self):
        '''

        :return:
        '''
        # import table
        data = pd.read_table(self.parent_folder +'resultants_df.txt')[[
                           'JobNum', 'Dataset_x', 'Dataset_y', 'Scan', 'Protein', 'Peptide', \
                           'NTT', 'DelM', 'DelM_PPM', 'StatMomentsArea', 'PeakMaxIntensity', \
                           'MSGFDB_SpecEValue', 'EValue', 'QValue', 'PepQValue']]

        # Add protein type, clean peptide sequence, specID
        data.rename(columns={'Dataset_x': 'Dataset', 'Dataset_y': 'Dataset_ID', 'JobNum': 'Job'}, inplace=True)
        data['Protein_Type'] = data['Protein'].apply(self.findproteinname)
        data['Clean Peptide Sequence'] = data['Peptide'].apply(self.cleansequence)
        data["SpecID"] = data.apply(lambda row: str(row["Dataset_ID"]) + "_" + str(row["Scan"]), axis=1)

        # Save protein to peptide mapping of forward peptide identification
        df_mapping = data[data['Protein_Type'] == 'Forward'] [['Protein', 'Peptide', 'Clean Peptide Sequence']].copy()
        df_mapping.drop_duplicates(inplace=True)

        # Calculate the redundancy by clean peptide sequence and protein_type
        df_new = df_mapping[['Protein', 'Clean Peptide Sequence']].copy()
        df_new.drop_duplicates(inplace=True)
        df_redundancy = df_new.groupby(['Clean Peptide Sequence']).count()
        df_redundancy.reset_index(inplace=True)
        df_redundancy.rename(columns={'Protein': 'Clean Peptide Sequence Redundancy'}, inplace=True)

        # merge redundancy to protein peptide map
        df_mapping = df_mapping.merge(df_redundancy, how='left', on=['Clean Peptide Sequence'])
        del df_redundancy

        # save to file
        df_mapping.to_csv("Results/protein_peptide_map.csv", index=False)
        del df_mapping

        # save the map of dataset id, dataset name and job id
        df_ids = data[['Dataset', 'Job', 'Dataset_ID']].copy()
        df_ids.drop_duplicates(inplace=True)
        df_ids.to_csv("Results/dataset_job_map.csv", index=False)
        del df_ids

        # Remove 'Protein', 'dataset_name', 'job_id'
        del data['Protein']
        del data['Dataset']
        del data['Job']
        data.drop_duplicates(inplace=True)  # this is so important to remove any duplicated rows

        # Make a table for just forward peptides
        data_cleaned_forward = data[data['Protein_Type'] == 'Forward'].copy()
        del data_cleaned_forward['Protein_Type']
        data_cleaned_forward.drop_duplicates(inplace=True)

        # Make a table for just reversed peptides
        data_cleaned_reversed = data[data['Protein_Type'] == 'Reversed'].copy()
        del data_cleaned_reversed['Protein_Type']
        data_cleaned_reversed.drop_duplicates(inplace=True)

        # remove data file
        del data

        # Dataset list
        dataset_list = data_cleaned_forward['Dataset_ID'].unique().tolist()

        # Export forward and reverse for individual dataset
        for i in range(len(dataset_list)):
            df_ff = data_cleaned_forward[data_cleaned_forward['Dataset_ID'] == dataset_list[i]].copy()
            df_ff.to_csv("Results/Data/" + str(dataset_list[i]) + "_forward_peptide_identification.csv", \
                         index=False)
            df_rr = data_cleaned_reversed[data_cleaned_reversed['Dataset_ID'] == dataset_list[i]].copy()
            df_rr.to_csv("Results/Data/" + str(dataset_list[i]) + "_reversed_peptide_identification.csv", \
                         index=False)
            del df_ff
            del df_rr

        del data_cleaned_forward
        del data_cleaned_reversed


# optimize the filtering criteria and filter the data
def parameter_optimization(dataset_ID):
    '''

    :param dataset_ID:
    :return:
    '''
    data_f = pd.read_csv("Results/Data/" + str(dataset_ID) + \
                         "_forward_peptide_identification.csv")
    data_r = pd.read_csv("Results/Data/" + str(dataset_ID) + \
                         "_reversed_peptide_identification.csv")
    
    # Only fully tryptic and maximum of 2 missed cleavages
    data_f["MissedCleavage"] = data_f["Clean Peptide Sequence"].apply(MissedCleavages)
    data_r["MissedCleavage"] = data_r["Clean Peptide Sequence"].apply(MissedCleavages)
    data_f = data_f[(data_f["NTT"]==2)&(data_f["MissedCleavage"]<=2)].copy()
    data_r = data_r[(data_r["NTT"]==2)&(data_r["MissedCleavage"]<=2)].copy()

    # Fit a 1-D data of DelM_PPM and get the peak ppm_shift
    ppm_shift_df = data_f[(data_f["DelM_PPM"] < 10) & (data_f["DelM_PPM"] > -10)].copy()
    ax = sns.distplot(ppm_shift_df["DelM_PPM"])
    ppm_shift = ax.get_lines()[0].get_xdata()[np.argmax(ax.get_lines()[0].get_ydata())]

    def PepFDR(Params):
        '''

        :param Params:
        :return:
        '''
        # function to minimize
        delppm1, delppm2, log10_specprob = Params  # use log10 value so that it is managable for the computer
        ### The FDR function ###
        pep_total = data_r['Clean Peptide Sequence'].unique().size
        df_r = data_r[(data_r["DelM_PPM"] < ppm_shift + delppm1) & (data_r["DelM_PPM"] > ppm_shift - delppm2) & \
                      (data_r["MSGFDB_SpecEValue"] < 10 ** log10_specprob)].copy()
        df_f = data_f[(data_f["DelM_PPM"] < ppm_shift + delppm1) & (data_f["DelM_PPM"] > ppm_shift - delppm2) & \
                      (data_f["MSGFDB_SpecEValue"] < 10 ** log10_specprob)].copy()
        f_pep = df_f['Clean Peptide Sequence'].unique().size
        r_pep = df_r['Clean Peptide Sequence'].unique().size
        ### fdr_pep ###
        if (f_pep == 0) & (r_pep == 0):
            fdr_pep = 1
        else:
            fdr_pep = r_pep / (f_pep + r_pep)
        return 1 / (0.050001 - fdr_pep) * (-f_pep)
        # Note:
        # in this function, if fdr_pep>0.01, it returns +,
        # if fdr_pep < 0.01, it returns -, the closest to 0.01, the smaller the return

    # Add constraint
    def constraint1(Params):
        delppm1, delppm2, log10_specprob = Params  # use log10 value so that it is managable for the computer
        return 20 - delppm1 - delppm2

    def constraint2(Params):
        delppm1, delppm2, log10_specprob = Params  # use log10 value so that it is managable for the computer
        return delppm1 - 5

    def constraint3(Params):
        delppm1, delppm2, log10_specprob = Params  # use log10 value so that it is managable for the computer
        return delppm2 - 5

    con1 = {'type': 'ineq', 'fun': constraint1}
    con2 = {'type': 'ineq', 'fun': constraint2}
    con3 = {'type': 'ineq', 'fun': constraint3}

    # Miminize PepFDR
    initial_guess = [min(10, max(data_f["DelM_PPM"]) - ppm_shift), \
                     min(10, ppm_shift - min(data_f["DelM_PPM"])), \
                     -15]
    # print(dataset_ID, initial_guess)
    result = minimize(PepFDR, initial_guess, method='COBYLA', constraints=[con1, con2, con3])
    if result.success:
        fitted_params = result.x
    else:
        n = 0
        while (result.success == False) & (n < 40):
            initial_guess = [initial_guess[0] + 0.2, initial_guess[1] + 0.2, initial_guess[2]]
            result = minimize(PepFDR, initial_guess, method='COBYLA', constraints=[con1, con2, con3])
            n = n + 1
        if result.success:
            fitted_params = result.x
        else:
            print("Failed optimization @ Dataset ID: %s, Database: %s" % (dataset_ID, DB))
            print(initial_guess)
            raise ValueError(result.message)

    ### plots ###
    plt.close('all')
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(7, 15))
    # plot the DelM_PPM
    sns.distplot(data_f["DelM_PPM"], label='forward', ax=axes[0], bins=200)
    axes[0].axvline(x=ppm_shift, color='r', linestyle='--')
    axes[0].axvline(x=ppm_shift - fitted_params[1], color='g', linestyle='-')
    axes[0].axvline(x=ppm_shift + fitted_params[0], color='g', linestyle='-')
    sns.distplot(data_r["DelM_PPM"], label='reversed', ax=axes[0], bins=200)
    axes[0].set_xlabel('DelM_PPM')
    axes[0].set_ylabel('Density')
    axes[0].set_title(r'left bound: %.2f; right bound: %.2f' % \
                      (ppm_shift - fitted_params[1], ppm_shift + fitted_params[0]))
    # plot MSGFDB_SpecEValue
    sns.distplot(data_f['MSGFDB_SpecEValue'], label='forward', kde=False, ax=axes[1], bins=5000)
    sns.distplot(data_r['MSGFDB_SpecEValue'], label='reversed', kde=False, ax=axes[1], bins=5000)
    axes[1].axvline(x=10 ** fitted_params[2], color='g', linestyle='-')
    axes[1].set_title('MSGFDB_SpecEValue: %.2e' % (10 ** fitted_params[2]))
    axes[1].set_xlabel("log(MSGFDB_SpecEValue)")
    axes[1].set_ylabel('Density')
    axes[1].set_xscale('log')
    # save
    plt.tight_layout()
    plt.savefig("Results/Plots/" + str(dataset_ID) + ".jpg")
    plt.close('all')

    # #Filter the data
    df_f = data_f[
        (data_f["DelM_PPM"] < ppm_shift + fitted_params[0]) & (data_f["DelM_PPM"] > ppm_shift - fitted_params[1]) & \
        (data_f["MSGFDB_SpecEValue"] < 10 ** fitted_params[2])].copy()
    del data_f

    df_r = data_r[
        (data_r["DelM_PPM"] < ppm_shift + fitted_params[0]) & (data_r["DelM_PPM"] > ppm_shift - fitted_params[1]) & \
        (data_r["MSGFDB_SpecEValue"] < 10 ** fitted_params[2])].copy()
    del data_r

    # Calculate FDR
    f_spec = df_f['SpecID'].unique().size  # Modified on Apr23 to only count for unique dataset-scan
    r_spec = df_r['SpecID'].unique().size
    if (f_spec == 0) & (r_spec == 0):
        fdr_spec = 1
    else:
        fdr_spec = r_spec / (f_spec + r_spec)

    # Modified on Aug 02, 2019 to use Clean Peptide Sequence instead of peptide
    f_pep = df_f['Clean Peptide Sequence'].unique().size
    r_pep = df_r['Clean Peptide Sequence'].unique().size
    if (f_pep == 0) & (r_pep == 0):
        fdr_pep = 1
    else:
        fdr_pep = r_pep / (f_pep + r_pep)

    ######################################Dataset FDR_table################################################
    df_Metadata = pd.DataFrame({'PPM center': ppm_shift, 'PPM cutoff left': ppm_shift - fitted_params[1], \
                                'PPM cutoff right': ppm_shift + fitted_params[0], \
                                'SpecProb cutoff': 10 ** fitted_params[2], \
                                'Spectra forward': f_spec, 'Spectra reverse': r_spec, \
                                'Peptide forward': f_pep, 'Peptide reverse': r_pep, \
                                'PeptideFDR': [fdr_pep], 'SepctraFDR': [fdr_spec], 'Dataset_ID': dataset_ID, \
                                'Initial ppm cutoff left': ppm_shift - initial_guess[1], \
                                'Initial ppm cutoff right': ppm_shift + initial_guess[0], \
                                'Initial SpecProb cutoff': 10 ** initial_guess[2]})

    ######################################SpecID table################################################
    df_SpecID_f = df_f[['SpecID', 'Dataset_ID', 'Scan', 'Peptide', 'Clean Peptide Sequence', 'MSGFDB_SpecEValue', \
                        'StatMomentsArea', 'DelM_PPM']].copy()
    df_SpecID_f.drop_duplicates(inplace=True)

    df_SpecID_r = df_r[['SpecID', 'Dataset_ID', 'Scan', 'Peptide', 'Clean Peptide Sequence', 'MSGFDB_SpecEValue', \
                        'StatMomentsArea', 'DelM_PPM']].copy()
    df_SpecID_r.drop_duplicates(inplace=True)

    ######################################Spetra Count table##########################################
    df_spec = df_f[['Clean Peptide Sequence', 'SpecID']].copy()
    df_spec.drop_duplicates(inplace=True)
    df_SpectraCount = df_spec.groupby(['Clean Peptide Sequence']).size().rename('Peptide Spectra Count').reset_index()
    del df_spec
    df_SpectraCount['Dataset_ID'] = dataset_ID

    ######################################Spetra Count table##########################################
    df_Intensity = df_f[['Clean Peptide Sequence', 'StatMomentsArea']].groupby(['Clean Peptide Sequence']).max()
    df_Intensity.rename(columns={'StatMomentsArea': 'Peptide Peak Area'}, inplace=True)
    df_Intensity.reset_index(inplace=True)
    df_Intensity['Dataset_ID'] = dataset_ID

    return df_Metadata, df_SpecID_f, df_SpecID_r, df_SpectraCount, df_Intensity
