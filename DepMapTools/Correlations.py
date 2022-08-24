# Import packages
import pandas as pd
import numpy as np
from scipy import stats
# Ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
#%%
# CREATE DF OF CORRELATIONS


class GeneCorrelations:
    """
    Class to calculate gene dependency correlations
    """

    @staticmethod
    def pearson(df, gene):
        """
        Calculate Pearson's Correlation Coefficient for a gene of interest (GOI) and all other genes in the DF
        with Scipy.Stats pearsonr

        :param df: (DF) DF of gene essentiality scores
        :param gene: (STR) Gene of interest name
        :return df_pr: (DF) Pandas DF of Pearson gene dependency correlations for the GOI
        """
        dict1 = {'Gene': [],
                 'Pearson_r': [],
                 'Pearson_p': []
                 }
        for col in df:
            (r, p) = stats.pearsonr(df[gene],
                                    df[col]
                                    )
            dict1['Gene'].append(col)
            dict1['Pearson_r'].append(r)
            dict1['Pearson_p'].append(p)
        df_pr = pd.DataFrame.from_dict(dict1)
        df_pr = df_pr.sort_values('Pearson_p')
        df_pr.reset_index(drop=True)
        return df_pr

    @staticmethod
    def spearman(df, gene):
        """
        Calculate Spearman's Rank Correlation Coefficient for a gene of interest (GOI)
        and all other genes in the DF with Scipy.Stats spearmanr

        :param df: (DF) DF of gene essentiality scores
        :param gene: (STR) Gene of interest name
        :return df_sr: (DF) Pandas DF of Spearman gene dependency correlations for the GOI
        """

        dict1 = {'Gene': [],
                 'Spearman_r': [],
                 'Spearman_p': []
                 }
        for col in df:
            (r, p) = stats.spearmanr(df[gene],
                                     df[col]
                                     )
            dict1['Gene'].append(col)
            dict1['Spearman_r'].append(r)
            dict1['Spearman_p'].append(p)
        df_sr = pd.DataFrame.from_dict(dict1)
        df_sr = df_sr.sort_values('Spearman_p')
        df_sr.reset_index(drop=True)
        return df_sr

    @staticmethod
    def diff_means(df, gene):
        """
        Calculate the statistical significance (independent T-Test) of the difference between the mean essentiality
        score of the top 50 (sensitive) genes and the bottom 50 (resistant) genes for a GOI
        :param df: (DF) DF of gene essentiality scores
        :param gene: (STR) Gene of interest name
        :return diff_mean: (DF) Pandas DF of p-values for difference in means for the GOI
        """

        # Select Achilles data on GOI
        df1 = df.loc[:, gene]
        df1 = df1.to_frame().reset_index()
        df1.columns = ['CCLE Name',
                       'gene_dep_score'
                       ]
        # Sort gene dependency scores from lowest (sensitive) to highest (resistant)
        df1 = df1.sort_values(by='gene_dep_score')
        df_res = df1.iloc[-50:]
        df_sen = df1.iloc[0:50]
        df = df.reset_index()
        df['sensitive'] = df['CCLE_Name'].isin(df_sen['CCLE Name']).astype(np.int8)
        df['resistant'] = df['CCLE_Name'].isin(df_res['CCLE Name']).astype(np.int8)
        # Compute the difference of means: diff_means_exp
        diff_mean = df[df['resistant'] == 1].iloc[:, 1:-2].mean() - df[df['sensitive'] == 1].iloc[:, 1:-2].mean()
        diff_mean = pd.DataFrame(diff_mean)
        diff_mean.columns = ['Mean_Diff']
        diff_mean.index.name = 'Gene'
        # Compute independent T-Test to calculate statistical significance (p)
        t, p = stats.ttest_ind(df[df['sensitive'] == 1].iloc[:, 1:-2],
                               df[df['resistant'] == 1].iloc[:, 1:-2],
                               nan_policy='omit')
        diff_mean['Mean_Diff_p'] = p
        diff_mean = diff_mean.sort_values('Mean_Diff_p')
        diff_mean = diff_mean.reset_index()
        return diff_mean

    @staticmethod
    def common_key_3(dict_a, dict_b, dict_c):
        """
        Find the genes (common keys) present in three dictionaries of gene essentiality correlations

        :param dict_a: (DICT) Dictionary of genes (keys) and correlations (values) for a GOI
        :param dict_b: (DICT) Dictionary of genes (keys) and correlations (values) for a GOI
        :param dict_c: (DICT) Dictionary of genes (keys) and correlations (values) for a GOI
        :return dict1: (DICT) Dictionary of common genes (keys) and correlations (values) for a GOI
        """

        list1 = []
        dict1 = {}
        for i in dict_a.keys():
            for j in dict_b.keys():
                for k in dict_c.keys():
                    if i == j and j == k:
                        list1.append(i)
        for k in list1:
            dict1.update({k: [dict_a[k],
                              dict_b[k],
                              dict_c[k]
                              ]
                          }
                         )
        return dict1
