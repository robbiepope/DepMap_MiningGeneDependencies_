# Import packages
import time
import pandas as pd
import networkx as nx
from DepMapTools.Correlations import GeneCorrelations
from DepMapTools.Networks import NetworkAnalysis
from DepMapTools.Communities import CommunityAnalysis
from DepMapTools.GeneOntology import OntologyAnalysis
# Ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
#%%
# SINGLE GENE OF INTEREST ANALYSIS


class GOIAnalysis:
    """
    Class for single gene correlation analysis and network creation
    """

    def __init__(self, df, gene):
        """
        Constructor class for GOIAnalysis class

        :param df: (DF) Achilles DF of gene essentiality scores
        :param gene: (STR) Name of gene of interest
        """
        self.df = df
        self.gene = gene
        self.gc = GeneCorrelations()
        self.ntwk = NetworkAnalysis()
        self.ca = CommunityAnalysis()
        self.oa = OntologyAnalysis()

    def calculate_correlation(self):
        """
        Calculate pearson and spearman correlation coefficients and mean difference p-value for gene of interest

        :return df_pr: (DF) DF of pearson correlation coefficients
        :return df_sr: (DF) DF of spearman correlation coefficients
        :return df_dm: (DF) DF of mean differences
        """

        # Compute Pearson, Spearman and Mean differences for specified gene
        df_pr = self.gc.pearson(self.df, self.gene)
        df_sr = self.gc.spearman(self.df, self.gene)
        df_dm = self.gc.diff_means(self.df, self.gene)
        return df_pr, df_sr, df_dm

    @staticmethod
    def merge_dfs(df_pr, df_sr, df_dm, threshold):
        """
        Concatenate by column (axis=1) the top results of the individual correlation DFs for a gene of interest

        :param df_pr: (DF) DF of pearson correlation coefficients for a gene of interest
        :param df_sr: (DF) DF of spearman correlation coefficients for a gene of interest
        :param df_dm: (DF) DF of mean differences for a gene of interest
        :param threshold: (INT) Number of top correlations from DFs to merge (e.g., 400)
        :return df_con2: (DF) DF of all three DFs concatented on axis 1
        """

        # Concatenate Pearson and Spearman correlations then merge with mean difference
        df_con1 = pd.concat([df_pr.iloc[0:threshold, :],
                             df_sr.iloc[0:threshold, :]],
                            axis=1
                            )
        df_con2 = pd.concat([df_con1,
                             df_dm.iloc[0:threshold, :]],
                            axis=1
                            )
        df_con2.columns = ['gene_p', 'Pearson_r', 'Pearson_p',
                           'gene_s', 'Spearman_r', 'Spearman_p',
                           'gene_md', 'Mean_Diff', 'Mean_Diff_p'
                           ]
        return df_con2

    def correlation_df(self, df_merge, top_genes=100):
        """
        Compute the overlap between correlation measures and find the top common genes across the three measures
        for a gene of interest.

        :param df_merge: (DF) Concatenated DF of all three correlation measures for a gene of interest
        :param top_genes: (INT) Top number of gene correlations to return in the DF, default=100
        :return df_corr: (DF) DF of top correlations for a gene of interest
        """

        # Calculate the overlap between correlation measures
        pr_dict = dict(zip(df_merge['gene_p'],
                           df_merge['Pearson_r'])
                       )
        sr_dict = dict(zip(df_merge['gene_s'],
                           df_merge['Spearman_r'])
                       )
        md_dict = dict(zip(df_merge['gene_md'],
                           df_merge['Mean_Diff'])
                       )
        # Compute the common keys present in all three dictionaries
        dict_all = self.gc.common_key_3(pr_dict,
                                        sr_dict,
                                        md_dict
                                        )
        # Convert to pandas DF
        df_all = pd.DataFrame.from_dict(dict_all,
                                        orient='index',
                                        dtype=None
                                        )
        df_all = df_all.reset_index()
        df_all.columns = [self.gene + '_Correlations',
                          self.gene + '_Pearson_R',
                          self.gene + '_Spearman_R',
                          self.gene + '_Mean_Diff'
                          ]
        df_all[self.gene + '_Pearson_p'] = df_merge['Pearson_p']
        df_all[self.gene + '_Spearman_p'] = df_merge['Spearman_p']
        df_all[self.gene + '_Mean_Diff_p'] = df_merge['Mean_Diff_p']
        # Create DF of correlations for gene
        df_corr = pd.DataFrame()
        df_corr = pd.concat([df_corr,
                             df_all.iloc[0:top_genes, :]],
                            axis=1
                            )
        return df_corr

    def generate_network(self, df_corr):
        """
        Generate the STRING networks for a gene of interest using the top correlated genes

        :param df_corr: (DF) DF of top correlated genes to generate the STRING based network
        :return g: (G) Netowrkx undirected graph
        :return  interactions: (DF) Dataframe with network data + weights as STRING combined score
        """

        interactions = self.ntwk.generate_interactions(df_corr.iloc[:, 0].to_list())
        g = self.ntwk.network_import(interactions)
        g.remove_edges_from(nx.selfloop_edges(g))
        return g, interactions

    def goi_analysis(self, threshold):
        """
        Compute correlations and generate STRING network for a gene of interest

        :param threshold: (INT) Number of top correlations from DFs to merge (e.g., 400)
        :return result_dict: (DICT) Dictionary of results with the following keys + values:
                              Correlation: (DF) Correlation DF
                              Interaction: (DF) Dataframe with network data + weights as STRING combined score
                              Network: (G) Networkx undirected graph
        """

        result_dict = {}
        df_pr, df_sr, df_dm = self.calculate_correlation()
        df_merge2 = self.merge_dfs(df_pr, df_sr, df_dm, threshold)
        df_corr = self.correlation_df(df_merge2)
        g, interactions = self.generate_network(df_corr)
        time.sleep(1)
        result_dict['Correlation'] = df_corr
        result_dict['Interaction'] = interactions
        result_dict['Network'] = g
        return result_dict
