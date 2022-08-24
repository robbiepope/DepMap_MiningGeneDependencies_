# Import packages
import requests
import random
import time
import pandas as pd
import numpy as np
import networkx as nx
from math import comb
from statsmodels.stats import multitest
from DepMapTools.Correlations import GeneCorrelations
from DepMapTools.GeneOntology import OntologyAnalysis
# Ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
#%%
# GENERATE STRING BASED NETWORK FOR CORRELATED GENES


class NetworkAnalysis:
    """
    Class to generate STRING based networks and centrality measures
    """

    @staticmethod
    def generate_interactions(gene_list):
        """
        Request STRING network data via STRING API for list of genes and store as pandas DF

        :param gene_list: (LIST of STR) List of genes generated from correlation analysis
        :return interactions: (DF) Dataframe with network data + weights as STRING combined score
        """

        # Construct URL + set parameters
        string_api_url = "https://string-db.org/api"
        output_format = "tsv-no-header"
        method = "network"
        request_url = "/".join([string_api_url, output_format, method])
        params = {"identifiers": "%0d".join(gene_list),
                  "species": 9606,
                  "caller_identity": "http://www.sussex.ac.uk/lifesci/hocheggerlab/",
                  "add_nodes": 0
                  }

        # Call STRING + re-format as pandas DF
        r = requests.post(request_url, data=params)
        cols = ['stringId_A', 'stringId_B', 'preferredName_A', 'preferredName_B', 'ncbiTaxonId', 'score', 'nscore',
                'fscore', 'pscore', 'ascore', 'escore', 'dscore', 'tscore'
                ]

        # Pull the text from the response object and split based on new lines + split each line into its components
        lines = r.text.split('\n')
        data = [ln.split('\t') for ln in lines]

        # Convert to dataframe using the first row as the column names; drop empty, final row
        # Add unconnected nodes from original list
        df = pd.DataFrame(data[1:-1], columns=cols)
        for entry in gene_list:
            if entry not in df['preferredName_A']:
                df.loc[len(df.index)] = ['NaN', 'NaN', entry, entry, 'NaN', 'NaN', 'NaN',
                                         'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'
                                         ]

        # DF with the preferred names of the two proteins and the combined score of the interaction
        interactions = df[['preferredName_A',
                           'preferredName_B',
                           'score'
                           ]]
        return interactions

    @staticmethod
    def network_import(df):
        """
        Generate networkx graph using interaction data imported from STRING database

        :param df: (DF) DF of protein interactions imported from STRING database.
                   Column[0] = Protein A (node)
                   Column[1] = Protein B (node)
                   Column[2] = Combined confidence score of interaction (weight)
        :return g: (G) Networkx undirected graph with weighted edges
        """

        g = nx.Graph(name='Protein Interaction Graph')
        interactions = np.array(df)
        for i in range(len(interactions)):
            interaction = interactions[i]
            a = interaction[0]
            b = interaction[1]
            w = float(interaction[2])
            g.add_weighted_edges_from([(a, b, w)])
        return g


#%%
# CENTRALITY ANALYSIS FOR FUNCTIONAL INTERACTION NETWORKS


class CentralityAnalysis:
    """
    Class to calculate network centrality measures
    """

    def __init__(self):
        """
        Constructor for CentralityAnalysis Class
        """

        self.oa = OntologyAnalysis()

    @staticmethod
    def harmonic_centrality(g):
        """
        Calculate normalised harmonic centrality for a network

        :param g: (G) Networkx undirected graph
        :return har: (DICT) Dictionary of nodes (keys) and harmonic centrality measure as values
        """

        har = nx.harmonic_centrality(g)
        # Normalise results
        har = {k: v / len(g) - 1 for k, v in har.items()}
        return har

    @staticmethod
    def degree_centrality(g):
        """
        Calculate degree centrality for a network

        :param g: (G) Networkx undirected graph
        :return deg_cen: (DICT) Dictionary of nodes (keys) and degree centrality measure as values
        """

        deg_cen = nx.degree_centrality(g)
        return deg_cen

    @staticmethod
    def closeness(g):
        """
        Calculate closeness centrality for a network

        :param g: (G) Networkx undirected graph
        :return clo: (DICT) Dictionary of nodes (keys) and closeness centrality measure as values
        """

        clo = nx.closeness_centrality(g)
        return clo

    @staticmethod
    def eigen_centrality(g, iterations):
        """
        Calculate eigen centrality for a network

        :param g: (G) Networkx undirected graph
        :param iterations: (INT) Number of iterations to calculate eigen centrality
        :return eig: (DICT) Dictionary of nodes (keys) and eigen centrality measure as values
        """

        eig = nx.eigenvector_centrality(g, max_iter=iterations)
        return eig

    @staticmethod
    def betweenness_centrality(g):
        """
        Calculate betweenness centrality for a network

        :param g: (G) Networkx undirected graph
        :return bet: (DICT) Dictionary of nodes (keys) and betweenness centrality measure as values
        """

        bet = nx.betweenness_centrality(g)
        return bet

    @staticmethod
    def ntwk_desnity(g):
        """
        Calculate density of network

        :param g: (G) Networkx undirected graph
        :return: (INT) Density measure of network
        """

        return nx.density(g)

    def calculate_centralities(self, g, eigen_iter=1000):
        """
        Calculate harmonic, degree, closeness, eigen and betweenness centrality measures for an undirected network

        :param g: (G) Undirected networkx graph
        :param eigen_iter: (INT) Number of iterations for eigen centrality
        :return centralities: (DF) Pandas DF of centrality measures for each gene in the network
        """

        degree = self.degree_centrality(g)
        close = self.closeness(g)
        harm = self.harmonic_centrality(g)
        eigen = self.eigen_centrality(g, eigen_iter)
        between = self.betweenness_centrality(g)
        centralities = pd.concat([pd.Series(c) for c in (degree, close, harm, eigen, between)], axis=1)
        centralities.columns = ("Degree", "Closeness", "Harmonic", "Eigen", "Betweenness")
        return centralities

    def goi_centralities(self, gene, g, eigen_iter=1000):
        """
        Calculate harmonic, degree, closeness, eigen and betweenness centrality measures for a gene of interest
        in an undirected network

        :param gene: (STR) Name of gene of interest
        :param g: (G) Undirected networkx graph
        :param eigen_iter: (INT) Number of iterations for eigen centrality
        :return centralities: (DF) Pandas DF of centrality measures for the gene of interest in the network
        """

        deg = self.degree_centrality(g)
        degree = deg[gene]
        clo = self.closeness(g)
        close = clo[gene]
        ha = self.harmonic_centrality(g)
        harm = ha[gene]
        ei = self.eigen_centrality(g, eigen_iter)
        eigen = ei[gene]
        bet = self.betweenness_centrality(g)
        between = bet[gene]
        density = self.ntwk_desnity(g)
        centralities = pd.concat([pd.Series(c) for c in (degree, close, harm, eigen, between, density)], axis=1)
        centralities.columns = ("Degree", "Closeness", "Harmonic", "Eigen", "Betweenness", 'Density')
        centralities = centralities.reset_index()
        centralities.rename(columns={'index': 'Gene'}, inplace=True)
        centralities.loc[0, 'Gene'] = gene
        return centralities

    def gene_centralities(self, dict_genes):
        """
        Compute centrality measures for a gene of interest for a body of networks

        :param dict_genes: (DICT) Dictionary of dictionaries with genes as keys and networkx graphs stored for each gene
                           under the key 'Network'
        :return dict_results: (DICT) Dictionary with genes as keys and centrality DF for gene of interest as values
        """

        dict_results = {}
        for gene in dict_genes:
            g = dict_genes[gene]["Network"]
            centralities = self.calculate_centralities(g)
            dict_results[gene] = centralities
        return dict_results

#%%
# RANDOM NETWORK GENERATION AND ANALYSIS FOR STATISTICAL TESTING


class RandomNetworks:
    """
    Class to generate random networks and perform single gene correlation analysis on the networks
    """
    def __init__(self, df, gene_dict):
        """
        Constructor class for RandomNetworks

        :param df: (DF) Achilles clean CRISPR dataset as pandas DF
        :param gene_dict: (DICT) Single gene analysis results for well characterised genes as dictionary of
                           dictionaries with genes
        """

        self.df = df
        self.gene_dict = gene_dict
        self.gc = GeneCorrelations()
        self.oa = OntologyAnalysis()
        self.ntwk = NetworkAnalysis()

    def random_genelist(self, gene_dict, number_genes=100, ran_seed=16):
        """
        Generate lists of genes randomly sampled from the Achilles dataset. Generate the same number of lists as the
        number of well characterised genes analysed in 'gene_dict'

        :param gene_dict: (DICT) Single gene analysis results for well characterised genes as dictionary of
                           dictionaries with genes as keys
        :param number_genes: (INT) Number of genes to randomly sample, default=100
        :param ran_seed: (INT) Random seed for reproducibility, default=16
        :return random_lists: (LIST of LIST) List of lists of randomly sampled genes
        """

        genes = list(self.df.columns)
        number = len(gene_dict.keys())
        random.seed(ran_seed)
        random_lists = []
        for i in range(0, number, 1):
            # Make list of random genes
            rand = random.sample(genes,
                                 number_genes
                                 )
            random_lists.append(rand)
            ran_seed += 1
            random.seed(ran_seed)
        return random_lists

    def generate_random_network(self, gene_list):
        """
        For a specified gene list, generate the STRING based network

        :param gene_list: (LIST) List of genes
        :return g: (G) Netowrkx undirected graph
        :return  interactions: (DF) Dataframe with network data + weights as STRING combined score
        """

        interactions = self.ntwk.generate_interactions(gene_list)
        g = self.ntwk.network_import(interactions)
        g.remove_edges_from(nx.selfloop_edges(g))
        return g, interactions

    def random_goi_analysis(self, gene_list):
        """
        For a lit of randomly sampled genes, perform a single gene analysis to generate correlations, network
        and interactions
        :param gene_list: (LIST) List of randomly smapled genes
        :return result_dict: (DICT) Dictionary of results with the following keys + values:
                              Correlation: (DF) Correlation DF
                              Interaction: (DF) Dataframe with network data + weights as STRING combined score
                              Network: (G) Networkx undirected graph
        """

        result_dict = {}
        g, interactions = self.generate_random_network(gene_list)
        time.sleep(1)
        result_dict['Correlation'] = gene_list
        result_dict['Interaction'] = interactions
        result_dict['Network'] = g
        return result_dict

    def random_network_analysis(self):
        """
        Generate and analyse a body of random networks the same size as the single gene analysis results dictionary

        :return random_dict: (DICT) Dictionary of dictionaries with each gene as a key and the gene's dictionary
                              with the following keys + values:
                              Correlation: (DF) Correlation DF
                              Interaction: (DF) Dataframe with network data + weights as STRING combined score
                              Network: (G) Networkx undirected graph
        """

        random_dict = {}
        random_lists = self.random_genelist(self.gene_dict)
        for lst in random_lists:
            gene = lst[0]
            result_dict = self.random_goi_analysis(lst)
            random_dict[gene] = result_dict
        return random_dict

#%%
# NETWORK PERMUTATION AND EDGE DENSITY FOR COMMUNITY STATISTICAL ANALYSIS


class Permutations:
    """
    Class to randomly permute networks and calcularte an empirical P value. Methodology implemented as described
    by by Pan et al. (2018): https://pubmed.ncbi.nlm.nih.gov/29778836/
    """
    def __init__(self):
        """
        Constructor for CentralityAnalysis Class
        """

    @staticmethod
    def rewire(g, ran_seed=16):
        """
        Randomly rewire a network maintaining the same node degree
        :param g: (G) Networkx functional interaction undirected graph
        :param ran_seed: (INT) Random seed for random rewire, default=16
        :return: (G) Randomly rewired undirected graph with node degree maintained
        """

        network = g.copy()
        rewired = nx.generators.expected_degree_graph([deg for (name, deg) in network.degree()],
                                                      selfloops=False,
                                                      seed=ran_seed
                                                      )
        names = [node for (node, val) in network.degree]
        values = [node for (node, val) in rewired.degree]
        mapping = dict(zip(values,
                           names)
                       )
        rewired = nx.relabel_nodes(rewired,
                                   mapping,
                                   copy=False
                                   )
        return rewired

    @staticmethod
    def edge_den_ratio(g, community):
        """
        Compute the internal:external edge density ratio.
        :param g: (G) Networkx rewired undirected graph
        :param community: (LIST) Genes in the community being tested for significance
        :return: (FLOAT) Internal:external edge density ratio
        """

        # Compute internal edges and total possible internal edges of the community
        nodes_set = set(community)
        internal_e = len([(a, b) for a, b in g.edges(nodes_set) if a in nodes_set and b in nodes_set])
        internal_poss = comb(len(community), 2)
        internal_density = internal_e / internal_poss
        # Compute external edges and total possible external edges of the community
        number_edges = []
        for node in community:
            edges = g.edges(node)
            number_edges.append(len(edges))
        external_e = sum(number_edges) - internal_e
        external_poss = (len(g.nodes()) - len(community)) * len(community)
        external_density = external_e / external_poss
        try:
            return internal_density / external_density
        except ZeroDivisionError:
            return 'Error'

    def permute_network(self, g, community, permutations=10000, seed=16):
        """
        For a set number of permutations, rewire the network of interest and compute the community internal:external
        edge density ratio. Rewire network maintaining node degree for fair analysis
        :param g: (G) Networkx rewired undirected graph
        :param community: (LIST) Genes in the community being tested for significance
        :param permutations: (INT) Number of times to randomly rewire the network, default=10,000
        :param seed: (INT) Random seed for reproducible results, default=16
        :return: (LIST) List of edge density ratios from network permutations
        """

        ratios = []
        seed_start = seed
        for i in range(0, permutations, 1):
            ran_net = self.rewire(g, ran_seed=seed_start)
            ratio = self.edge_den_ratio(ran_net, community)
            if ratio == 'Error':
                continue
            else:
                ratios.append(ratio)
                seed_start += 1
        return ratios

    def empirical_p(self, g, community, ratios):
        """
        Calculate the empirical p value from the edge density ratio values from the permutation analysis
        :param g: (G) Networkx rewired undirected graph
        :param community: (LIST) Genes in the community being tested for significance
        :param ratios: (LIST of FLOAT) List of edge density ratio values for the permuted networks
        :return: (FLOAT) Empirical P value computed with bias
        """

        g_ratio = self.edge_den_ratio(g, community)
        n_above = [x for x in ratios if x >= g_ratio]
        return (len(n_above)+1) / (len(ratios)+1)

    def permutation_analysis(self,
                             gene_dict,
                             com_dict,
                             num_permute=10000,
                             ran_seed=16,
                             fdr_alpha=0.05,
                             fdr_method='indep'):
        """
        Permute each network in the gene_dict and compute the empirical p-value based on the edge density ratio.
        Adjust the final p-values for false discovery rate using the statsmodels multitest fdr correction. Default
        methodology set to 'indep' Benjamini-Hochberg
        :param gene_dict: (DICT) Single gene analysis results for well characterised genes as dictionary of
                           dictionaries with genes as the keys
        :param com_dict: (DICT) Community detection results for well characterised genes as dictionary of
                           dictionaries with genes as the keys
        :param num_permute: (INT) Number of network permutations, default=10,000
        :param ran_seed: (INT) Random seed number, default=16
        :param fdr_alpha: (FLOAT) False discovery rate cut-off
        :param fdr_method: (STR) False discovery rate methodology, default='indep'
        :return: (DICT) P-value and adjusted p-value (values) of each gene (keys)
        """

        # Perform permutation analysis to get empirical p-values
        genes = []
        p_vals = []
        for key in gene_dict:
            if len(com_dict[key]['Community']) >= 6:
                g = gene_dict[key]['Network']
                community = list(com_dict[key]['Community']['gene'])
                ratios = self.permute_network(g,
                                              community,
                                              permutations=num_permute,
                                              seed=ran_seed
                                              )
                p = self.empirical_p(g,
                                     community,
                                     ratios
                                     )
                p_vals.append(p)
                genes.append(key)
            else:
                continue
        # Adjust p-values for FDR
        rej, adj_p = multitest.fdrcorrection(p_vals,
                                             alpha=fdr_alpha,
                                             method=fdr_method
                                             )
        # Create final dictionary of results with genes (keys) and p-value + adjusted p-value (values)
        final_dict = dict(zip(genes, zip(p_vals,
                                         adj_p)
                              )
                          )
        return final_dict

    @staticmethod
    def make_sig_dict(com_dict, emp_dict):
        """
        Filter the community analysis dictionary for only significant communities identified by the network
        permutation analysis.
        :param com_dict: (DICT) Community detection results for well characterised genes as dictionary of
                           dictionaries with genes as the keys
        :param emp_dict: (DICT) Dictionary of network permutation results. P-value and adjusted p-value (values) of
                          each gene (keys)
        :return: (DICT) Final dictionary filtered for only significant communities
        """

        df = pd.DataFrame.from_dict(emp_dict, orient='index', columns=['p_val', 'adj_p'])
        df = df.reset_index(level=0)
        df = df.rename(columns={'index': 'gene_name', 'p_val': 'p_val', 'adj_p': 'adj_p'})
        # Get list of genes with significant communities
        df_sig = df[df['adj_p'] <= 0.01]
        sig_genes = df_sig['gene_name'].tolist()
        final_k = {}
        for gene in sig_genes:
            temp = {'Community': com_dict[gene]['Community'], 'GSEA': com_dict[gene]['GSEA']}
            final_k[gene] = temp
        return final_k
