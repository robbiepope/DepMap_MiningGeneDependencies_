# Import packages
import time
import pandas as pd
import networkx as nx
from DepMapTools.GeneOntology import OntologyAnalysis
from DepMapTools.Networks import NetworkAnalysis
# Ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
#%%
# COMMUNITY ANALYSIS


class CommunityAnalysis:
    """
    Class to detect communities via various algorithms and to perform community analysis
    """

    def __init__(self):
        """
        Constructor for CommunityAnalysis class
        """
        self.oa = OntologyAnalysis()
        self.nn = NetworkAnalysis()

    @staticmethod
    def k_clique_community(g, k=3):
        """
        K-clique percolation community detection algorithm. Detect communities if a clique >= k exists in the network
        and return a list of gene in the largest community. If no clique exists, return 'Community too Small'

        :param g: (G) Networkx undirected graph
        :param k: (INT) Number of minimum nodes that constitutes a cliqe, default=3
        :return genelist: (LIST of STR) List of genes in largest network community
        """

        # Only perform k-clique community detection if the a clique >= k exists
        max_clique = nx.cliques = nx.graph_clique_number(g)
        if max_clique >= k:
            community = sorted(nx.algorithms.community.k_clique_communities(g, k),
                               key=lambda x: len(x)
                               )
        else:
            community = []
        # If the length of the community is greater than 0 return community genes as list
        if len(list(community)) > 0:
            genelist = list(list(community)[-1])
        else:
            genelist = ['Community too Small']
        return genelist

    def k_clique_com_analysis(self, gene_dict):
        """
        For a selection of genes of interest detect the k-clique communities and perform a GSEA for any viable
        communities
        Due to limitations of the Enrichr server, error mitigation implemented via a while loop to ensure a GSEA can be
        performed for all community gene lists

        :param gene_dict: (DICT) Dictionary of dictionaries with genes as the keys and the following keys for each gene:
                           Correlation: (DF) Correlation DF
                           Interaction: (DF) Dataframe with network data + weights as STRING combined score
                           Network: (G) Networkx undirected graph

        :return k_dict: (DICT) Dictionary of dictionaries with genes as the keys and the following keys for each gene:
                         Community: (LIST of STR) Names of genes in the largest community of the network
                         GSEA: (DF) Results of the GSEA as a pandas DF
        """

        list_of_genes = list(gene_dict.keys())
        k_dict = {}
        # While length of gene list is greater than 0 call Enrichr GSEA
        while len(list_of_genes) > 0:
            for key in list_of_genes:
                # For each gene try to generate the community and GSEA
                try:
                    result_dict = {}
                    g = gene_dict[key]["Network"]
                    com = self.k_clique_community(g)
                    if len(com) >= 6:
                        df_go = self.oa.gsea(com)
                    # If community length too small, return 'No GO Enrichment'
                    else:
                        df_go = ['No GO Enrichment']
                    # Create dictionary of results for community and GSEA
                    result_dict['Community'] = pd.DataFrame(com,
                                                            columns=['gene']
                                                            )
                    result_dict['GSEA'] = df_go
                    # Add dictionary of results to final dictionary of dictionaries
                    k_dict[key] = result_dict
                    # Remove the key (gene) from the gene list
                    list_of_genes.remove(key)
                    time.sleep(1)
                # If the 'Error analyzing gene list' is raised, skip the gene and start again
                except Exception as e:
                    if str(e) == 'Error analyzing gene list':
                        continue
        else:
            print(f'{len(k_dict.keys())} Networks Analysed')
            print('Complete')
            print('-'*20)
            return k_dict

    @staticmethod
    def edge_to_remove(g):
        """
        Calculate edge betweenness centrality score for network and extract the edge with the highest score

        :param g: (G) Networkx undirected graph
        :return edge: (INT) Edge with the highest betweenness centrality score to be removed
        """
        ntwk_dict = nx.edge_betweenness_centrality(g)
        edge = ()
        for key, value in sorted(ntwk_dict.items(),
                                 key=lambda item: item[1],
                                 reverse=True
                                 ):
            edge = key
            break
        return edge

    def girvan_newman(self, g):
        """
        Iteratively remove edges based on betweenness score until the graph is deivided into >1 sub graphs
        :param g: (G) Networkx undirected graph
        :return sg: (SET) Frozen sets of nodes of each sub graph, representing the detected communities
        """

        # Find the number of connected sub graphs
        sg = nx.connected_components(g)
        sg_count = nx.number_connected_components(g)
        while sg_count == 1:
            g.remove_edge(self.edge_to_remove(g)[0],
                          self.edge_to_remove(g)[1]
                          )
            sg = nx.connected_components(g)
            sg_count = nx.number_connected_components(g)
        return sg

    def girvan_community(self, g):
        """
        Girvan Newman community detection algorithm. Detect communities based on betweenness score and sub graph
        connectivity and return a list of gene in the largest community. If no community exists, return
        'Community too Small'

        :param g: (G) Networkx undirected graph
        :return genelist: (LIST of STR) List of genes in largest network community
        """

        gc = g.copy()
        # Remove isolated nodes to generate a connected subgraph
        to_be_removed = [x for x in gc.nodes() if gc.degree(x) < 1]
        for x in to_be_removed:
            gc.remove_node(x)
        community = sorted(self.girvan_newman(gc),
                           key=lambda x: len(x)
                           )
        largest_com = list(max(community,
                               key=len)
                           )
        # If the length of the community is greater than 0 return community genes as list
        if len(largest_com) > 0:
            genelist = largest_com
        else:
            genelist = ['Community too Small']
        return genelist

    def girvan_com_analysis(self, gene_dict):
        """
        For a selection of genes of interest using the Girvan Newman community detection algorithm, find communities
        and perform a GSEA for any viable communities
        Due to limitations of the Enrichr server, error mitigation implemented via a while loop to ensure a GSEA can be
        performed for all community gene lists

        :param gene_dict: (DICT) Dictionary of dictionaries with genes as the keys and the following keys for each gene:
                           Correlation: (DF) Correlation DF
                           Interaction: (DF) Dataframe with network data + weights as STRING combined score
                           Network: (G) Networkx undirected graph

        :return girv_dict: (DICT) Dictionary of dictionaries with genes as keys and the following keys for each gene:
                           Community: (LIST of STR) Names of genes in the largest community of the network
                           GSEA: (DF) Results of the GSEA as a pandas DF
        """

        list_of_genes = list(gene_dict.keys())
        girv_dict = {}
        # While length of gene list is greater than 0 call Enrichr GSEA
        while len(list_of_genes) > 0:
            for key in list_of_genes:
                # For each gene try to generate the community and GSEA
                try:
                    result_dict = {}
                    g = gene_dict[key]["Network"]
                    com = self.girvan_community(g)
                    if len(com) >= 6:
                        df_go = self.oa.gsea(com)
                    # If community length too small, return 'No GO Enrichment'
                    else:
                        df_go = ['No GO Enrichment']
                    # Create dictionary of results for community and GSEA
                    result_dict['Community'] = pd.DataFrame(com,
                                                            columns=['gene']
                                                            )
                    result_dict['GSEA'] = df_go
                    # Add dictionary of results to final dictionary of dictionaries
                    girv_dict[key] = result_dict
                    # Remove the key (gene) from the gene list
                    list_of_genes.remove(key)
                    time.sleep(1)
                # If the 'Error analyzing gene list' is raised, skip the gene and start again
                except Exception as e:
                    if str(e) == 'Error analyzing gene list':
                        continue
        else:
            print(f'{len(girv_dict.keys())} Networks Analysed')
            print('Complete')
            print('-'*20)
            return girv_dict

    @staticmethod
    def louvain_community(g, resolution=1, threshold=0.0000001, seed=16):
        """
        Louvain modularity based community detection algorithm. Detect communities based on modularity gain and return
        a list of gene in the largest community. If no community exists, return 'Community too Small'

        :param g: (G) Networkx undirected graph
        :param resolution: (FLT) Size of communities favoured, <1=small, >1=large. default=1
        :param threshold: (FLT) Modularity gain threshold, below threshold and the algorithm stops, default=1e-7
        :param seed: (INT) Random seed for the algorithm to begin at, default=16
        :return genelist: (LIST of STR) List of genes in largest network community
        """
        community = nx.community.louvain_communities(g,
                                                     resolution=resolution,
                                                     threshold=threshold,
                                                     seed=seed
                                                     )
        largest_com = list(max(community,
                               key=len
                               )
                           )
        if len(largest_com) > 0:
            genelist = largest_com
        else:
            genelist = ['Community too 1_Small']
        return genelist

    def louvain_com_analysis(self, gene_dict, resolution=1, threshold=0.0000001, seed=16):
        """
        For a selection of genes of interest using the Louvain modularity community detection algorithm, find
        communities and perform a GSEA for any viable communities
        Due to limitations of the Enrichr server, error mitigation implemented via a while loop to ensure a GSEA can be
        performed for all community gene lists

        :param gene_dict: (DICT) Dictionary of dictionaries with genes as the keys and the following keys for each gene:
                           Correlation: (DF) Correlation DF
                           Interaction: (DF) Dataframe with network data + weights as STRING combined score
                           Network: (G) Networkx undirected graph
        :param resolution: (FLT) Size of communities favoured, <1=small, >1=large. default=1
        :param threshold: (FLT) Modularity gain threshold, below threshold and the algorithm stops, default=1e-7
        :param seed: (INT) Random seed for the algorithm to begin at, default=16
        :return louv_dict: (DICT) Dictionary of dictionaries with genes as keys and the following keys for each gene:
                           Community: (LIST of STR) Names of genes in the largest community of the network
                           GSEA: (DF) Results of the GSEA as a pandas DF
        """

        list_of_genes = list(gene_dict.keys())
        louv_dict = {}
        # While length of gene list is greater than 0 call Enrichr GSEA
        while len(list_of_genes) > 0:
            # For each gene try to generate the community and GSEA
            for key in list_of_genes:
                try:
                    result_dict = {}
                    g = gene_dict[key]["Network"]
                    com = self.louvain_community(g,
                                                 resolution=resolution,
                                                 threshold=threshold,
                                                 seed=seed
                                                 )
                    if len(com) >= 6:
                        df_go = self.oa.gsea(com)
                    # If community length too small, return 'No GO Enrichment'
                    else:
                        df_go = ['No GO Enrichment']
                    # Create dictionary of results for community and GSEA
                    result_dict['Community'] = pd.DataFrame(com,
                                                            columns=['gene']
                                                            )
                    result_dict['GSEA'] = df_go
                    # Add dictionary of results to final dictionary of dictionaries
                    louv_dict[key] = result_dict
                    # Remove the key (gene) from the gene list
                    list_of_genes.remove(key)
                    time.sleep(10)
                # If the 'Error analyzing gene list' is raised, skip the gene and start again
                except Exception as e:
                    if str(e) == 'Error analyzing gene list':
                        continue
        else:
            print(f'{len(louv_dict.keys())} Networks Analysed')
            print('Complete')
            print('-'*20)
            return louv_dict

    @staticmethod
    def labelprop_community(g):
        """
        Label propagation community detection algorithm. Semi supervised approach to detect communities and return
        a list of gene in the largest community. If no community exists, return 'Community too Small'

        :param g: (G) Networkx undirected graph
        :return genelist: (LIST of STR) List of genes in largest network community
        """
        community = nx.community.label_propagation_communities(g)
        largest_com = list(max(community,
                               key=len)
                           )
        if len(largest_com) > 0:
            genelist = largest_com
        else:
            genelist = ['Community too 1_Small']
        return genelist

    def label_com_analysis(self, gene_dict):
        """
        For a selection of genes of interest using the Label propagation community detection algorithm, find
        communities and perform a GSEA for any viable communities
        Due to limitations of the Enrichr server, error mitigation implemented via a while loop to ensure a GSEA can be
        performed for all community gene lists

        :param gene_dict: (DICT) Dictionary of dictionaries with genes as the keys and the following keys for each gene:
                           Correlation: (DF) Correlation DF
                           Interaction: (DF) Dataframe with network data + weights as STRING combined score
                           Network: (G) Networkx undirected graph
        :return label_dict: (DICT) Dictionary of dictionaries with genes as keys and the following keys for each gene:
                             Community: (LIST of STR) Names of genes in the largest community of the network
                             GSEA: (DF) Results of the GSEA as a pandas DF
        """

        list_of_genes = list(gene_dict.keys())
        label_dict = {}
        # While length of gene list is greater than 0 call Enrichr GSEA
        while len(list_of_genes) > 0:
            # For each gene try to generate the community and GSEA
            for key in list_of_genes:
                try:
                    result_dict = {}
                    g = gene_dict[key]["Network"]
                    com = self.labelprop_community(g)
                    if len(com) >= 6:
                        df_go = self.oa.gsea(com)
                    # If community length too small, return 'No GO Enrichment'
                    else:
                        df_go = ['No GO Enrichment']
                    # Create dictionary of results for community and GSEA
                    result_dict['Community'] = pd.DataFrame(com, columns=['gene'])
                    result_dict['GSEA'] = df_go
                    # Add dictionary of results to final dictionary of dictionaries
                    label_dict[key] = result_dict
                    # Remove the key (gene) from the gene list
                    list_of_genes.remove(key)
                    time.sleep(1)
                # If the 'Error analyzing gene list' is raised, skip the gene and start again
                except Exception as e:
                    if str(e) == 'Error analyzing gene list':
                        continue
        else:
            print(f'{len(label_dict.keys())} Networks Analysed')
            print('Complete')
            print('-'*20)
            return label_dict
