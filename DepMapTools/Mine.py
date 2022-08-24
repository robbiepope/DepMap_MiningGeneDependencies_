# Import packages
import pandas as pd
from Bio import Entrez
# Ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
#%%
# STRING BASED NETWORK FOR CORRELATED GENES


class MineData:
    """
    Class to mine the DepMap dataset as well as the scientific literature
    """

    @staticmethod
    def mine_data(term, gsea_dict, p_val=0.01, enrich_score=4):
        """
        Search community GSEA Gene Ontologies for GOs containing a specific term and return list of genes
        enriched for that term

        :param term: (STR) Term to search GOs
        :param gsea_dict: (DICT) Dictionary of GSEA for communities (values) for genes of interest (keys)
        :param p_val: (FLOAT) P-value threshold, only search GO terms that are below specified p-value, default=0.01
        :param enrich_score: (INT) Number of GO terms for that community that contain the search term, default=4
        :return gene_list: (LIST of STR) List of genes that have communities enriched in the search term
        """

        gene_dict = {}
        # For each gene's GSEA, if valid, extract the GO terms
        for gene in gsea_dict:
            gsea = gsea_dict[gene]['GSEA']
            if isinstance(gsea, list):
                continue
            else:
                gsea = gsea.sort_values(by='Adjusted P-value')
                gsea = gsea[gsea['Adjusted P-value'] <= p_val]
                gsea_vals = gsea['Term'].to_list()
                score = 0
            # For each extracted GO term, if the search term is present add 1 to the score
            for val in gsea_vals:
                if term in val:
                    score += 1
                else:
                    continue
            # If the score is greater or equal to the enrichment score, append gene name to gene_list
            if score >= enrich_score:
                gene_dict[gene] = score
        return gene_dict

    def search_terms(self, term_dict, gsea_dict, p_val=0.01):
        """
        Mine the larger DepMap dataset for a selection of terms and return a list of genes with enriched communities
        for each term

        :param term_dict: (DICT) Dictionary of terms (keys) with enrichment scores (values)
        :param gsea_dict: (DICT) GSEA results for a body of genes as a dictionary of dictionary with each gene as a key
                          and subsequent keys ['Community'] and ['GSEA']
        :param p_val: (FLOAT) P-value threshold, only search GO terms that are below specified p-value, default=0.01
        :return gene_dict: (DICT) Dictionary of terms (keys) with list of genes that have communities enriched in
                            the search term (values)
        """

        final_dict = {}
        # For each term in term dictionary, mine dataset. Take term_dict value as enrichment score
        for term in term_dict:
            gene_dict = self.mine_data(term,
                                       gsea_dict,
                                       p_val=p_val,
                                       enrich_score=term_dict[term]
                                       )
            gene_list = list(gene_dict.keys())
            values = list(gene_dict.values())
            final_dict[term] = {'genes': gene_list, 'score': values}
        return final_dict

    def mine_terms(self, term_dict, gsea_dict, p_val=0.01):
        """
        Mine the larger DepMap dataset for a selection of terms and return a finl DF of genes with enriched communities
        for each term

        :param term_dict: (DICT) Dictionary of terms (keys) with enrichment scores (values)
        :param gsea_dict: (DICT) GSEA results for a body of genes as a dictionary of dictionary with each gene as a key
                          and subsequent keys ['Community'] and ['GSEA']
        :param p_val: (FLOAT) P-value threshold, only search GO terms that are below specified p-value, default=0.01
        :return final_df: (DF) Pandas DF of genes and their enriched term
        """

        print('Searching Terms')
        print('-'*20)
        gene_dict = self.search_terms(term_dict, gsea_dict, p_val=p_val)
        for key in gene_dict:
            print(f'{key} : {len(gene_dict[key])} genes')
        print('-'*20)
        print('Compiling DF')
        final_df = pd.DataFrame()
        for key in gene_dict:
            gene_list = gene_dict[key]
            df = pd.DataFrame(gene_list).assign(Term=key)
            final_df = final_df.append(df, ignore_index=True)
        print('-'*20)
        print('Complete')
        return final_df

    @staticmethod
    def get_publications(gene_list, search_term, email, database='pubmed', max_retain='40'):
        """
        Using BioPython Entrez package mine a specified database for publications that contain a specified search term

        :param gene_list: (LIST of STR) List of genes to search in the database
        :param search_term: (STR) Additional term to search the literature for as well as the gene name
        :param email: (STR) User's email address so that NCBI can identify the caller
        :param database: (STR) Database to search, default='pubmed'
        :param max_retain: (STR) Number of records to retain in memory, default='40'
        :return df_results: (DF) Indexed by gene with count of publications containing gene name and gene name + term
        """

        Entrez.email = email
        final_results = {}
        # For every gene search the specified database for gene and gene + search_term
        for gene in gene_list:
            dict_results = {}
            search = gene+" AND "+search_term
            handle = Entrez.esearch(db=database,
                                    term=gene,
                                    retmax=max_retain
                                    )
            record = Entrez.read(handle)
            dict_results['gene'] = int(record['Count'])
            handle1 = Entrez.esearch(db=database,
                                     term=search,
                                     retmax=max_retain
                                     )
            record1 = Entrez.read(handle1)
            dict_results['gene_term'] = int(record1['Count'])
            final_results[gene] = dict_results
        # Return final results as a pandas DF
        df_results = pd.DataFrame.from_dict(final_results,
                                            orient='index'
                                            )
        df_results = df_results.sort_values('gene')
        return df_results
