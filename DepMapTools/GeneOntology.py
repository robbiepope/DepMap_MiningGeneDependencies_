# Import packages
import os
import wget
import gseapy as gp
from goatools.mapslim import mapslim
from goatools.obo_parser import GODag
from goatools.semsim.termwise.wang import SsWang
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.base import download_ncbi_associations, get_godag
from goatools.cli.ncbi_gene_results_to_python import ncbi_tsv_to_py
# Ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
#%%
# GENE ONTOLOGY ANALYSIS


class OntologyAnalysis:
    """
    Class to perform gene ontology analysis for networks
    """

    @staticmethod
    def download_ncbi(filename='gene_result.txt', output='genes_ncbi_9606_proteincoding.py'):
        """
        From the relevant NCBI protein coding txt file create a .py file containing the GENEID2NT dictionary

        :param filename: (STR) Name of NCBI file downloaded, default='gene_result.txt'
        :param output:  (STR) Name of output .py file, default='genes_ncbi_9606_proteincoding.py'
        :return: (PY) NCBI txt file converted to .py saved in the same directory
        """

        ncbi_tsv = filename
        output_py = output
        return ncbi_tsv_to_py(ncbi_tsv, output_py)

    @staticmethod
    def gene_mapper(ncbi_dict):
        """
        Map the name of a gene to its relevant NCBI gene ID

        :param ncbi_dict: (DICT) GENE2ID NCBI dictionary of gene names and IDs
        :return mapper: (DICT) Dictionary of gene names (keys) and gene IDs (values)
        """

        mapper = {}
        for key in ncbi_dict:
            mapper[ncbi_dict[key].Symbol] = ncbi_dict[key].GeneID
        return mapper

    @staticmethod
    def get_ncbi_associations(species=9606):
        """
        Download the NCBI association files to retrieve gene ontologies
        :param species: (INT) TaxId for desired species, default=9606 (human)
        :return ns2assc: (DICT) Given associations returned in three dictionaries (BP, MF, CC)
        """

        gene2go = download_ncbi_associations()
        annotation = Gene2GoReader(gene2go,
                                   taxids=[species]
                                   )
        ns2assc = annotation.get_ns2assc()
        return ns2assc

    @staticmethod
    def gsea(genelist, database=['GO_Biological_Process_2021']):
        """
        Perform a GSEA analysis for a specified gene list by querying the EnrichR database

        :param genelist: (LIST of STR) List of gene names
        :param database: (LIST of STR) List of database names to query, default='GO_Biological_Process_2021'
        :return df_go: (DF) Pandas DF of GSEA results
        """

        enr = gp.enrichr(gene_list=genelist,
                         gene_sets=database,
                         organism='Human',
                         cutoff=0.5)
        df_go = enr.results
        return df_go

    @staticmethod
    def get_obo(url, description):
        """
        Download specified OBO file and save in directory

        :param url: (STR) URL of OBO file to download
        :param description: (STR) Description of OBO file (basic or slim)
        :return go_obo: (OBO) OBO file saved in specified directory
        """

        data_folder = os.getcwd() + '/GOdata3'
        if not os.path.isfile(data_folder):
            try:
                os.mkdir(data_folder)
            except OSError as e:
                if e.errno != 17:
                    raise e
        else:
            raise Exception('Data path exists as a file')
        if not os.path.isfile(data_folder + f'/go-{description}.obo'):
            go_obo = wget.download(url, data_folder + f'/go-{description}.obo')
        else:
            go_obo = data_folder + f'/go-{description}.obo'
        return go_obo

    @staticmethod
    def get_goids(go_list):
        """
        From a list of gene ontologies extract the numerical GO ID for each GO

        :param go_list: (LIST of STR) List of GOs containing numerical GO IDs
        :return goids: (LIST of STR) List of isolated GO IDs in format 'GO: XXXXXXXX'
        """

        goids = []
        for go in go_list:
            s = go
            val = s.split('GO', 1)[1].split(')')[0]
            val = "GO"+val
            goids.append(val)
        return goids

    @staticmethod
    def map_slim_go(go_obo, go_slim, goids):
        """
        Map complex GOs to the relevant GO slim term
        :param go_obo: (OBO) basic OBO file as python object
        :param go_slim: (OBO) slim OBO file as python object
        :param goids: (LIST) List of GO IDs in format 'GO: XXXXXXXX'
        :return dict_slim: (DICT) Dictionary of slim GOs with GO IDs (keys) and descriptions (values)
        """

        go = GODag(go_obo)
        goslimdag = GODag(go_slim)
        direct_anc = []
        all_anc = []
        # Map basic GO ID to slim GO ID
        for gid in goids:
            try:
                anc, allanc = mapslim(gid, go, goslimdag)
                for term in anc:
                    if term not in direct_anc:
                        term2 = term
                        direct_anc.append(term2)
                for term in allanc:
                    if term not in all_anc:
                        term2 = term
                        all_anc.append(term2)
            except:
                pass
        # Package into dictionary of keys (GO Term) and values (Description)
        values = []
        for onto in direct_anc:
            go_term = go[onto]
            values.append(go_term.name)
        dict_slim = dict(zip(direct_anc, values))
        return dict_slim

    @staticmethod
    def set_godag(obo_file='go-basic.obo'):
        """
        Get the GO directed acyclic graph (DAG) for specified OBO file for similarity analysis. Optional attributes
        set to 'relationship'

        :param obo_file: (STR) Name of OBO file to use to create DAG
        :return godag: (DAG) DAG of specified OBO file with optional attribute set to 'relationship'
        """

        godag = get_godag(obo_file, optional_attrs={'relationship'})
        return godag

    @staticmethod
    def get_list_goids(goi_slim, com_slim):
        """
        Generate a list of total unique slim GO IDs for the community and gene of interest

        :param goi_slim: (DICT) Dictionary of slim GO terms for gene of interest (keys=GO IDs)
        :param com_slim: (DICT) Dictionary of slim GO terms for community (keys=GO IDs)
        :return goids: (LIST of STR) List of unique slim GO IDs for the community and gene of interest
        """

        goids = list(com_slim.keys())
        goi_list = list(goi_slim.keys())
        for gid in goi_list:
            if gid not in goids:
                goids.append(gid)
        return goids

    @staticmethod
    def wang_sim(goids, godag, goi_slim, com_slim, relationships={'part_of'}):
        """
        Compute Wangs Semantic Similarity Score for slim GOs for the gene of interest and community gene set

        :param goids: (LIST) List of total slim GO IDs present for gene of interest and community
        :param godag: (DAG) DAG of specified OBO file with optional attribute set to 'relationship'
        :param goi_slim: (LIST) List of slim GO IDs present for gene of interest and community
        :param com_slim: (LIST) List of slim GO IDs present for community
        :param relationships: (DICT of STR) Define GO relationship for analysis, default='part_of'
        :return score: (FLOAT) Wang's Semantic Similarity (0-1) for community and GOI
        :return error: (STR) 'Error' for similarity scores divisible by 0
        """

        wang = SsWang(goids, godag, relationships)
        # Perform pairwise comparison of GOs WANG (SGO)
        # Extract highest similarity for community and GOI
        sims_goi = {}
        sims_com = {}
        for go_a in goi_slim:
            value_a = 0
            for go_b in com_slim:
                val = wang.get_sim(go_a, go_b)
                if val > value_a:
                    value_a = val
                    sims_goi[go_a] = val
        for go_b in com_slim:
            value_b = 0
            for go_a in goi_slim:
                val = wang.get_sim(go_b, go_a)
                if val > value_b:
                    value_b = val
                    sims_com[go_b] = val
        # Calculate wang similarity between GOI and community
        goi_list = list(sims_goi.values())
        com_list = list(sims_com.values())
        try:
            return (sum(goi_list) + sum(com_list)) / (len(goi_list) + len(com_list))
        except ZeroDivisionError:
            return "Error"

    def make_go_dict(self, ncbi_dict, gsea_dict, sortby="Adjusted P-value", threshold=0.01):
        """
        Create a dictionary of dictionaries of a body of genes with goi_go and com_go as keys containing the GO terms
        as values for the GOI and community specifically looking at "Biological Process" GO terms

        :param ncbi_dict: (DICT) NCBI associations dictionary
        :param gsea_dict: (DICT) GSEA results for a body of genes as a dictionary of dictionary with each gene as a key
                          and subsequent keys ['Community'] and ['GSEA']
        :param sortby: (STR) GSEA column to sort results by, default=Adjusted P-Value' with most significant first
        :param threshold: (FLOAT) Only consider GO terms below the Adjusted P-value threshold (most significant terms)
        :return go_dict: (DICT) Dictionry of dictionaries with genes as the keys. Each key then has the keys ['goi_go']
                          and ['com_gsea] with list GO terms as values
        :return small_com: (LIST of STR) List of genes with communities too small to analyse
        :return missing_genes: (LIST of STR) List of genes not present in the NCBI associations dict
        """

        mapper = self.gene_mapper(ncbi_dict)
        associations = self.get_ncbi_associations()
        list_ids = list(associations["BP"].keys())
        go_dict = {}
        small_com = []
        missing_genes = []
        # For each gene get the GO terms fot that gene and extract the GSEA of the ocmmunity
        for key in gsea_dict:
            results_dict = {}
            try:
                gene = mapper[key]
                # If gene is missing from the NCBI associations dict then skip
                if gene not in list_ids:
                    missing_genes.append(key)
                    continue
                else:
                    go_goi = associations['BP'][gene]
                    gsea_com = gsea_dict[key]["GSEA"]
                    # If community is too small then skip
                    if isinstance(gsea_com, list):
                        small_com.append(key)
                        continue
                    # Create results dict of GOI GO terms and community GSEA GO terms
                    else:
                        gsea_com = gsea_com.sort_values(by=sortby)
                        gsea_com = gsea_com[gsea_com['Adjusted P-value'] < threshold]
                        results_dict['goi_go'] = list(go_goi)
                        results_dict['com_gsea'] = list(gsea_com["Term"])
                        go_dict[key] = results_dict
            except Exception as e:
                if str(e) == key:
                    missing_genes.append(key)
                    print(f'{key} error')
                    continue
        return go_dict, small_com, missing_genes

    def wang_sim_analysis(self, go_dict, url_basic, url_slim):
        """
        Compute Wang's Semantic Similarity Score for a dictionary of dictionaries of genes (keys) with GOI GO terms
        and community GO terms as the values.

        :param go_dict: (DICT) Dictionary of dictionaries with genes of interest as the keys then for each gene the GO
                        terms for the GOI and community stored under the keys ['goi_go'] and ['com_gsea']
        :param url_basic: (STR) URL of basic OBO file to download
        :param url_slim: (STR) URL of slim OBO file to download
        :return dict_results: (DICT) Dictionary of results with genes (keys) and wang similarity score (values)
        :return no_go: (LIST of STR) List of genes with no GO enrichment
        """

        dict_results = {}
        godag = self.set_godag()
        go_obo = self.get_obo(url_basic, 'basic')
        go_slim = self.get_obo(url_slim, 'slim')
        no_go = []
        # For each gene collect the GO IDs and map the terms to slim terms
        for key in go_dict:
            print(key)
            goi_ids = self.get_goids(go_dict[key]["goi_go"])
            com_ids = self.get_goids(go_dict[key]["com_gsea"])
            goi_slim = self.map_slim_go(go_obo, go_slim, goi_ids)
            com_slim = self.map_slim_go(go_obo, go_slim, com_ids)
            # If no GO enrichment for community and GOI, then skip gene, otherwise compute similarity
            if len(goi_slim.keys()) and len(com_slim.keys()) == 0:
                no_go.append(key)
                continue
            else:
                goids = self.get_list_goids(goi_slim, com_slim)
                final_sim = self.wang_sim(goids, godag, goi_slim, com_slim)
            # If error raised, append gene to no_go list and skip, otherwise add gene to dictionary of results
            if final_sim == "Error":
                no_go.append(key)
                continue
            else:
                dict_results[key] = final_sim
        return dict_results, no_go

    @staticmethod
    def get_gsea_terms(gene_dict):
        """
        From dictionary of dictionaries with genes as keys get community GSEA results and extract GO terms

        :param gene_dict: (DICT) Dictionary of dictionaries with genes as keys and GSEA analysis stored as value under
                           the key ['GSEA']
        :return go_dict: (DICT) Dictionary with genes as keys and values as a list of the gene's GO terms
        :return no_go_enrichment: (LIST of STR) List of genes with no GO enrichment
        """

        go_dict = {}
        no_go_enrichment = []
        for gene in gene_dict:
            gsea = gene_dict[gene]['GSEA']
            if isinstance(gsea, list):
                no_go_enrichment.append(gene)
                continue
            else:
                gsea = gsea.sort_values(by='Adjusted P-value')
                gsea = gsea[gsea['Adjusted P-value'] < 0.01]
                go_list = gsea['Term'].to_list()
                go_dict[gene] = go_list
        return go_dict, no_go_enrichment

    def get_gsea_slim(self, go_dict, url_basic, url_slim):
        """
        From dictionary of genes (keys) and community GSEA GO terms, map to slim GO terms

        :param go_dict: (DICT) Dictionary of genes (keys) and a list of community GSEA GO terms (values) for each gene
        :param url_basic: (STR) URL of basic OBO file to download
        :param url_slim: (STR) URL of slim OBO file to download
        :return dict_results: (DICT) Dictionary of results with genes (keys) and list of slim GO terms (values)
        """

        dict_results = {}
        go_obo = self.get_obo(url_basic, 'basic')
        go_slim = self.get_obo(url_slim, 'slim')
        for key in go_dict:
            goids = self.get_goids(go_dict[key])
            goids_slim = self.map_slim_go(go_obo, go_slim, goids)
            if len(goids_slim.keys()) == 0:
                final_ids = go_dict[key]
                dict_results[key] = final_ids
            else:
                dict_results[key] = goids_slim
        return dict_results
