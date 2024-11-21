import requests
import urllib.parse
import json
import pandas as pd
from tqdm import tqdm
from requests.packages.urllib3.exceptions import InsecureRequestWarning

requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
#HPO_SEARCH_URL = 'https://hpo.jax.org/api/hpo/search'

class Hpo:
    def __init__(self, HPO_SEARCH_URL):
        self.HPO_SEARCH_URL = HPO_SEARCH_URL


    def create_overall_list(self, path_to_phenotypesList):
        """_summary_

            This module can create list of candidate genes.
            
            Args:
                path_to_phenotypesList (path): 

            Returns:
                dictionary : The return value. 
                            {'Inputed phenotype 1': 
                                {'Candidate phenotype 1': 
                                    {'Ontology ID': HP:******,
                                    'Genes': ['Gene 1', 'Gene 2', ...]
                                    },
                            'Inputed phenotype 2': ...
                            }
        """

        with open(path_to_phenotypesList, 'r', encoding='utf_8') as phenotypes_file:
            inputed_phenotypes = phenotypes_file.read().splitlines()
        
        phenotypes_to_genes = {} # Final output dictionary including inputed phenotypes and candidate ontology ID and genes

        for phenotype in tqdm(inputed_phenotypes):
            
            search_params = {
                'max': '-1',
                'offest': '0',
                'category': 'terms', 
                'q': phenotype
            }
        
            candidate_PhenotypesDictionary_list = []
            res_search_raw = requests.get(self.HPO_SEARCH_URL, verify=False, params=search_params)
            res_search_json = res_search_raw.json()
            totalNamesCount = len(res_search_json['terms'])
            phenotypeNames = [res_search_json['terms'][i]['name'] for i in range(totalNamesCount)]
                            
            for j in tqdm(range(totalNamesCount), leave=False):
                id = res_search_json['terms'][j]['ontologyId']
                ontologyId = urllib.parse.quote(id)
                HPO_TERM_URL = f'https://hpo.jax.org/api/hpo/term/{ontologyId}/genes'

                term_params = {
                    'max': '-1',
                    'offest': '0'
                }

                res_term_raw = requests.get(HPO_TERM_URL, verify=False, params = term_params)
                res_term_json = res_term_raw.json()
                totalGeneSymbolCount = len(res_term_json['genes'])
                genes_list = [res_term_json['genes'][k]['geneSymbol'] for k in range(totalGeneSymbolCount)]

                candidatePhenotyps_dic = dict(genes = genes_list, ontology_ID = id)
                candidate_PhenotypesDictionary_list.append(candidatePhenotyps_dic)
                
            matched_phenotypes_dic = dict(zip(phenotypeNames, candidate_PhenotypesDictionary_list))
            phenotypes_to_genes[phenotype] = matched_phenotypes_dic

        return phenotypes_to_genes

    def to_tsv(self, path_to_phenotypesList, GeneList_json, outputFilePath):
        with open(path_to_phenotypesList, 'r', encoding='utf_8') as phenotypes_file:
            inputed_phenotypes = phenotypes_file.read().splitlines()
        
        interval_lists = []

        for inputed_phenotype in inputed_phenotypes:
            for phenotype in GeneList_json[inputed_phenotype]:
                for gene in GeneList_json[inputed_phenotype][phenotype]['genes']:
                    id = GeneList_json[inputed_phenotype][phenotype]['ontology_ID']
                    interval_list = [inputed_phenotype, phenotype, id, gene]
                    interval_lists.append(interval_list)
        
        if len(interval_lists) != 0:
            df_without_columns = pd.concat([pd.Series(x) for x in interval_lists], axis=1).T
            df_without_columns.to_csv(outputFilePath, header=False, index=False, sep='\t')
            pass
        else:
            print('No HPO candidate genes')
            pass

    
    def create_gene_list_for_analysis(self, geneList_json):
        results = {}
        input_phenotypes = list(geneList_json.keys())
        for i in range(len(input_phenotypes)):
            candidate_phenotypes = list(geneList_json[input_phenotypes[i]])
            for j in range(len(candidate_phenotypes)):
                candidate_genes = list(geneList_json[input_phenotypes[i]][candidate_phenotypes[j]]['genes'])
                for candidate_gene in candidate_genes:
                    results.setdefault(candidate_gene, []).append(candidate_phenotypes[j])
    
        return results
    
    def annotation(self, dataframe, dic_candidate_genes):
        dataframe['HPO.Matched'] = dataframe['Gene.refGene'].map(dic_candidate_genes)

        return dataframe


