__author__ = "Nestor Bermudez"
__license__ = "MIT"
__version__ = "1.0.0"
__email__ = "nab6@illinois.edu"
__status__ = "Development"


import numpy as np
import pandas as pd
from collections import defaultdict


ORG = 'human'
BASE_DIR = 'data/mashup/raw/annotations/' + ORG
TERMS = 'go_human_ref_mf_terms.txt'
MIN_SIZE = 30


id2GO = {}
hierarchy = defaultdict(lambda: {'id': 0, 'children': [], 'genes': set()})


gene_go_map = pd.read_csv(BASE_DIR + '/go_human_ref_mf_adjacency.txt',
                          header=None, names=['Gene', 'GO'], sep='\t')
terms = pd.read_csv(BASE_DIR + '/' + TERMS, header=None, sep='\t').values
links = pd.read_csv(BASE_DIR + '/graph/go_mf.links', header=None,
                    names=['ParentGO', 'ChildGO'], sep='\t')


parents = {}
for row in links.values:
    parents[row[1]] = row[0]


def build_tree():
    for gene, go in gene_go_map.values:
        hierarchy[go]['genes'].add(gene)
    for parent, child in links.values:
        hierarchy[parent]['id'] = parent
        hierarchy[parent]['children'].append(child)
        genes = gene_go_map[gene_go_map['GO'] == child]['Gene'].values
        hierarchy[parent]['genes'].update(genes.tolist())


build_tree()
selected_functions = []


def select_function(root):
    if len(hierarchy[root]['genes']) > MIN_SIZE:
        return root
    else:
        if root in parents:
            return select_function(parents[root])
        return root


leaves = np.unique(gene_go_map.values[:, 1])
for leaf in leaves:
    mf = select_function(leaf)
    if mf in selected_functions:
        continue
    selected_functions.append(mf)


with open('data/mashup/raw/annotations/human/reduced_adjacency.txt', 'w') as f:
    for mf in selected_functions:
        for gene in hierarchy[mf]['genes']:
            f.write('{}\t{}\n'.format(gene, mf))

