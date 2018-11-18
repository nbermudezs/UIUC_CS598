__author__ = "Nestor Bermudez"
__license__ = "MIT"
__version__ = "1.0.0"
__email__ = "nab6@illinois.edu"
__status__ = "Development"


import json
import networkx as nx
import numpy as np
import pandas as pd
import pickle
import random
import sys
from collections import defaultdict
from networkx.readwrite import json_graph
from sklearn.preprocessing import LabelEncoder

N_CLASSES = 100
BASE_DIR = 'data/ppi_nab6/'
WALK_LEN = 5
N_WALKS = 50
N_FEATURES = 50


def all_genes():
    print('Identifying all genes across networks')
    data = pd.read_csv(
        BASE_DIR + 'raw/interactions.txt',
        keep_default_na=False,
        delimiter='\t').values
    genes = np.unique(data)[1:]
    for file in ['HomoSapiens_binary_hq.txt', 'HomoSapiens_cocomp_hq.txt']:
        data = pd.read_csv(BASE_DIR + '/raw/' + file, usecols=['Gene_A', 'Gene_B'], sep='\t')
        genes = np.unique(np.concatenate((genes, ['test_' + x for x in np.unique(data.values)])))

    with open(BASE_DIR + '/onto/c5.bp.v6.2.symbols.gmt') as f:
        for line in f:
            parts = line.replace('\n', '').split('\t')
            if len(parts[2:]) > 10:
                genes = np.unique(np.concatenate((genes, parts[2:])))
    import pdb; pdb.set_trace()
    return genes


def add_test_graphs(G, encoder):
    print('Adding testing graphs')
    for file in ['HomoSapiens_binary_hq.txt', 'HomoSapiens_cocomp_hq.txt']:
        data = pd.read_csv(BASE_DIR + '/raw/' + file, usecols=['Gene_A', 'Gene_B'], sep='\t')
        us = encoder.transform(['test_' + x for x in data.values[:, 0]])
        vs = encoder.transform(['test_' + x for x in data.values[:, 1]])
        for u, v in zip(us, vs):
            G.add_node(u, test=True, val=False)
            G.add_node(v, test=True, val=False)
            G.add_edge(u, v)
    return G


def add_train_graphs(G, encoder):
    print('Adding more training graphs')
    data = pd.read_csv(BASE_DIR + '/raw/Multinet.interactions.network_presence.txt', sep='\t')
    networks = ['METABOLIC', 'PPI', 'GENETIC', 'SIGNALING']
    for net in networks:
        print('Processing {} network'.format(net))
        network = data[data[net] == 1]['INTERACTION_NAME'].str.split('_').values
        for i, row in enumerate(network):
            sys.stdout.write('{} out of {}\r'.format(i, network.shape[0]))
            sys.stdout.flush()
            u, v = encoder.transform(row)
            G.add_node(u, test=False, val=False)
            G.add_node(v, test=False, val=False)
            G.add_edge(u, v)
    return G


def run_random_walks(G, nodes, num_walks=N_WALKS):
    pairs = []
    for count, node in enumerate(nodes):
        if G.degree(node) == 0:
            continue
        for i in range(num_walks):
            curr_node = node
            for j in range(WALK_LEN):
                next_node = random.choice(G.neighbors(curr_node))
                # self co-occurrences are useless
                if curr_node != node:
                    pairs.append((node,curr_node))
                curr_node = next_node
        if count % 1000 == 0:
            print("Done walks for", count, "nodes")
    return pairs


def create_features(id_map, encoder):
    map = defaultdict(lambda: [])
    u_feats = []
    with open(BASE_DIR + '/onto/c7.all.v6.2.symbols.gmt') as f:
        for i, line in enumerate(f):
            if i == N_FEATURES:
                break
            parts = line.replace('\n', '').split('\t')
            label = parts[0]
            u_feats.append(label)
            for gene in parts[2:]:
                map[gene].append(label)
    feat_encoder = LabelEncoder()
    feat_encoder.fit(u_feats)
    N = np.array(len(id_map), dtype=int).max() + 1
    ndim = len(u_feats)

    feats = np.zeros((N, ndim))
    for key, gene_id in id_map.items():
        gene_symbol = encoder.inverse_transform(gene_id)
        labels = np.zeros(ndim)
        if len(map[gene_symbol]) > 0:
            labels[feat_encoder.transform(map[gene_symbol])] = 1
        feats[gene_id, :] = labels
    return feats


if __name__ == '__main__':
    nodes = all_genes()
    classes = {}
    with open(BASE_DIR + '/onto/c5.bp.v6.2.symbols.gmt') as f:
        for line in f:
            parts = line.replace('\n', '').split('\t')
            gene_fn = parts[0]
            genes = parts[2:]
            if len(genes) > 10:
                classes[gene_fn] = genes
    N = len(classes.keys())
    selected = np.random.choice(range(N), N_CLASSES)

    fns = list(classes.keys())
    fns = [fns[idx] for idx in selected]

    annotated_genes = list(classes.values())
    annotated_genes = [annotated_genes[idx] for idx in selected]
    
    print('Creating gene function encoder')
    gene_fn_encoder = LabelEncoder()
    fn_ids = gene_fn_encoder.fit_transform(fns)
    with open(BASE_DIR + '/fn_encoder.pkl', 'wb') as f:
        pickle.dump(gene_fn_encoder, f)

    data = pd.read_csv(
        BASE_DIR + 'raw/interactions.txt',
        keep_default_na=False,
        delimiter='\t').values

    print('Creating gene encoder')
    encoder = LabelEncoder()
    encoder.fit(nodes)
    with open(BASE_DIR + '/gene_encoder.pkl', 'wb') as f:
        pickle.dump(encoder, f)

    print('Creating first graph')
    G = nx.Graph()
    for i in range(data.shape[0]):
        sys.stdout.write('i={}\r'.format(i))
        sys.stdout.flush()
        if len(data[i, 0]) == 0 or len(data[i, 1]) == 0:
            continue
        u, v = encoder.transform(data[i, :])
        G.add_node(u, test=False, val=False)
        G.add_node(v, test=False, val=False)
        G.add_edge(u, v)

    G = add_train_graphs(G, encoder)
    G = add_test_graphs(G, encoder)

    data = json_graph.node_link_data(G)
    with open(BASE_DIR + '/ppi-G.json', 'w') as f:
        json.dump(data, f)

    class_map = defaultdict(lambda: np.zeros(N_CLASSES))
    for fn, gene_set in zip(fn_ids, annotated_genes):
        for gene in encoder.transform(gene_set):
            class_map[gene][fn] = 1

    id_map = {}
    for i, node in enumerate(G.nodes()):
        id_map[node] = i

    for key, val in class_map.items():
        class_map[key] = val.astype(int).tolist()

    with open(BASE_DIR + '/ppi-class_map.json', 'w') as f:
        json.dump(class_map, f)

    with open(BASE_DIR + '/ppi-id_map.json', 'w') as f:
        json.dump(id_map, f)

    train_nodes = [n for n in G.nodes() if not G.node[n]["test"]]
    train_G = G.subgraph(train_nodes)
    pairs = run_random_walks(train_G, train_nodes)
    with open(BASE_DIR + 'ppi-walks.txt', 'w') as fp:
        fp.write("\n".join([str(p[0]) + "\t" + str(p[1]) for p in pairs]))

    features = create_features(id_map, encoder)
    features.dump(BASE_DIR + '/ppi-feats.npy')