from __future__ import print_function

__author__ = "Nestor Bermudez"
__license__ = "MIT"
__version__ = "1.0.0"
__email__ = "nab6@illinois.edu"
__status__ = "Development"


import json
import networkx as nx
import numpy as np
import pandas as pd
import random
import tensorflow as tf
from networkx.readwrite import json_graph
from sklearn.preprocessing import LabelEncoder
from sklearn.externals.joblib import dump


assert(nx.__version__ <= (1, 11))


flags = tf.flags
FLAGS = flags.FLAGS

flags.DEFINE_string("networks_dir", "data/mashup/raw/networks/", "")
flags.DEFINE_string("org", "human", "")
flags.DEFINE_string("annotations_dir", "data/mashup/raw/annotations/", "")
flags.DEFINE_string("output_dir", "data/graphsage/from_mashup/human2/", "")


WALK_LEN = 5
N_WALKS = 100


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


def get_all_genes():
    return pd.read_csv(FLAGS.networks_dir + FLAGS.org + '/human_string_genes.txt', header=None, names=['gene'])


def read_adjacency_file(filepath):
    data = pd.read_csv(filepath, sep='\t', header=None, names=['u', 'v', 'w'], dtype={'u': int, 'v': int, 'w': float})
    data['u'] = data['u'] - 1
    data['v'] = data['v'] - 1

    return data


def read_annotations(filepath):
    data = pd.read_csv(filepath, sep='\t', header=None, names=['gene', 'GO'])
    data['gene'] = data['gene'] - 1
    label_encoder = LabelEncoder()
    data['GO'] = label_encoder.fit_transform(data['GO'])
    return data, label_encoder


def gene_class_map(gene, annotations, n_classes):
    indices = annotations[annotations['gene'] == gene]['GO'].values
    class_map = np.zeros(n_classes, dtype=int)
    class_map[indices] = 1
    return class_map.tolist()


if __name__ == '__main__':
    import os
    G = nx.Graph()
    id_map = {}
    class_map = {}

    genes = get_all_genes()
    labels, label_encoder = read_annotations(
        FLAGS.annotations_dir + FLAGS.org + '/reduced_adjacency.txt')
    n_classes = label_encoder.classes_.size
    n_genes = genes.size
    n_networks = 6

    network_id = 0
    node_id = 0
    node_id_set = set()
    for file in os.listdir(FLAGS.networks_dir + FLAGS.org):
        if not file.endswith('adjacency.txt'):
            continue
        for i, gene in enumerate(range(n_genes)):
            node_id = '{}_{}'.format(network_id, gene)
            node_id_set.add(node_id)
            val = 'cooccurence' in file
            G.add_node('{}_{}'.format(network_id, gene), test=False, val=val)
            class_map[node_id] = gene_class_map(gene, labels, n_classes)
            id_map[node_id] = network_id * n_genes + i
        network_id += 1
    network_id = 0

    for file in os.listdir(FLAGS.networks_dir + FLAGS.org):
        if not file.endswith('adjacency.txt'):
            continue
        print('Processing {} graph'.format(file))
        full_data = read_adjacency_file(FLAGS.networks_dir + FLAGS.org + '/' + file)
        for data in np.array_split(full_data, full_data.shape[0] / 20000 + 1):
            for i, (u, v, w) in enumerate(data.values):
                u = int(u)
                v = int(v)
                print('Processing link {0:6d} out of {1:6d}'.format(i, data.values.shape[0]), end='\r')
                u_id = '{}_{}'.format(network_id, u)
                v_id = '{}_{}'.format(network_id, v)
                G.add_edge(u_id, v_id, weight=w)
        network_id += 1

    print('=' * 120)
    print('Saving files...')
    walks = run_random_walks(G, G.nodes())

    if not os.path.exists(FLAGS.output_dir):
        os.makedirs(FLAGS.output_dir)

    with open(FLAGS.output_dir + '/ppi-G.json', 'w') as f:
        json.dump(json_graph.node_link_data(G), f)

    with open(FLAGS.output_dir + '/ppi-id_map.json', 'w') as f:
        json.dump(id_map, f)

    with open(FLAGS.output_dir + '/ppi-class_map.json', 'w') as f:
        json.dump(class_map, f)

    with open(FLAGS.output_dir + '/ppi-walks.txt', 'w') as f:
        for walk in walks:
            f.write('\t'.join(walk) + '\n')

    dump(label_encoder, FLAGS.output_dir + '/ppi-label_encoder.pkl')
    print('DONE')
