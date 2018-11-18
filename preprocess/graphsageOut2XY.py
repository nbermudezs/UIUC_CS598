__author__ = "Nestor Bermudez"
__license__ = "MIT"
__version__ = "1.0.0"
__email__ = "nab6@illinois.edu"
__status__ = "Development"


import numpy as np
import os
import pandas as pd
import tensorflow as tf
from sklearn.preprocessing import LabelEncoder
from sklearn.externals.joblib import dump


flags = tf.flags
FLAGS = flags.FLAGS

flags.DEFINE_string("embeddings_filepath", "", "")
flags.DEFINE_integer("n_graphs", 6, "")
flags.DEFINE_boolean("average", False, "")
flags.DEFINE_string("output_dir", "data/processed", "")
flags.DEFINE_string("prefix", "graphsage-mean", "")


labels = {
    'human': 'data/mashup/raw/annotations/human/reduced_adjacency.txt',
    'yeast': 'data/mashup/raw/annotations/yeast/yeast_mips_level2_adjacency.txt'
}


ORG = 'human'


def main(_):
    raw_data = np.load(FLAGS.embeddings_filepath)
    if FLAGS.average:
        n_nodes = raw_data.shape[0] / FLAGS.n_graphs
        new_shape = (FLAGS.n_graphs, n_nodes, raw_data.shape[1])
        reshaped = raw_data.reshape(new_shape)
        feats = np.mean(reshaped, axis=0)
    else:
        n_nodes = raw_data.shape[0]
        feats = raw_data
    assert(feats.shape == (n_nodes, raw_data.shape[1]))
    print('Feature matrix shape: ', feats.shape)

    y = pd.read_csv(labels[ORG], header=None, names=['gene', 'mf'], sep='\t')
    print('# links between genes and MFs: ', y.values.shape[0])

    label_encoder = LabelEncoder()
    label_encoder.fit(y.values[:, 1])
    label_size = len(label_encoder.classes_)

    new_features = []
    new_labels = []
    genes = []
    for gene, mfs in y.groupby('gene'):
        genes.append(gene)
        gene_feats = feats[gene - 1, :]
        gene_y = np.zeros(label_size)
        gene_y[label_encoder.transform(mfs['mf'].values)] = 1
        new_features.append(gene_feats)
        new_labels.append(gene_y)

    new_features = np.vstack(new_features)
    new_labels = np.vstack(new_labels)

    print('Shape of selected gene matrix: ', new_features.shape)
    print('Shape of selected label matrix: ', new_labels.shape)

    if not os.path.exists(FLAGS.output_dir):
        os.makedirs(FLAGS.output_dir)

    new_features.dump('{}/{}.{}-features.npy'.format(FLAGS.output_dir, FLAGS.prefix, ORG))
    new_labels.dump('{}/{}.{}-labels.npy'.format(FLAGS.output_dir, FLAGS.prefix, ORG))

    np.array(genes).dump('{}/{}.{}-genes.npy'.format(FLAGS.output_dir, FLAGS.prefix, ORG))
    dump(label_encoder, '{}/{}.{}-labelencoder.pkl'.format(FLAGS.output_dir, FLAGS.prefix, ORG))


if __name__ == '__main__':
    main(0)