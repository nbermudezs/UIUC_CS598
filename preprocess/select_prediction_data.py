__author__ = "Nestor Bermudez"
__license__ = "MIT"
__version__ = "1.0.0"
__email__ = "nab6@illinois.edu"
__status__ = "Development"


import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.externals.joblib import dump


labels = {
    'human': 'data/mashup/raw/annotations/human/reduced_adjacency.txt',
    'yeast': 'data/mashup/raw/annotations/yeast/yeast_mips_level2_adjacency.txt'
}


features = {
    'human': '~/Downloads/mashup_features.csv',
    'yeast': ''
}


ORG = 'human'


feats = pd.read_csv(features[ORG], header=None).values.T
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


new_features.dump('data/processed/mashup.human-features.npy')
new_labels.dump('data/processed/mashup.human-labels.npy')

np.array(genes).dump('data/processed/mashup.human-genes.npy')
dump(label_encoder, 'data/processed/mashup.human-labelencoder.pkl')
