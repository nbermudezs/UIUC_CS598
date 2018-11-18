# Comparison of embedding methods for gene function prediction

In this repo you will find the code necessary to replicate the results presented
in the CS598SS: Advance Bioinformatics class at UIUC in Fall 2018.

## Setup

This codebase uses Python 3.x, for best results we recommend using `virtualenv`.
```bash
virtualenv env --python=$(which python3)
source env/bin/activate
pip install -r requirements
```
will create a new environment and setup all necessary dependencies for running the models.


The datasets used for this work can be downloaded from this Box [folder](https://uofi.box.com/s/5dyilux65ntviye1o1f9e5miti9cn01c).
Please take all the folder structure and place it in the `data` folder at the root of this repo.

Once the data is in your machine you can decide whether to preprocess it yourself or use the preprocessed files.

The preprocessing itself is different for each of the methods evaluated in this work.
For GraphSAGE the data is expected to be in a JSON representation of a NetworkX graph
accompanied with corresponding files to map each node to its id and classes as well as features.
To do this prepprocessing we created the `preprocess/graphsage.py` script.

Run
```bash
python preprocess/graphsage.py
python preprocess/mashup.py
python preprocess/gat.py
python preprocess/yagcn.py
```

to perform all the necessary preprocessing. Beware, this may take a while!
