import pandas as pd
import scipy.io
import csv
import os
import numpy as np

def get_dataframe(lanes):
    data_dir = 'C:/Users/niekb/OneDrive/Bachelor Eindproject/count_matrices_per_lane'

    data = pd.DataFrame()

    for i in range(1,lanes+1):
        lane_dir = f'{data_dir}/lane_{i}'
        mat = scipy.io.mmread(f'{lane_dir}/matrix.mtx')
        B = mat.todense()
        df = pd.DataFrame(B)

        barcodes_path = f'{lane_dir}/barcodes.tsv'
        with open(barcodes_path, newline='\n') as csvfile:
            barcodes = [f'{row[0][0:16]}_lane{i}' for row in csv.reader(csvfile, delimiter="\t")]

        features_path = f'{lane_dir}/genes.tsv'
        with open(features_path, newline='\n') as csvfile:
            gene_names = [row[1] for row in csv.reader(csvfile, delimiter="\t")]

        df.columns = barcodes
        df.index = gene_names

        data = pd.concat([data,df],axis=1)

    return data


def filtering_dataframe(df):
    # filter genes with <3 cells
    df['Total cells'] = (df[df.columns].gt(0).sum(axis=1))+(df[df.columns].lt(0).sum(axis=1))
    df_3 = df.loc[df['Total cells']>3]
    df_3 = df_3.drop('Total cells', axis=1)

    # filter cells with >5% MT-genes
    df_3_mt = df_3.transpose()
    df_3_mt['percent mito'] = ((df_3_mt['MT-ND1']+df_3_mt['MT-ND2']+df_3_mt['MT-CO1']+df_3_mt['MT-CO2']+df_3_mt['MT-ATP8']+
                             df_3_mt['MT-ATP6']+df_3_mt['MT-CO3']+df_3_mt['MT-ND3']+df_3_mt['MT-ND4L']+df_3_mt['MT-ND4']+
                             df_3_mt['MT-ND5']+df_3_mt['MT-ND6']+df_3_mt['MT-CYB'])/df_3_mt.sum(axis=1))*100
    df_3_mt = df_3_mt.loc[df_3_mt['percent mito']<5]
    df_3_mt = df_3_mt.drop('percent mito',axis=1)

    #filter cells with >3500 genes
    df_3_mt['Total genes'] = (df_3_mt[df_3_mt.columns].gt(0).sum(axis=1))+(df_3_mt[df_3_mt.columns].lt(0).sum(axis=1))
    df_3_mt_3500 = df_3_mt.loc[df_3_mt['Total genes']<3500]
    filtered_dataframe = df_3_mt_3500.drop('Total genes', axis=1)
    filtered_dataframe = filtered_dataframe.transpose()

    shapes = [df_3.shape, df_3_mt.shape, filtered_dataframe.shape]

    return filtered_dataframe


def demultiplex_dataframe(df):
    #filter doublet or inconclusive cells (not in cell_barcodes file)
    singlets_path = 'C:/Users/niekb/Documents/BEP datasets/cell_barcodes/pilot3_persons.tsv'
    with open(singlets_path, newline='\n') as csvfile:
        singlets = np.array([row[0] for row in csv.reader(csvfile, delimiter="\t")])

    os.chdir('C:/Users/niekb/Documents/BEP datasets/demultiplexed_lanes')

    filtered_df = pd.DataFrame()

    for sample in df.columns:
        if sample in singlets:
            filtered_df[f'{sample}'] = df[f'{sample}']

    features_path = 'C:/Users/niekb/OneDrive/Bachelor Eindproject/count_matrices_per_lane/lane_1/genes.tsv'
    with open(features_path, newline='\n') as csvfile:
        gene_names = [row[1] for row in csv.reader(csvfile, delimiter="\t")]

    filtered_df.to_csv('demultiplexed_alllanes.csv')
    # filtered_df.to_pickle('demultiplexed_alllanes.pkl')


def get_labels():
    labels_dir = 'C:/Users/niekb/Documents/BEP datasets/barcodes_to_cell_types'
    labels_path = f'{labels_dir}/barcodes_to_cell_types.tsv'
    labels_df = pd.read_csv(labels_path, sep='\t', header=0)

    return labels_df


def get_barcodes_to_labels(data, labels):

    os.chdir('C:/Users/niekb/Documents/BEP datasets/labels_lanes')

    lst = []
    for sample_name in data.columns:
        for i in range(len(labels)):
            if sample_name == labels.iloc[i,0]:
                lst.append(labels.iloc[i,1])
                break
        else:
            lst.append('no_label')

    df = pd.DataFrame(lst)
    df.to_pickle('labels_alllanes.pkl')


def normalize_data(df):
    result = pd.DataFrame(index=df.index.copy())
    for feature_name in df.columns:
        total = df[feature_name].sum()
        result[feature_name] = (df[feature_name]/total)*10000
    return result




data = pd.read_pickle('C:/Users/niekb/Documents/BEP datasets/demultiplexed_lanes/demultiplexed_alllanes.pkl')
normalized_data = normalize_data(data)
normalized_data.to_csv('C:/Users/niekb/Documents/BEP datasets/normalized_lanes/normalized_demultiplexed_alllanes.csv')
normalized_data.to_pickle('C:/Users/niekb/Documents/BEP datasets/normalized_lanes/normalized_demultiplexed_alllanes.pkl')

