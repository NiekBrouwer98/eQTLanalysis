import pandas as pd
from csv import writer, reader

features_path = 'C:/Users/niekb/OneDrive/Bachelor Eindproject/count_matrices_per_lane/lane_1/genes.tsv'
genes = pd.read_csv(features_path, delimiter="\t")

chr_location_path = 'C:/Users/niekb/Documents/BEP datasets/chr_locations_R.csv'
chr_location = pd.read_csv(chr_location_path, delimiter=',')

print(chr_location.shape)
print(genes.shape)
