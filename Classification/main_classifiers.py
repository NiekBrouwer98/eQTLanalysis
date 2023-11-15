import pandas as pd
from wijst_SVM_classification import Wijst_SVM
from wijst_SVMrejection_classification import Wijst_SVMrejection
from wijst_data_preprocessing import get_dataframe, filtering_dataframe, get_labels


'''DataPath, LabelsPath, CV_RDataPath, OutputDir, GeneOrderPath = "", NumGenes = 0)'''

# #10xv2 PBMC dataset:
Data_folder = 'C:/Users/niekb/Documents/BEP datasets/scRNAseq_Benchmark_datasets/Inter-dataset/PbmcBench/10Xv2'
DataPath  = f'{Data_folder}/10Xv2_pbmc1.csv'
LabelsPath = f'{Data_folder}/10Xv2_pbmc1Labels.csv'
# CV_RDataPath = f'{Data_folder}/10Xv2_pbmc1_folds.RData'
# outputdir = 'C:/Users/niekb/Documents/BEP datasets/Classification_results/SVMrejection_10xv2'
#
data_10xv2 = pd.read_csv(DataPath,index_col=0,sep=',')
labels_10xv2 = pd.read_csv(LabelsPath, header=0,index_col=None, sep=',')


# Wijst data:
data = pd.read_pickle('C:/Users/niekb/Documents/BEP datasets/demultiplexed_lanes/demultiplexed_alllanes.csv')
data = data.transpose()
labels = pd.read_pickle('C:/Users/niekb/Documents/BEP datasets/labels_lanes/labels_alllanes.pkl')
labels = labels[1:]
output = 'C:/Users/niekb/Documents/BEP datasets/Classification_results/SVM_wijst'
output_rejection = 'C:/Users/niekb/Documents/BEP datasets/Classification_results/SVMrejection_wijst'

SVM = Wijst_SVM(data, labels, output)
SVM.run_SVM_crossval()
SVM.run_SVM_trainset(data_10xv2, labels_10xv2)

print(data_10xv2.head())
print(data.head())
