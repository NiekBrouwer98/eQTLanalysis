import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix
import seaborn as sn
import matplotlib.pyplot as plt

class output_analysis:

    def __init__(self, output_dir, prediction_name, true_name):
        self.output_dir = output_dir
        self.prediction = pd.read_csv(f'{output_dir}/{prediction_name}')
        self.true = pd.read_csv(f'{output_dir}/{true_name}')

    def count(self):
        self.prediction['comparison'] = np.where(self.prediction['0']==self.true['0'],True, False)
        return self.prediction.groupby('comparison').count()

    def count_rejected(self):
        self.prediction['rejected'] = np.where(self.prediction['0']=='Unknown', True, False)
        return self.prediction.groupby('rejected').count()

    def confusion_matrix(self, labels):
        print('do not forget to add plt.show() after executing this function.')
        array = confusion_matrix(self.true, self.prediction, labels)

        df_cm = pd.DataFrame(array, range(len(array)), range(len(array)))
        # plt.subplots(figsize=(10,10))
        sn.set(font_scale=1.4)
        return sn.heatmap(df_cm, annot=True, annot_kws={"size": 6}, xticklabels=labels, yticklabels=labels)

output_dir = 'C:/Users/niekb/Documents/BEP datasets/Classification_results/SVM_wijst'
prediction_name = 'crossval_SVM_pred_labels.csv'
true_name = 'crossval_SVM_true_labels.csv'
labels = pd.read_pickle('C:/Users/niekb/Documents/BEP datasets/labels_lanes/labels_alllanes.pkl')
labels = labels[0].unique()

SVM_crossval = output_analysis(output_dir, prediction_name, true_name)
SVM_crossval.confusion_matrix(labels)
plt.show()

