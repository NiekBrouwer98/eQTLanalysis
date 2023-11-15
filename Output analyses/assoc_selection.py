import pandas as pd
from pandas_plink import read_rel

path = "C:/Users/niekb/plink_workspace/eqtl_ct1/"

for i in range(1, 2):
    df = pd.DataFrame
    file = read_rel('{}eQTL_analysis_ct1_run1.Gene{}.assoc.linear.adjusted'.format(path, i))

    # file_line = file.readlines()
    # for j in file_line[1:3]:
    #     column = j.split('  ')
    #     print(column[14])
        # if int(column[14]) < 0.09:
        #     print(column[2])

p

