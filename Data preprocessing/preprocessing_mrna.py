import pandas as pd
import numpy as np

# database_path = 'C:/Users/niekb/Documents/BEP datasets/Montgomery datasets/raw_files'
#
# pheno_df = pd.read_csv(f'{database_path}/RNASEQ60_array_rep_expr.txt', sep='\t')
# snp_df = pd.read_csv(f'{database_path}/RNASEQ60_array_rep_snps.full.txt', sep='\t')


def pheno_files(df):
    for chr_number in range(1, 23):
        chr_df = pd.DataFrame(columns=df.columns)
        chr_df = chr_df.append(df.loc[df['chr'] == chr_number])
        chr_df.to_csv(f'C:/Users/niekb/Documents/BEP datasets/Montgomery datasets/phenotype_files/phenotypefile_chr{chr_number}.txt',sep='\t', index=False)

def snp_files(df):
    for chr_number in range(1, 23):
        chr_df = pd.DataFrame(columns=df.columns)
        chr_df = chr_df.append(df.loc[df['chromosome'] == chr_number])
        chr_df.to_csv(f'C:/Users/niekb/Documents/BEP datasets/Montgomery datasets/genotype_files/genotypefile_chr{chr_number}.txt', sep='\t', index=False)

def genotype_file_to_number(filepath):
    for i in range(1,2):
        df = pd.read_csv(f'{filepath}/genotypefile_chr{i}.txt', sep='\t')
        for index, row in df.iterrows():
            marker = row[3].replace('/', '')
            for j in df.columns[4:64]:
                if row[j] == marker:
                    df[row[j]] = '0'
                elif row[j][0] == marker[0] or row[j][1] == marker[1]:
                    df[row[j]] = '1'
                else:
                    df[row[j]] = '2'

        df.to_csv(f'C:/Users/niekb/Documents/BEP datasets/Montgomery datasets/genotype_files/number_genotypefile_chr{i}.txt', sep='\t', index=False)

genotype_file_to_number('C:/Users/niekb/Documents/BEP datasets/Montgomery datasets/genotype_files')



# import sys
#
# infile_name = sys.argv[1]
#
# pedDict = {
# "0" : "R R",
# "1" : "R A",
# "2" : "A A",
# "3" : "0 0"
# }
#
# def convertToPlink(infile_name):
#     with open(infile_name, 'r') as infile:
#         header = infile.readline().rstrip().split()
#         chromosome = infile.readline().split()[0]
#         with open('chr_' + chromosome + '.map', 'w') as mapfile:
#             for POS in header[1:]:
#                 mapfile.write("\t".join([chromosome, chromosome+"_"+POS+"_SNP", "0", POS])+"\n")
#     with open(infile_name, 'r') as infile:
#         with open('chr_' + chromosome + '.ped', 'w') as pedfile:
#             id_index = 0
#             for line in infile:
#                 if not line.startswith("CHR"):
#                     id_index += 1
#                     IID = "ind_" + str(id_index)
#                     line = line.rstrip().split()
#                     pedfile.write(" ".join([IID, IID, "0", "0", "0", "-9"]+[pedDict[genotype] for genotype in line][1:])+ "\n")
#
#
# convertToPlink(infile_name)
