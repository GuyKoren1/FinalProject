from pathlib import Path
import pandas as pd
from Bio.Seq import Seq
from SysEvalOffTarget_src import general_utilities
import numpy as np
from Bio import SeqIO
input_file = general_utilities.HOME_DIR + "GRCh38_latest_genomic.fna"
output_file = "GRCH38222222.txt"
genom = ""
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
real_fasta_sequences = [0 for i in range(26)]
count = 0
# print(fasta_sequences)
dataset_df = pd.read_excel(general_utilities.CHANGE_SEQ_PATH)
SixBefore = [0 for i in range(len(dataset_df.index))]
dataset_df['SixBefore'] = SixBefore
dataset_df['All32'] = SixBefore
dataset_df['SixAfter'] = SixBefore
print(dataset_df)
# dataset_df
# print(dataset_df.sort_values('chrom'))
print(dataset_df.loc[dataset_df['chrom'] == 'chr1'])

# print(len(fasta_sequences))
l = fasta_sequences
# for i in range(26):
ListOfIds = ["NC_00000"+str(i) for i in range(1,10)] + ["NC_0000"+str(i) for i in range(10,27)]
ListOfChrs = ["chr"+str(i) for i in range(1,27)]
ChromStartVector = dataset_df['chromStart']
chromEndVector = dataset_df['chromEnd']
print(ListOfIds)
for i in fasta_sequences:
    # print(i.id[0:2])
    if i.id[0:2] == "NC":
        if i.id[0:9] in ListOfIds:
            # TempDataSet = dataset_df.loc[dataset_df['chrom'] == ListOfChrs[int(i.id[8])]]
            print(i.id[7:9])
            print(ListOfChrs[int(i.id[7:9])-1])
            for j in dataset_df.loc[dataset_df['chrom'] == ListOfChrs[int(i.id[7:9])-1]].index:
                if dataset_df['strand'][j] == '+':
                    # print("strand = " + dataset_df['strand'][j])
                    # print("off target from table:")
                    # print(dataset_df['offtarget_sequence'][j])
                    # print("seq from genome by start and end inexes:")
                    print(i.seq[ChromStartVector[j]:chromEndVector[j]])
                    print("-----------------------------------")
                    if dataset_df['offtarget_sequence'][j] != i.seq[ChromStartVector[j]:chromEndVector[j]]:
                        print(j)
                        print(dataset_df['offtarget_sequence'][j])
                        print(i.seq[ChromStartVector[j]:chromEndVector[j]])
                        raise ValueError
                    dataset_df['SixAfter'][j] = str(i.seq[chromEndVector[j]:chromEndVector[j] + 6])
                    dataset_df['SixBefore'][j] = str(i.seq[ChromStartVector[j] - 6:ChromStartVector[j]])
                    dataset_df['All32'][j] = str(i.seq[ChromStartVector[j]-6:chromEndVector[j] + 6])
            dataset_df.to_excel( general_utilities.DATASETS_PATH + "New" + "CHANGE-seq.xlsx" )

# print(dataset_df)
# print(dataset_df['SixAfter'])
# print(dataset_df['SixBefore'])
# print(dataset_df['All32'])
