
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
# dataset_df['All32'] = SixBefore
dataset_df['SixAfter'] = SixBefore
# print(dataset_df)
# dataset_df

# print(dataset_df.sort_values('chrom'))
# print(dataset_df.loc[dataset_df['chrom'] == 'chr1'])
def ReverseAndChange(OrigSeq):
    OrigSeq = OrigSeq[::-1]
    ReverseSeq=""
    for i in range(len(OrigSeq)):
        if OrigSeq[i] == 'A':
            ReverseSeq += 'T'
        elif OrigSeq[i] == 'T':
            ReverseSeq += 'A'
        elif OrigSeq[i] == 'G':
            ReverseSeq += 'C'
        elif OrigSeq[i] == 'C':
            ReverseSeq += 'G'
        else:
            continue
    return ReverseSeq
# print(len(fasta_sequences))
l = fasta_sequences
# for i in range(26):
ListOfIds = ["NC_00000"+str(i) for i in range(1,10)] + ["NC_0000"+str(i) for i in range(10,27)]
ListOfChrs = ["chr"+str(i) for i in range(1,27)]
ChromStartVector = dataset_df['chromStart']
chromEndVector = dataset_df['chromEnd']
print(ListOfIds)
def CheckSequence(seq,index):
    if isinstance(seq, float):
        print("found ",seq,"  in index =",index)
        return True
    if len(str(seq)) > 23:
        print("found seq in len() = ",len(str(seq))," in index = ",index)
        return True
for i in fasta_sequences:
    # print(i.id[0:2])
    if i.id[0:2] == "NC":
        if i.id[0:9] in ListOfIds:
            # TempDataSet = dataset_df.loc[dataset_df['chrom'] == ListOfChrs[int(i.id[8])]]
            print("Starting adding 6 nuc before and after in %s", ListOfChrs[int(i.id[7:9]) - 1])
            for j in dataset_df.loc[dataset_df['chrom'] == ListOfChrs[int(i.id[7:9])-1]].index:
                if CheckSequence(dataset_df['offtarget_sequence'][j],j):
                    continue
                TempSeqLen = len(dataset_df['offtarget_sequence'][j])
                if dataset_df['strand'][j] == '+':
                    if dataset_df['offtarget_sequence'][j] != i.seq[ChromStartVector[j]:chromEndVector[j]].upper() and "-" not in dataset_df['offtarget_sequence'][j]:
                        print(j)
                        print(dataset_df['offtarget_sequence'][j])
                        print(i.seq[ChromStartVector[j]:chromEndVector[j]])
                        raise ValueError
                    temp35 = i.seq[ChromStartVector[j] - 6:chromEndVector[j] + 6]
                    dataset_df['SixAfter'][j] = str(temp35[TempSeqLen+5:TempSeqLen+11]).upper()
                    dataset_df['SixBefore'][j] = str(temp35[0:6]).upper()
                    # dataset_df['All32'][j] = str(temp35).upper()
                else:
                    orig_str = str(i.seq[ChromStartVector[j] - 6:chromEndVector[j] + 6]).upper()
                    temp35 = ReverseAndChange(orig_str)
                    if dataset_df['offtarget_sequence'][j] != temp35[6:TempSeqLen+6] and "-" not in dataset_df['offtarget_sequence'][j]:
                        print(j)
                        print(dataset_df['offtarget_sequence'][j])
                        print(temp35[6:TempSeqLen+6])
                        print(orig_str)
                        raise ValueError
                    dataset_df['SixAfter'][j] = str(temp35[TempSeqLen+5:TempSeqLen+11])
                    dataset_df['SixBefore'][j] = str(temp35[0:6])
                    # dataset_df['All32'][j] = str(temp35)
            if i.id[0:9] == "NC_000023":
                for j in dataset_df.loc[dataset_df['chrom'] == 'chrX'].index:
                    if CheckSequence(dataset_df['offtarget_sequence'][j],j):
                        continue
                    TempSeqLen = len(dataset_df['offtarget_sequence'][j])
                    if dataset_df['strand'][j] == '+':
                        if dataset_df['offtarget_sequence'][j] != i.seq[ChromStartVector[j]:chromEndVector[j]].upper() and "-" not in dataset_df['offtarget_sequence'][j]:
                            print(j)
                            print(dataset_df['offtarget_sequence'][j])
                            print(i.seq[ChromStartVector[j]:chromEndVector[j]])
                            raise ValueError
                        temp35 = i.seq[ChromStartVector[j] - 6:chromEndVector[j] + 6]
                        dataset_df['SixAfter'][j] = str(temp35[TempSeqLen+5:TempSeqLen+11]).upper()
                        dataset_df['SixBefore'][j] = str(temp35[0:6]).upper()
                        # dataset_df['All32'][j] = str(temp35).upper()
                    else:
                        orig_str = str(i.seq[ChromStartVector[j] - 6:chromEndVector[j] + 6]).upper()
                        temp35 = ReverseAndChange(orig_str)
                        if dataset_df['offtarget_sequence'][j] != temp35[6:TempSeqLen+6] and "-" not in dataset_df['offtarget_sequence'][j]:
                            print(j)
                            print(dataset_df['offtarget_sequence'][j])
                            print(temp35[6:29])
                            print(orig_str)
                            raise ValueError
                        dataset_df['SixAfter'][j] = str(temp35[TempSeqLen+5:TempSeqLen+6])
                        dataset_df['SixBefore'][j] = str(temp35[0:6])
                        # dataset_df['All32'][j] = str(temp35)
            dataset_df.to_excel( general_utilities.DATASETS_PATH + "CHANGE-seqNew.xlsx" )

# print(dataset_df)
# print(dataset_df['SixAfter'])
# print(dataset_df['SixBefore'])
# print(dataset_df['All32'])
