import pandas as pd
from SysEvalOffTarget_src import general_utilities
from Bio import SeqIO
import warnings
warnings.filterwarnings("ignore")

def ReverseAndChange(OrigSeq):
    OrigSeq = OrigSeq[::-1].upper()
    ReverseSeq = ""
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


def CheckSequence(seq, index):
    if isinstance(seq, float):
        # print("found ", seq, "  in index =", index)
        return True


def AddingSixBeforeAndSixAfterToExcel(GenomePath, OrigExcelPath, SixBeforeName, SixAfterName, ChromStartName, ChromEndName, ResultExcelPath):
    input_file = GenomePath
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
    dataset_df = pd.read_excel(OrigExcelPath)
    SixBefore = [0 for i in range(len(dataset_df.index))]
    dataset_df[SixBeforeName] = SixBefore
    # dataset_df['All32'] = SixBefore # TODO put in comment
    dataset_df[SixAfterName] = SixBefore
    ListOfIds = ["NC_00000" + str(i) for i in range(1, 10)] + ["NC_0000" + str(i) for i in range(10, 25)]
    ListOfChrs = ["chr" + str(i) for i in range(1, 25)]
    ChromStartVector = dataset_df[ChromStartName]
    chromEndVector = dataset_df[ChromEndName]

    for i in fasta_sequences:
        # print(i.id[0:2])
        if i.id[0:2] == "NC":
            if i.id[0:9] in ListOfIds:
                print("Starting adding 6 nuc before and after in ", ListOfChrs[int(i.id[7:9]) - 1])
                for j in dataset_df.loc[dataset_df['chrom'] == ListOfChrs[int(i.id[7:9]) - 1]].index:
                    if CheckSequence(dataset_df['offtarget_sequence'][j], j):
                        # print(dataset_df.loc[[j]])
                        continue
                    TempSeqLen = len(dataset_df['offtarget_sequence'][j])
                    if dataset_df['strand'][j] == '+':
                        if dataset_df['offtarget_sequence'][j] != i.seq[ChromStartVector[j]:chromEndVector[j]].upper() and "-" not in dataset_df['offtarget_sequence'][j]:
                            print(j)
                            print(dataset_df['offtarget_sequence'][j])
                            print(i.seq[ChromStartVector[j]:chromEndVector[j]])
                            raise ValueError
                        temp35 = i.seq[ChromStartVector[j] - 6:chromEndVector[j] + 6]
                        dataset_df.loc[j,SixAfterName] = str(temp35[TempSeqLen + 6:TempSeqLen + 12]).upper()
                        dataset_df.loc[j,SixBeforeName] = str(temp35[0:6]).upper()
                        # dataset_df['All32'][j] = str(temp35).upper()# TODO put in comment
                    else:
                        orig_str = str(i.seq[ChromStartVector[j] - 6:chromEndVector[j] + 6]).upper()
                        temp35 = ReverseAndChange(orig_str)
                        if dataset_df['offtarget_sequence'][j] != temp35[6:TempSeqLen + 6] and "-" not in dataset_df['offtarget_sequence'][j]:
                            print(j)
                            print(dataset_df['offtarget_sequence'][j])
                            print(temp35[6:TempSeqLen + 6])
                            print(orig_str)
                            raise ValueError
                        dataset_df.loc[j,SixAfterName] = str(temp35[TempSeqLen + 6:TempSeqLen + 12])
                        dataset_df.loc[j,SixBeforeName] = str(temp35[0:6])
                        # dataset_df['All32'][j] = str(temp35)  # TODO put in comment

                if i.id[0:9] == "NC_000023":
                    for j in dataset_df.loc[dataset_df['chrom'] == 'chrX'].index:
                        if CheckSequence(dataset_df['offtarget_sequence'][j], j):
                            continue
                        TempSeqLen = len(dataset_df['offtarget_sequence'][j])
                        if dataset_df['strand'][j] == '+':
                            if dataset_df['offtarget_sequence'][j] != i.seq[ChromStartVector[j]:chromEndVector[j]].upper() and "-" not in dataset_df['offtarget_sequence'][j]:
                                print(j)
                                print(dataset_df['offtarget_sequence'][j])
                                print(i.seq[ChromStartVector[j]:chromEndVector[j]])
                                raise ValueError
                            temp35 = i.seq[ChromStartVector[j] - 6:chromEndVector[j] + 6]
                            dataset_df.loc[j,SixAfterName] = str(temp35[TempSeqLen + 6:TempSeqLen + 12]).upper()
                            dataset_df.loc[j,SixBeforeName] = str(temp35[0:6]).upper()
                            # dataset_df['All32'][j] = str(temp35).upper()
                        else:
                            orig_str = str(i.seq[ChromStartVector[j] - 6:chromEndVector[j] + 6]).upper()
                            temp35 = ReverseAndChange(orig_str)
                            if dataset_df['offtarget_sequence'][j] != temp35[6:TempSeqLen + 6] and "-" not in dataset_df['offtarget_sequence'][j]:
                                print(j)
                                print(dataset_df['offtarget_sequence'][j])
                                print(temp35[6:29])
                                print(orig_str)
                                raise ValueError
                            dataset_df.loc[j,SixAfterName] = str(temp35[TempSeqLen + 6:TempSeqLen + 12])
                            dataset_df.loc[j,SixBeforeName] = str(temp35[0:6])
                            # dataset_df['All32'][j] = str(temp35)
                if i.id[0:9] == "NC_000024":
                    for j in dataset_df.loc[dataset_df['chrom'] == 'chrY'].index:
                        if CheckSequence(dataset_df['offtarget_sequence'][j], j):
                            continue
                        TempSeqLen = len(dataset_df['offtarget_sequence'][j])
                        if dataset_df['strand'][j] == '+':
                            if dataset_df['offtarget_sequence'][j] != i.seq[ChromStartVector[j]:chromEndVector[j]].upper() and "-" not in dataset_df['offtarget_sequence'][j]:
                                print(j)
                                print(dataset_df['offtarget_sequence'][j])
                                print(i.seq[ChromStartVector[j]:chromEndVector[j]])
                                raise ValueError
                            temp35 = i.seq[ChromStartVector[j] - 6:chromEndVector[j] + 6]
                            dataset_df.loc[j,SixAfterName] = str(temp35[TempSeqLen + 6:TempSeqLen + 12]).upper()
                            dataset_df.loc[j,SixBeforeName] = str(temp35[0:6]).upper()
                            # dataset_df['All32'][j] = str(temp35).upper()
                        else:
                            orig_str = str(i.seq[ChromStartVector[j] - 6:chromEndVector[j] + 6]).upper()
                            temp35 = ReverseAndChange(orig_str)
                            if dataset_df['offtarget_sequence'][j] != temp35[6:TempSeqLen + 6] and "-" not in dataset_df['offtarget_sequence'][j]:
                                print(j)
                                print(dataset_df['offtarget_sequence'][j])
                                print(temp35[6:29])
                                print(orig_str)
                                raise ValueError
                            dataset_df.loc[j,SixAfterName] = str(temp35[TempSeqLen + 6:TempSeqLen + 12])
                            dataset_df.loc[j,SixBeforeName] = str(temp35[0:6])
                            # dataset_df['All32'][j] = str(temp35)
                dataset_df.to_excel(ResultExcelPath)


AddingSixBeforeAndSixAfterToExcel(GenomePath=general_utilities.HOME_DIR + "GRCh38_latest_genomic.fna", OrigExcelPath=general_utilities.CHANGE_SEQ_PATH,
                                  SixBeforeName='SixBefore', SixAfterName='SixAfter', ChromStartName='chromStart', ChromEndName='chromEnd',
                                  ResultExcelPath=general_utilities.DATASETS_PATH + "CHANGE-seqNew.xlsx")

AddingSixBeforeAndSixAfterToExcel(GenomePath=general_utilities.HOME_DIR + "GRCh38_latest_genomic.fna", OrigExcelPath=general_utilities.GUIDE_SEQ_PATH,
                                  SixBeforeName='SixBefore', SixAfterName='SixAfter', ChromStartName='chromStart', ChromEndName='chromEnd',
                                  ResultExcelPath=general_utilities.DATASETS_PATH + "GUIDE-seqNew.xlsx")
