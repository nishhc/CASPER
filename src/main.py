from input_handling import *

TARGET_FASTA = "target.fasta"

if __name__ == '__main__':
    sequence_data = SequenceData(TARGET_FASTA)
    sequence_data.preprocess()