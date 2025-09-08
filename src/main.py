from utils.tgtinputhandler import *
from primergen import *
from filtering import *

TARGET_FASTA = "target.fasta"

if __name__ == '__main__':
    sequence_data = SequenceData(TARGET_FASTA)
    sequence_data.preprocess()

    primer_generator = PrimerGenerator(sequence_data)
    primer_generator.generate_primers()