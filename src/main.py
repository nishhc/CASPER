from utils.tgtinputhandler import *
from pairgen import *
from filtering import *
from features import *
from ranking.ranking import *
import csv


TARGET_FASTA = "target.fasta"

sequence_data = SequenceData(TARGET_FASTA)
sequence_data.preprocess()
'''
primer_generator = PairGenerator(sequence_data)
primer_generator.generate_primers_pairs()
print("Base sequences generated")

primer_filter = PrimerFilter("pairs.csv")
filtered = primer_filter.run()
filtered.to_csv("filtered_sequences.csv")   
print("Filtered sequences saved to filtered_sequences.csv")


primer_features = FeatureCalculator("filtered_sequences.csv")
primer_features.compute_all_features()
primer_features.to_csv("output_with_features.csv")

'''
ranking = Ranker("output_with_features.csv")
ranked = ranking.rank()
ranked.to_csv("ranked_output.csv")
