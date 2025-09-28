from utils.tgtinputhandler import *
from pairgen import *
from filtering import *
from features import *
from ranking.rankingtwo import *
import csv


TARGET_FASTA = "target.fasta"
'''
sequence_data = SequenceData(TARGET_FASTA)
sequence_data.preprocess()

primer_generator = PairGenerator(sequence_data)
primer_generator.generate_primers_pairs()
print("Base sequences generated")
'''
primer_filter = PrimerFilter("ordered.csv")
filtered = primer_filter.run()
filtered.to_csv("filtered_sequences.csv")   
print("Filtered sequences saved to filtered_sequences.csv")

primer_features = FeatureCalculator("filtered_sequences.csv")
primer_features.compute_all_features()
primer_features.to_csv("output_with_features.csv")

ranking = RankerVibe("output_with_features.csv")
print("Finished ranking")
ranked = ranking.rank()
ranked.to_csv("rankedsequences.csv", index=False)