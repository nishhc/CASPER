from casper.utils.target_input_handler import *
from casper.set_generation.base_set_generator import *
from casper.set_generation.base_filters import *
from casper.set_generation.complex_features import *
from casper.ranking.ranker import *
import shutil

TARGET_FASTA = "target.fasta"

sequence_data = SequenceData(TARGET_FASTA)
sequence_data.preprocess()

'''
primer_generator = PairGenerator(sequence_data)
primer_generator.generate_primers_pairs(sequence_data.sequence)
primer_generator.to_csv("pairs.csv")
shutil.move("pairs.csv", "output/")
print("Base sequences generated")
'''

primer_filter = PrimerFilter("output/pairs.csv")
filtered = primer_filter.run()
filtered.to_csv("filtered_sequences.csv") 
shutil.move("filtered_sequences.csv", "output/")
print("Filtered sequences saved to filtered_sequences.csv")

primer_features = FeatureCalculator("output/filtered_sequences.csv")
primer_features.compute_all_features()
primer_features.to_csv("output_with_features.csv")
shutil.move("output_with_features.csv", "output/")

ranking = Ranker("output/output_with_features.csv")
print("Finished ranking")
ranked = ranking.rank()
ranked.to_csv("ranked.csv", index=False)
shutil.move("ranked.csv", "output/")
