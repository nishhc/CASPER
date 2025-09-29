from casper.utils.target_input_handler import *
from casper.set_generation.base_set_generator import *
from casper.set_generation.base_filters import *
from casper.set_generation.complex_features import *
from casper.ranking.ranker import *
import shutil
import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A script to generate, filter, and rank primer pairs for a given target sequence."
    )

    parser.add_argument("--target-fasta", 
                        required=True,
                        help="Path to the target FASTA file."
                        )

    parser.add_argument(
        "--input-csv",
        help="Path to an existing primer CSV file. If not provided, primers will be generated."
    )

    args = parser.parse_args()
    output_dir = "output"

    if os.path.exists(output_dir):
        for filename in os.listdir(output_dir):
            file_path = os.path.join(output_dir, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    TARGET_FASTA = args.target_fasta
    OPTIONAL_INPUT_CSV = args.input_csv
    GENERATION = OPTIONAL_INPUT_CSV is None 

    sequence_data = SequenceData(TARGET_FASTA)
    sequence_data.preprocess()

    if GENERATION:
        primer_generator = PairGenerator(sequence_data)
        primer_generator.generate_primers_pairs(sequence_data.sequence)
        primer_generator.to_csv("pairs.csv")
        shutil.move("pairs.csv", "output/")
        print("Base sequences generated")

    primer_filter = PrimerFilter("output/pairs.csv" if GENERATION else OPTIONAL_INPUT_CSV)
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
