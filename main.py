from casper.utils.target_input_handler import *
from casper.set_generation.base_set_generator import *
from casper.set_generation.base_filters import *
from casper.set_generation.complex_features import *
from casper.utils import *
from casper.ranking.ranker import *
import shutil
import argparse
import os
import configparser

def main():
    config = configparser.ConfigParser()
    config.read("config.ini")
    run_headless = int(config["FEATURES"]["MIN_PRIMER_LENGTH"])
    if run_headless:
        headless(config)

def headless(config):
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
    output_dir = "output/"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    
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

    try:
        min_primer_length = int(config["FEATURES"]["MIN_PRIMER_LENGTH"])
        max_primer_length = int(config["FEATURES"]["MAX_PRIMER_LENGTH"])

        min_amplicon_length = int(config["FEATURES"]["MIN_AMPLICON_LENGTH"])
        max_amplicon_length = int(config["FEATURES"]["MAX_AMPLICON_LENGTH"])

        min_crrna_length = int(config["FEATURES"]["MIN_CRRNA_LENGTH"])
        max_crrna_length = int(config["FEATURES"]["MAX_CRRNA_LENGTH"])

        min_gc_content = int(config["FEATURES"]["MIN_GC_CONTENT"])
        max_gc_content = int(config["FEATURES"]["MAX_GC_CONTENT"])

        sets = int(config["FEATURES"]["NUM_SETS"])
    except:
        clog("File \"config.ini\" not found. Running with defaults.", "CONFIG PARSER")

    sequence_data = SequenceData(TARGET_FASTA)
    sequence_data.preprocess()
    
    if GENERATION:
        primer_generator = PairGenerator(sequence_data)
        primer_generator.generate_primers_pairs(sequence_data.sequence, [min_primer_length, max_primer_length], [min_amplicon_length, max_amplicon_length], crrnalen=[min_crrna_length, max_crrna_length])
        primer_generator.to_csv("pairs.csv")
        shutil.move("pairs.csv", "output/")
        print("Base sequences generated")

    primer_filter = PrimerFilter("output/pairs.csv" if GENERATION else OPTIONAL_INPUT_CSV)
    filtered = primer_filter.run([min_gc_content, max_gc_content])
    filtered.to_csv("filtered_sequences.csv") 
    shutil.move("filtered_sequences.csv", "output/")
    print("Filtered sequences saved to filtered_sequences.csv\n")

    primer_features = FeatureCalculator("output/filtered_sequences.csv")
    primer_features.compute_all_features()
    primer_features.to_csv("output_with_features.csv")
    shutil.move("output_with_features.csv", "output/")
    
    ranking = Ranker("output/output_with_features.csv")
    ranking.rank()
    print("Finished ranking")
    ranking.to_csv("ranked.csv", sets)
    shutil.move("ranked.csv", "output/")

if __name__ == "__main__":
    main()