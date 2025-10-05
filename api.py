from casper.utils.target_input_handler import *
from casper.set_generation.base_set_generator import *
from casper.set_generation.base_filters import *
from casper.set_generation.complex_features import *
from casper.utils import *
from casper.ranking.ranker import *
import shutil
import os
import json
from fastapi import FastAPI, UploadFile, Form
from pydantic import BaseModel
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware 
from fastapi.responses import JSONResponse
import pandas as pd

app = FastAPI()

origins = [
    "http://localhost:3000",
    "http://localhost:5173",
    "http://localhost:8080",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class PrimerConfig(BaseModel):
    min_primer_length: int
    max_primer_length: int
    min_amplicon_length: int
    max_amplicon_length: int
    min_crrna_length: int
    max_crrna_length: int
    min_gc_content: int
    max_gc_content: int
    num_sets: int
    generation: bool

@app.post("/process")
async def process_files(
    config: str = Form(...),  
    target_fasta: UploadFile = None,
    input_csv: UploadFile = None
):
    print("Received config:", config)
    print("Received FASTA:", target_fasta.filename if target_fasta else "None")
    print("Received CSV:", input_csv.filename if input_csv else "None")
    config = PrimerConfig(**json.loads(config))

    output_dir = "output/"
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

    fasta_path = None
    if target_fasta:
        fasta_path = os.path.join(output_dir, "input.fasta")
        with open(fasta_path, "wb") as f:
            f.write(await target_fasta.read())

    csv_path = None
    if input_csv:
        csv_path = os.path.join(output_dir, "input.csv")
        with open(csv_path, "wb") as f:
            f.write(await input_csv.read())

    seq = SequenceData(fasta_path)
    seq.preprocess()

    GENERATION = csv_path is None and True

    if GENERATION:
        gen = PairGenerator(seq)
        gen.generate_primers_pairs(
            seq.sequence,
            [config.min_primer_length, config.max_primer_length],
            [config.min_amplicon_length, config.max_amplicon_length],
            crrnalen=[config.min_crrna_length, config.max_crrna_length],
        )
        gen.to_csv(f"{output_dir}/pairs.csv")
    

    ftr = PrimerFilter(f"{output_dir}/pairs.csv" if GENERATION else csv_path)
    filtered = ftr.run([config.min_gc_content, config.max_gc_content])
    filtered.to_csv(f"{output_dir}/filtered_sequences.csv")

    features = FeatureCalculator(f"{output_dir}/filtered_sequences.csv")
    features.compute_all_features()
    features.to_csv(f"{output_dir}/output_with_features.csv")

    ranker = Ranker(f"{output_dir}/output_with_features.csv")
    ranker.rank()
    print("Finished ranking")
    ranker.to_csv(f"{output_dir}/ranked.csv", config.num_sets)

    ranked_csv_path = f"{output_dir}/ranked.csv"

    df = pd.read_csv(ranked_csv_path)
    json_data = df.to_dict(orient='records')

    with open(ranked_csv_path, 'r') as f:
        csv_string = f.read()

    return JSONResponse(content={
        "jsonData": json_data,
        "csvData": csv_string
    })