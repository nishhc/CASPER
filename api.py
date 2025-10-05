from casper.utils.target_input_handler import *
from casper.set_generation.base_set_generator import *
from casper.set_generation.base_filters import *
from casper.set_generation.complex_features import *
from casper.utils import *
from casper.ranking.ranker import *
import shutil
import os
import json
from fastapi import FastAPI, UploadFile, Form, Request
from pydantic import BaseModel
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware 
from fastapi.responses import JSONResponse
import pandas as pd
import asyncio
from fastapi.responses import StreamingResponse

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
    request: Request,
    config: str = Form(...),  
    target_fasta: UploadFile = None,
    input_csv: UploadFile = None
):
    async def stream_generator():
        try:
            print("Received config:", config)
            print("Received FASTA:", target_fasta.filename if target_fasta else "None")
            print("Received CSV:", input_csv.filename if input_csv else "None")

            yield json.dumps({"status": f"Received config:, {config}"}) + "\n"
            yield json.dumps({"status": f"Received FASTA:, {target_fasta.filename if target_fasta else 'None'}"}) + "\n"
            yield json.dumps({"status": f"Received CSV:, {input_csv.filename if input_csv else 'None'}"}) + "\n"
            
            if await request.is_disconnected():
                raise asyncio.CancelledError

            config_model = PrimerConfig(**json.loads(config))
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

            yield json.dumps({"status": "Preprocessing sequence data..."}) + "\n"
            if await request.is_disconnected():
                raise asyncio.CancelledError
            seq = SequenceData(fasta_path)
            await asyncio.to_thread(seq.preprocess)
            
            GENERATION = csv_path is None and True
            if GENERATION:
                yield json.dumps({"status": "Generating primer pairs. This will take a while..."}) + "\n"
                if await request.is_disconnected():
                    raise asyncio.CancelledError
                gen = PairGenerator(seq)
                await asyncio.to_thread(
                    gen.generate_primers_pairs,
                    seq.sequence,
                    [config_model.min_primer_length, config_model.max_primer_length],
                    [config_model.min_amplicon_length, config_model.max_amplicon_length],
                    crrnalen=[config_model.min_crrna_length, config_model.max_crrna_length],
                )
                await asyncio.to_thread(gen.to_csv, f"{output_dir}/pairs.csv")
            
            yield json.dumps({"status": "Filtering primers..."}) + "\n"
            if await request.is_disconnected():
                raise asyncio.CancelledError
            ftr = PrimerFilter(f"{output_dir}/pairs.csv" if GENERATION else csv_path)
            filtered = await asyncio.to_thread(ftr.run, [config_model.min_gc_content, config_model.max_gc_content])
            await asyncio.to_thread(filtered.to_csv, f"{output_dir}/filtered_sequences.csv")

            yield json.dumps({"status": "Calculating features..."}) + "\n"
            if await request.is_disconnected():
                raise asyncio.CancelledError
            features = FeatureCalculator(f"{output_dir}/filtered_sequences.csv")
            await asyncio.to_thread(features.compute_all_features)
            await asyncio.to_thread(features.to_csv, f"{output_dir}/output_with_features.csv")

            yield json.dumps({"status": "Ranking generated sets..."}) + "\n"
            if await request.is_disconnected():
                raise asyncio.CancelledError
            ranker = Ranker(f"{output_dir}/output_with_features.csv")
            await asyncio.to_thread(ranker.rank)
            print("Finished ranking")
            await asyncio.to_thread(ranker.to_csv, f"{output_dir}/ranked.csv", config_model.num_sets)

            ranked_csv_path = f"{output_dir}/ranked.csv"
            df = pd.read_csv(ranked_csv_path)
            json_data = df.to_dict(orient='records')
            with open(ranked_csv_path, 'r') as f:
                csv_string = f.read()
                
            final_payload = {
                "jsonData": json_data,
                "csvData": csv_string
            }
            yield json.dumps({"status": "complete", "data": final_payload}) + "\n"

        except asyncio.CancelledError:
            print("Client disconnected, job cancelled.")
        finally:
            print("Stream finished.")

    return StreamingResponse(stream_generator(), media_type="application/x-ndjson")