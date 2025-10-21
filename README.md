# CASPER

**CASPER (Combined Amplification & Spacer Engine for RPA-Cas12a)**  
A combined RPA primer and CRISPR-Cas12a crRNA designer optimized to create pairs that work together effectively.  

Given a target sequence, CASPER generates sequences and produces rankings for candidate RPA amplicons, primers, and crRNAs using a **composite ranking framework**.  

This framework integrates multiple features important for RPA + Cas12a design and applies **weights trained on 100 experimental primer–crRNA pairs** from previous research.

---

# Features

Read 2025.igem.wiki/lambert-ga/software

##  Implementation

```
git clone https://github.com/nishhc/CASPER.git
pip install -r requirements.txt
pip install -e .
```
To use the generator model to generate new RPA primers and crRNAs, add the fasta file of your target fasta into the main CASPER directory.
Then, run
```
casper --target-fasta target.fasta
```

To use the ranking model to rank existing RPA primers and crRNAs, add the fast file of your target fasta into the main CASPER directory and a csv file called existing_primers.csv. Label each header as so 
| forward_primer                         | backward_primer                      | amplicon                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | crrna                 |
|----------------------------------------|--------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------|
| GGTTGGTCTGGTTGGCCCGTGTGTCATTAC         | ACAAACCCAACCCATCCCATCCCGCCACCA      | GGTTGGTCTGGTTGGCCCGTGTGTCATTACGGGTTGGATAAGATAGTAAGTGCAATCTGGGGGTTGTGTTTGAGCGGCGTTTCAGTTGTTTATTTCCCTTTGTTATTCCCTTTGGGGTTGTTGTTTGGTTGTGTGTTTATACCAGCTTATTCAATTCACTTGGTGGTGGTGGCGGGATGGGATGGGTTGGGTTTGT | AGCGGCGTTTCAGTTGTTTA |
| CGTGTGTCATTACGGGTTGGATAAGATAGTA        | ACCACCACCAAGTGAATTGAATAAGCTGGTA      | CGTGTGTCATTACGGGTTGGATAAGATAGTAAGTGCAATCTGGGGGTTGTGTTTGAGCGGCGTTTCAGTTGTTTATTTCCCTTTGTTATTCCCTTTGGGGTTGTTGTTTGGTTGTGTGTTTATACCAGCTTATTCAATTCACTTGGTGGTGGT | AGCGGCGTTTCAGTTGTTTA |
| CGTGTGTCATTACGGGTTGGATAAGATAGTA        | ATAAGCTGGTATAAACACACAACCAAACAAC      | CGTGTGTCATTACGGGTTGGATAAGATAGTAAGTGCAATCTGGGGGTTGTGTTTGAGCGGCGTTTCAGTTGTTTATTTCCCTTTGTTATTCCCTTTGGGGTTGTTGTTTGGTTGTGTGTTTATACCAGCTTAT | AGCGGCGTTTCAGTTGTTTA |
Then, run
```
casper --target-fasta target.fasta --input-csv existing_primers.csv
```

The script will run through multiple steps and end with all the ranked sequences with the scores in output/ranked.csv

The resulting sequences will be in ranked.csv.
