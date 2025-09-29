# CASPER

**CASPER (Combined Amplification & Spacer Engine for RPA-Cas12a)**  
A combined RPA primer and CRISPR-Cas12a crRNA designer optimized to create pairs that work together effectively.  

Given a target (and optional background genomes), CASPER generates sequences and produces rankings for candidate RPA amplicons, primers, and crRNAs using a **composite ranking framework**.  
This framework integrates multiple features important for RPA + Cas12a design and applies **weights trained on 100 experimental primer–crRNA pairs** from previous research.

---

## Features

CASPER evaluates a wide range of features to rank candidate primer–crRNA pairs.  

### 1. GC Content
- **Forward primer GC%**: 30–70%  
- **Reverse primer GC%**: 30–70%  
- **Amplicon GC%**: 30–70%  
- Ensures stable but not overly strong binding.  

---

### 2. Primer 3′ End Stability
- `fp_3p_self_run` and `rp_3p_self_run` ≤ 2  
- Prevents self-priming and nonspecific elongation.  

---

### 3. Homopolymer Runs
- Forward/reverse primers: max run ≤ 4  
- Amplicon: max run ≤ 5  
- Avoids slippage or polymerase stalling.  

---

### 4. Dimerization Checks
- **Self-dimers**: `fp_self_dimer_pct`, `rp_self_dimer_pct`  
- **Cross-dimers**: `fp_rp_cross_dimer_pct`  
- Ensures efficient amplification.  

---

### 5. Thermodynamic Balance
- Primer melting temperature (**Tm**): `fp_tm_C`, `rp_tm_C`  
- ΔTm (`delta_tm_C`) kept within an acceptable range  
- Reduces imbalance in primer annealing.  

---

### 6. Protospacer/PAM Overlap
- Ensures primers do not overlap Cas12a protospacer or PAM (`overlap_protospacer`, `overlap_pam`)  
- Prevents competition between amplification and Cas12a targeting.  

---

### 7. crRNA Off-Target Risk
- `crrna_offtarget_mm_imp`: mismatch impact scores based on genome background scans  
- Penalizes guides likely to bind unintended loci.  

---

### 8. Primer Off-Target Risk
- Forward primer: `fp_offtarget_mm_imp`  
- Reverse primer: `rp_offtarget_mm_imp`  
- Ensures specific amplification.  

---

### 9. Guide Seed Region Stability
- Seed **GC%** (`seed_gc_pct`)  
- Seed **max run** (`seed_max_run`)  
- Seed **unpaired fraction** (`guide_seed_unpaired_frac`)  
- Ensures guide accessibility and prevents overly rigid binding.  

---

### 10. Guide RNA Folding
- Predicted minimum free energy: `guide_mfe_kcal`  
- Avoids crRNAs with highly stable hairpins that hinder Cas12a binding.  

---
## Scoring
CASPER’s scoring framework is built on weights derived from experimental RPA–CRISPR-Cas12a designs reported in prior research. Each feature is assigned an importance value based on how strongly it influenced successful amplification and detection in these studies.

Currently, CASPER uses these predetermined weights to evaluate candidate primer–crRNA pairs, ensuring that the scoring reflects patterns validated in published literature.
Our predetermined weights are as follows
```
    "crrna_offtarget_mm_imp": 0.1096,
    "fp_offtarget_mm_imp":    0.1005,
    "rp_offtarget_mm_imp":    0.1005,

    "fp_3p_self_run":         0.0593,
    "rp_3p_self_run":         0.0593,
    "fp_rp_3p_cross_run":     0.0593,

    "fp_self_dimer_pct":      0.0365,
    "rp_self_dimer_pct":      0.0365,
    "fp_rp_cross_dimer_pct":  0.0365,
    "fp_gc_pct":              0.0228,
    "rp_gc_pct":              0.0228,
    "amplicon_gc_pct":        0.0228,

    "seed_max_run":           0.0319,
    "seed_gc_pct":            0.0228,
    "guide_mfe_kcal":         0.0136,
    "guide_seed_unpaired_frac": 0.0274,

    "delta_tm_C":             0.0913,
    "fp_tm_C":                0.0274,
    "rp_tm_C":                0.0274,

    "overlap_protospacer":    0.0182,
    "overlap_pam":            0.0182,
```

We are also developing a feature that will allow users to customize or adjust weights based on their own experimental datasets—making CASPER adaptable to diverse targets, conditions, and research goals.
---
##  Implementation

```
git clone https://github.com/nishhc/CASPER.git
conda create --name CASPER python=3.12
conda activate CASPER
cd CASPER
pip install -r requirements.txt
```
To use the generator model to generate new RPA primers and crRNAs, add the fasta file of your target fasta into the main CASPER directory.
Then, run
```
python src/main.py
```
The script will run through multiple steps and end with all the ranked sequences with the scores in ranked.csv

To use the predictor model to predict the score of existing RPA primers and crRNAs, add the sequences to target.fasta. Then, run
```
python src/orderedsequences.csv
```
The resulting sequences will be in rankedsequences.csv.
---
