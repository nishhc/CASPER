

import pandas as pd
import re

class PrimerFilter:
    def __init__(self, csv_filename):
        self.csv_filename = csv_filename

    @staticmethod
    def pam_cas12a_match(pam):
        return bool(re.match(r"TTT[ACG]$", pam))

    @staticmethod
    def gc_clamp_ok(seq):
        clamp = seq[-5:]
        gc_count = clamp.count('G') + clamp.count('C')
        return 1 <= gc_count <= 3

    @staticmethod
    def calc_gc_pct(seq):
        seq = seq.upper()
        gc = seq.count('G') + seq.count('C')
        return 100 * gc / len(seq) if len(seq) > 0 else 0

    @staticmethod
    def max_hpoly_run(seq):
        return max([len(m.group(0)) for m in re.finditer(r'(A+|T+|G+|C+)', seq.upper())] or [0])

    @staticmethod
    def max_run(seq):
        return max([len(m.group(0)) for m in re.finditer(r'(.)\1*', seq.upper())] or [0])

    @staticmethod
    def self_3p_run(seq):
        seq = seq.upper()
        last_base = seq[-1]
        run = 0
        for base in reversed(seq):
            if base == last_base:
                run += 1
            else:
                break
        return run

    
    def to_csv(self, output_csv: str):
        self.df.to_csv(output_csv, index=False)
        return output_csv
    def process_row(self, row):
        fp = row['forward_primer']
        rp = row['backward_primer']
        amplicon = row['amplicon']
        cr = row.get('crrna', '')
        guide_pam = cr[:4] if len(cr) >= 4 else ''
        overlap_pam = 0
        overlap_protospacer = 0
        return pd.Series({
            'fp_gc_pct': self.calc_gc_pct(fp),
            'rp_gc_pct': self.calc_gc_pct(rp),
            'amplicon_gc_pct': self.calc_gc_pct(amplicon),
            'fp_len': len(fp),
            'rp_len': len(rp),
            'amplicon_len': len(amplicon),
            'fp_3p_self_run': self.self_3p_run(fp),
            'rp_3p_self_run': self.self_3p_run(rp),
            'fp_max_hpoly_run': self.max_hpoly_run(fp),
            'rp_max_hpoly_run': self.max_hpoly_run(rp),
            'amplicon_max_run': self.max_run(amplicon),
            'overlap_pam': overlap_pam,
            'overlap_protospacer': overlap_protospacer
        })

    def run(self):
        df = pd.read_csv(self.csv_filename)
        features = df.apply(self.process_row, axis=1)
        data = pd.concat([df, features], axis=1)
        return data [
            (data['overlap_pam'] == 0) &
            (data['overlap_protospacer'] == 0) &   
            (data['fp_gc_pct'].between(35, 65)) &
            (data['rp_gc_pct'].between(35, 65)) &
            (data['amplicon_gc_pct'].between(35, 65)) &
            (data['fp_3p_self_run'] <= 2) &
            (data['rp_3p_self_run'] <= 2) &
            (data['fp_max_hpoly_run'] <= 4) &
            (data['rp_max_hpoly_run'] <= 4) &
            (data['amplicon_max_run'] <= 5)
        ]



