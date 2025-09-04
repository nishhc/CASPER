import sys

class SequenceData:
    def __init__(self, filepath):
        self.filepath = filepath
        self.sequence_name = None
        self.sequence = ""
        self.length = 0
        self.gc_content = 0.0
        self.max_homopolymer_run = 0

    def preprocess(self):
        try:
            with open(self.filepath, 'r') as f:
                lines = f.readlines()
                if not lines:
                    raise ValueError("Input file is empty.")
                
                if lines[0].startswith('>'):
                    self.sequence_name = lines[0].strip()[1:]
                else:
                    raise ValueError("FASTA file must start with a header line (>).")
                
                raw_sequence = "".join(line.strip() for line in lines[1:])
                self.sequence = self.cleanseq(raw_sequence)
                
            self.stat()
            
            print("[SYSTEM] Preprocessing complete.")
            print(f"Sequence Name: {self.sequence_name}")
            print(f"Sequence: {self.sequence}")
            print(f"Sequence Length: {self.length}")
            print(f"GC Content: {self.gc_content:.2f}%")
            print(f"Max Homopolymer Run: {self.max_homopolymer_run}")

        except FileNotFoundError:
            sys.exit(f"Error: The file at '{self.filepath}' was not found.")
        except ValueError as e:
            sys.exit(f"Error: {e}")
            
    def cleanseq(self, seq):
        seq_upper = seq.upper()
        valid_bases = {'A', 'T', 'C', 'G', 'U', 'N'}
        if not all(base in valid_bases for base in seq_upper):
            invalid_chars = set(c for c in seq_upper if c not in valid_bases)
            raise ValueError(f"Invalid characters found in sequence: {invalid_chars}")
        
        return seq_upper.replace('U', 'T').replace('N', '')
        
    def stat(self):
        self.length = len(self.sequence)
        
        if self.length > 0:
            gc_count = self.sequence.count('G') + self.sequence.count('C')
            self.gc_content = (gc_count / self.length) * 100
            
            max_run = 0
            for base in ['A', 'T', 'C', 'G']:
                run_length = 0
                for i in range(self.length):
                    if self.sequence[i] == base:
                        run_length += 1
                    else:
                        if run_length > max_run:
                            max_run = run_length
                        run_length = 0
                if run_length > max_run:
                    max_run = run_length
            self.max_homopolymer_run = max_run
            