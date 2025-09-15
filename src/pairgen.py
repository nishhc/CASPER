from utils import *
from utils.seqhandler import *
from Pair import *
import csv

class PairGenerator:

    def __init__(self, target_data):
        self.target_data = target_data
        self.sequence = target_data.sequence
        self.primersPairs = []

    # forward primer - bp
    # backward primer - complementary
    def generate_primers_pairs(self):
        print("Generating primer pairs...")
        seen = set()
        primerLenRange = [28,36]
        amplen = [100,200]
        crrnalen = [20,24]
        seq = self.sequence
        # len - 1 (?) - 2(smallest primer) - smallest amplicion
        lastPotentialPos = len(seq)-(2*primerLenRange[0])-amplen[0]
        for curPos in range(0, lastPotentialPos+1):
            fp = ""
            bp = ""
            amplicon = ""
            for i in range(primerLenRange[0], primerLenRange[1]+1):
                fp = seq[curPos:curPos + i]
                for ampEnd in range(amplen[0], amplen[1]+1):
                    amplicon = seq[curPos:curPos+i+ampEnd]
                    has_pam = False
                    pam_locs = []
                    for pam in ["TTTA", "TTTC", "TTTG", "TTTT"]:
                        editAmp = amplicon
                        offset = 0
                        amp = editAmp.find(pam)
                        while amp != -1:
                            has_pam = True
                            pam_locs.append(offset + amp+4)
                            offset += amp + 1
                            editAmp = amplicon[offset:]
                            amp = editAmp.find(pam)
                        
                    if (not has_pam):
                        break
                    for x in range(primerLenRange[0], primerLenRange[1]+1):
                        bp = amplicon[len(amplicon)-1-x:]
                        if len(bp) >= primerLenRange[0]:
                            for i in pam_locs:
                                for cr_len in range(crrnalen[0], crrnalen[1]+1):
                                    crRNA = find_complementary(self, amplicon[i:i+cr_len])[::-1]
                                    if len(crRNA) >= crrnalen[0]:
                                        pair = Pair(fp, bp, amplicon, crRNA)
                                        key = (fp, find_complementary(self, bp)[::-1], amplicon, crRNA)
                                        if key not in seen:
                                            seen.add(key)
                                            self.primersPairs.append(Pair(*key))

        pairs = self.primersPairs  

        with open('pairs.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['forward_primer', 'backward_primer', 'amplicon', 'crrna'])
            for pair in pairs:
                writer.writerow([pair.fp, pair.bp, pair.amp, pair.crRNA])