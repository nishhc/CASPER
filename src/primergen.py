from utils.tgtinputhandler import *

class PrimerGenerator:

    def __init__(self, target_data):
        self.target_data = target_data
        self.sequence = target_data.sequence

    def generate_primers(self):
        primers = []
        seq = self.sequence
        for i in range(1,len(seq)):
            primers.append(seq[i-1:i+1])
        