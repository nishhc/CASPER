def find_complementary(self, primer):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}        
        return ''.join(complement[base] for base in primer)