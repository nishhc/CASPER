class Pair:
    def __init__(self, fp, bp, amp, crRNA = ""):
        self.fp = fp
        self.bp = bp
        self.amp = amp
        self.crRNA = crRNA

    def __str__(self):
        return f"{self.fp}, {self.bp}, {self.amp}, {self.crRNA}"
    
    