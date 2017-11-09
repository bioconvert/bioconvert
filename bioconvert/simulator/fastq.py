



class FastqSim():
    def __init__(self, outfile):
        self.outfile = outfile
        self.nreads = 1000000
        # do we want a wrap or unwrap version ?
        # For now, unwrap
        self.wrap = False
        self.read_length = 250

    def simulate(self):
        RL = self.read_length
        with open(self.outfile, "w") as fout:
            for i in range(self.nreads):
                fout.write(">identifier whatever it means but long enough\n")
                fout.write("ACGT" * (RL // 4) + "A" * (RL % 4) + "\n")
                fout.write("+\n")
                fout.write("C" * RL + "\n")
        
