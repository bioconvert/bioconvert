
cdef char* name = NULL
cdef char* seq = NULL

cdef _fastq2fasta(infile, outfile):


    with open(outfile, "w") as fasta, open(infile, "r") as fastq:
        try:
            while True:
                name, seq, _, _ = next(fastq), next(fastq), next(fastq), next(fastq)
                fasta.write(">{}{}".format(name[1:], seq))
        except StopIteration:
            pass

    """with open(outfile, "w") as fasta, open(infile, "r") as fastq:
       for (name, seq, _) in readfq(fastq):
           fasta.write(">{}\n{}\n".format(name, seq))
    """

def fastq2fasta(infile, outfile):
    _fastq2fasta(infile, outfile)


cdef _fastq2fasta2(infile, outfile):


    with open(outfile, "w") as fasta, open(infile, "r") as fastq:
        try:
            while True:
                name2, seq2, _, _ = next(fastq), next(fastq), next(fastq), next(fastq)
                fasta.write(">{}{}".format(name2[1:], seq2))
        except StopIteration:
            pass

    """with open(outfile, "w") as fasta, open(infile, "r") as fastq:
       for (name, seq, _) in readfq(fastq):
           fasta.write(">{}\n{}\n".format(name, seq))
    """

def fastq2fasta2(infile, outfile):
    _fastq2fasta2(infile, outfile)


