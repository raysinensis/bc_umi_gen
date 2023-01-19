import argparse
import random
import gzip
from itertools import groupby

# read names from our sequencer
def sim_name(ran, r):
    return("@A00217:536:H3GHKDSX5:4:1101:"+ str(ran) + ":1031 " + str(r) + ":N:0:NCCGTTCTCG+NCAATCCGTC")
# get gene
def read_fasta(fa):
    with open(fa) as f:
        faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
        for header in faiter:
            headerStr = header.__next__()[1:].strip()
            seq = "".join(s.strip() for s in faiter.__next__())
            yield headerStr, seq
# get barcode
def read_bc(file_bc):
    if file_bc:
        with open(file_bc, "r") as f:
            bc = [line.rstrip() for line in f]
    return(bc)
# random umi
def generate_umi(length):
    s = []
    for x in range(length):
        s += random.sample(["A", "G", "C", "T"], 1)
    return(''.join(x for x in s))
# random piece of gene
def generate_mapped(length, fasta):
    seq = random.sample(fasta, 1)[0][1]
    idx = random.randrange(0, len(seq) - length + 1)
    return(seq[idx : (idx + length)])
# set up seqs for r1 and r2, gene portion is more complicated
def sim_seq(r, length, fasta, bc1, bc2, bc3, numi, tn, sp1, sp2):
    if r == 1:
        seq = random.sample(bc1, 1)[0] + sp1 + random.sample(bc2, 1)[0] + sp2 + random.sample(bc3, 1)[0] + generate_umi(numi) + "T"*tn
        if length > len(seq):
            seq = seq + random.sample(["A", "G", "C"], 1)[0] + generate_umi(length - len(seq) - 1)
        else:
            seq = seq[0:length]
        return(seq)
    else:
        if len(fasta) == 0:
            seq = generate_umi(length)
        else:
            seq = generate_mapped(length, fasta)
        return(seq)
# 3rd line in fastq is just "+"
line3_template = "+"
# random, high quality scores
def sim_score(length):
    #d = {"F" : 19, ":" : 1}
    #v = [[k]*v for k,v in d.items()]
    #v = [item for sublist in v for item in sublist]
    #s = []
    #for x in range(length):
    #    s += random.sample(v, 1)
    #return(''.join(x for x in s))
    return("F"*length)
# putting one entry together
def sim_fastq(n, r, length, fasta, bc1, bc2, bc3, numi ,tn, sp1, sp2):
    return("\n".join([sim_name(n, r), sim_seq(r, length, fasta, bc1, bc2, bc3, numi, tn, sp1, sp2), line3_template, sim_score(length)]))
# writing two files for r1 and 2
def sim_fastqs(file_name, n, length1, length2, fasta, bc1, bc2, bc3, numi, tn, sp1, sp2):
    with gzip.open(file_name + "_R1.fastq.gz", 'wt') as f1, gzip.open(file_name + "_R2.fastq.gz", 'wt') as f2:
        if n > 1000:
            ticks = [int(n/100) * x for x in list(range(100))]
            k = 0
        else:
            ticks = []
        for i in range(n):
            seq2 = sim_fastq(i + 1, 2, length2, fasta, bc1, bc2, bc3, numi, tn, sp1, sp2)
            f2.write(seq2)
            f2.write("\n")
            seq1 = sim_fastq(i + 1, 1, length1, fasta, bc1, bc2, bc3, numi, tn, sp1, sp2)
            f1.write(seq1)
            f1.write("\n")
            if i in ticks:
                print(str(k) + "%...")
                k += 1

""" Generate synthetic fastqs with custom barcode and UMI structures
"""

def main():
    parser = argparse.ArgumentParser(description = """
        Generate synthetic fastqs with custom barcode and UMI structures
        """)

    parser.add_argument('-n',
                        '--totalreads',
                        help = """
                        total number of read pairs to generate
                        """,
                        required = True)

    parser.add_argument('-o',
                        '--output',
                        help ="""
                        output fastq name prefix
                        """,
                        default = "out",
                        required = False)

    parser.add_argument('-f',
                        '--fasta',
                        help ="""
                        fasta file of target genes mapped, will generate random sequence if empty
                        """,
                        required = True)
    
    parser.add_argument('-bc1',
                        '--barcode1',
                        help ="""
                        txt file containing 1st whitelist barcodes
                        """,
                        required = True)

    parser.add_argument('-bc2',
                        '--barcode2',
                        help ="""
                        txt file containing 2nd whitelist barcodes
                        """,
                        required = True)
    
    parser.add_argument('-bc3',
                        '--barcode3',
                        help ="""
                        txt file containing 3rd whitelist barcodes
                        """,
                        required = True)
    
    parser.add_argument('-sp1',
                        '--spacer1',
                        help ="""
                        txt file containing 1st spacer sequence
                        """,
                        required = True)
        
    parser.add_argument('-sp2',
                        '--spacer2',
                        help ="""
                        txt file containing 2nd spacer sequence
                        """,
                        required = True)
    
    parser.add_argument('-u',
                        '--umilength',
                        help ="""
                        UMI length
                        """,
                        default = 12,
                        required = False)
        
    parser.add_argument('-t',
                        '--tlength',
                        help ="""
                        polyT length
                        """,
                        default = 30,
                        required = False)
    
    parser.add_argument('-l1',
                        '--length1',
                        help ="""
                        read1 length (default 150)
                        """,
                        default = "150",
                        required = False)
    
    parser.add_argument('-l2',
                        '--length2',
                        help ="""
                        read2 length (default 150)
                        """,
                        default = "150",
                        required = False)

    # read in arguments
    args = parser.parse_args()
    n = int(args.totalreads)
    output = args.output
    fa = args.fasta
    bc1 = args.barcode1
    bc2 = args.barcode2
    bc3 = args.barcode3
    sp1 = args.spacer1
    sp2 = args.spacer2
    numi = int(args.umilength)
    l1 = int(args.length1)
    l2 = int(args.length2)
    tn = int(args.tlength)
    # proc from arugments
    fasta = list(read_fasta(fa))
    bc1 = read_bc(bc1)
    bc2 = read_bc(bc2)
    bc3 = read_bc(bc3)
    sp1 = read_bc(sp1)[0]
    sp2 = read_bc(sp2)[0]
    
    print("generating " + str(n) + " fastq.gz read pairs...")
    sim_fastqs(output, n, l1, l2, fasta, bc1, bc2, bc3, numi, tn, sp1, sp2)
    print("done")
    
if __name__ == '__main__': main()