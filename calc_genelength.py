import os, sys
import numpy as np

# Calculate gene length by exons
# Useful for TPM

class ensembl:
    def __init__(self, name):
        self.name = name
        self.exons = []
        self.exonpos = []

    def add_exon(self, ex, x1pos, x2pos):
        self.exonpos.append((x1pos, x2pos))
        self.exons.append(ex)
    

if __name__ == '__main__':

    # read gtf file
    gtffile = sys.argv[1]
    dGenes = dict()
    
    with open(gtffile, 'r') as f:
        for line in f:
            if '##' not in line:
               # print(line)
                if 'gene' == line.split()[2]:
                    genename = line.split('\"')[1]
                    obj = ensembl(genename)
                    dGenes[genename] = obj

                elif 'exon' == line.split()[2]:
                    genename = line.split('\"')[1]
                    exonname = str(line.split('exon_id')[1].split()[0].split('\"')[1])
                    dGenes[genename].add_exon(exonname, int(line.split()[3]), int(line.split()[4]))

                else:
                    pass

    for k, v in dGenes.items():
        lpos = [element for tupl in v.exonpos for element in tupl]
        imax = max(lpos)
        imin = min(lpos)    
        ilen = imax-imin
        arr = np.zeros(ilen)
        
        for x1, x2 in v.exonpos:
            x1_rel = x1-imin
            x2_rel = x2-imin
            arr[x1_rel:x2_rel] = 1

        print(k, np.count_nonzero(arr)+1)
