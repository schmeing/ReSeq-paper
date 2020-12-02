#!/usr/bin/env python3

import csv
import getopt
import numpy as np
from scipy.stats import spearmanr
import sys

def find_id_of_minimum(y, interval=15):
    min_id = 1
    while min_id+1 < len(y) and y[min_id] < y[min_id + 1]: # Find the maximum of the systematic error peak (which might be already at 1)
        min_id += 1
    #print(min_id)
    while min_id+1 < len(y) and y[min_id] > np.min(y[min_id+1:min_id+2+interval]): # Find the minimum splitting systematic errors and signal peak
        #print(y[min_id], np.min(y[min_id+1:min_id+2+interval]))
        min_id = min_id+1+np.argmin(y[min_id+1:min_id+2+interval])
        #print(min_id, y[min_id])

    return min_id

def score_kmer(kmerFiles, title = ""):
    names = []
    for sf in kmerFiles:
        if 'gz' == sf.split('.')[-1]:
            names.append( sf.split('.',2)[0] )
        else:
            names.append( sf.rsplit('.',1)[0] )

    for i, name in enumerate(names):
        if '/' in name:
            names[i] = name.rsplit('/',1)[1].split('.')[0][0:15]
    
    names = names[1:]
    
    specs = []
    for sf in kmerFiles:
        with open(sf) as sFile:
            sp = csv.reader(sFile, delimiter=' ', skipinitialspace=True)
            x = []
            y = []
            for row in sp:
                if len(row) == 2:
                    x.append(int(row[0]))
                    y.append(int(row[1]))
                    
            if 0 == len(x):
                sFile.seek(0)
                sp = csv.reader(sFile, delimiter='\t', skipinitialspace=True)
                for row in sp:
                    if len(row) == 2:
                        x.append(int(row[0]))
                        y.append(int(row[-1]))
            
            specs.append( np.zeros(max(x)+1) )
            specs[-1][x] = y
    
    threshold = 15
    if len(title):
        print(title)
        
        if title == "Ec-Hi2000-TruSeq":
            threshold = 400
        elif title == "Ec-Hi2500-TruSeq" or "Ec-Hi2500-TruSeq-asm":
            threshold = 800
        elif title == "Ec-Hi4000-Nextera":
            threshold = 600
        elif title == "Bc-Hi4000-Nextera":
            threshold = 2000
        elif title == "Rs-Hi4000-Nextera":
            threshold = 3000
        elif title == "At-HiX-TruSeq":
            threshold = 50
        elif title == "Mm-HiX-Unknown":
            threshold = 60
        elif title == "Hs-HiX-TruSeq":
            threshold = 60
        elif title == "Hs-Nova-TruSeq":
            threshold = 150
        elif title == "Ec-Mi-TruSeq":
            threshold = 300
        elif title == "At-BGI":
            threshold = 200
        
        title = title + " "

    if 1 < len(specs):
        metric = title + "Kmer"
        
        max_len = max([ len(sp) for sp in specs ])
        for i in range(len(specs)):
            specs[i].resize(max_len)
        real = specs[0]
        sims = specs[1:]
        x_minimum = find_id_of_minimum(real, threshold//50)
        
        scores = []
        for name, sim in zip(names,sims):
            trunc_real = real[1:x_minimum+1]
            trunc_sim = sim[1:x_minimum+1]
            with np.errstate(divide='ignore', invalid='ignore'):
                    rel = np.where((0 == trunc_real) | (0 == trunc_sim), trunc_sim== trunc_real, np.where(trunc_sim > trunc_real, trunc_real/trunc_sim, trunc_sim/trunc_real))
            mean_rel = np.mean(rel)
            cor = spearmanr(trunc_real, trunc_sim)[0]
            abs_val = np.sum(np.abs(sim[x_minimum+1:]/np.sum(sim[x_minimum+1:])-real[x_minimum+1:]/np.sum(real[x_minimum+1:])))

            print( metric, name, mean_rel, cor, abs_val )
            
            if (cor > 0.98) and (mean_rel > 0.75) and (abs_val < 0.2):
                scores.append(1)
            elif (cor > 0.9) and (mean_rel > 0.6) and (abs_val < 0.37):
                scores.append(2)
            elif (cor > 0.75) and (mean_rel > 0.3) and (abs_val < 0.52):
                scores.append(3)
            elif (cor > 0.48) and (mean_rel > 0.05) and (abs_val < 0.75):
                scores.append(4)
            else:
                scores.append(5)
            
        print( metric, list(zip(names,scores)))


def usage():
    print( "Usage: python score_kmer.py [OPTIONS] File File2 [File3 ...]" )
    print( "Scores kmer spectrum of File2 [File3 ...] compared to File." )
    print( "  -h, --help            display this help and exit" )
    print( "  -t, --title           title added before each score" )
    pass

def main(argv):
    try:
        optlist, args = getopt.getopt(argv, 'ht:', ['help','title='])
    except getopt.GetoptError:
        print( "Unknown option\n" )
        usage()
        sys.exit(2)

    title = ""
    for opt, par in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-t", "--title"):
            title = par
            
    if 2 > len(args):
        print( "Minimum two files are required.\n" )
        usage()
        sys.exit(2)
        pass

    score_kmer(args,title)
    pass

if __name__ == "__main__":
    main(sys.argv[1:])
    pass

