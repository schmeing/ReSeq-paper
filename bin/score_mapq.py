#!/usr/bin/env python3

import csv
import getopt
import numpy as np
from scipy.stats import spearmanr
import sys

def score_mapq(mapqFiles, title = ""):
    names = []
    for sf in mapqFiles:
        if 'gz' == sf.split('.')[-1]:
            names.append( sf.split('.',2)[0] )
        else:
            names.append( sf.rsplit('.',1)[0] )

    for i, name in enumerate(names):
        if '/' in name:
            names[i] = name.rsplit('/',1)[1].split('.')[0].split('_')[3]
    
    names = names[1:]

    quals = []
    for mf in mapqFiles:
        with open(mf) as sFile:
            mq = csv.reader(sFile, delimiter=',', skipinitialspace=True)
            mapq = []
            count = []
            unmapped = 0
            next(mq, None)  # skip the headers
            for row in mq:
                if len(row) == 3:
                    if int(row[1]) == 0:
                        mapq.append(int(row[0]))
                        count.append(int(row[2]))
                    else:
                        unmapped = int(row[2])
            
            quals.append( np.zeros(max(mapq)+1) )
            quals[-1][mapq] = count
            quals[-1] = quals[-1] / ( np.sum(quals[-1]) + unmapped )
    
    if len(title):
        print(title)
        title = title + " "

    if 1 < len(quals):
        metric = title + "MapQ"
        
        real = quals[0]
        sims = quals[1:]
        
        scores = []
        for name, sim in zip(names,sims):
            abs_val = np.sum(np.abs(sim-real))

            #print( metric, name, abs_val )
            
            if abs_val < 0.05:
                scores.append(1)
            elif abs_val < 0.1:
                scores.append(2)
            elif abs_val < 0.25:
                scores.append(3)
            elif abs_val < 0.4:
                scores.append(4)
            else:
                scores.append(5)
            
        print( metric, list(zip(names,scores)))


def usage():
    print( "Usage: python score_mapq.py [OPTIONS] File File2 [File3 ...]" )
    print( "Scores mapping quality distribution in File2 [File3 ...] compared to File." )
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

    score_mapq(args,title)
    pass

if __name__ == "__main__":
    main(sys.argv[1:])
    pass

