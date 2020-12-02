#!/usr/bin/env python3

import getopt
import numpy as np
import os
from scipy.stats import spearmanr
import sys

if 'RESEQ_PYMODS' in os.environ and os.path.isdir(os.path.realpath(os.environ['RESEQ_PYMODS'])):
    sys.path.append(os.path.realpath(os.environ['RESEQ_PYMODS']))
import DataStats

def convert_to_numpy(stats):
    max_len = 0
    for st in stats:
        max_len = max(max_len, st[0] + len(st[1]))
        
    npstats = []
    for st in stats:
        npst = np.zeros(max_len)
        npst[st[0]:(st[0] + len(st[1]))] = st[1]
        npstats.append(npst)
    return npstats[0], npstats[1:]

def score_reseq(statsFiles, title = ""):
    names = []
    for sf in statsFiles:
        if 'gz' == sf.split('.')[-1]:
            names.append( sf.split('.',2)[0] )
        else:
            names.append( sf.rsplit('.',1)[0] )

    for i, name in enumerate(names):
        if '/' in name:
            names[i] = name.rsplit('/',1)[1].split('.')[0][0:15]
    
    names = names[1:]
    
    stats = [];
    for sf in statsFiles:
        st = DataStats.DataStatsInterface(None)
        if st.Load(sf):
            stats.append(st)
        else:
            print( "Error loading statistics from {0}".format(sf) )
    
    if len(title):
        print(title)
        title = title + " "

    if 1 < len(stats):
        # Systematic errors
        metric = title + "SysError"
        real, sims = convert_to_numpy([st.ErrorCoverage() for st in stats])
        
        max_x = len(real)
        for x in range(len(real)):
            if real[x] < 10:
                max_x = x
                break

        real_cov, dummy = convert_to_numpy([st.Coverage() for st in stats[:1]])
        stop = 0.99*np.sum(real_cov)
        tot = 0
        max_x_cov = 0
        while tot < stop:
            tot += real_cov[max_x_cov]
            max_x_cov += 1
        max_x = min(max_x, max_x_cov)
 
        real = real[:max_x]/np.sum(real)
        
        scores = []
        for sim in sims:
            sim = sim[:max_x]/np.sum(sim)
            with np.errstate(divide='ignore', invalid='ignore'):
                rel = np.where((0 == real) | (0 == sim), sim == real, np.where(sim > real, real/sim, sim/real))
            
            mean_rel = np.mean(rel)
            cor = spearmanr(real, sim)[0]
            abs_val = np.sum(np.abs(real-sim))
            
            if (cor > 0.99) and (mean_rel > 0.6) and (abs_val < 0.25):
                scores.append(1)
            elif (cor > 0.97) and (mean_rel > 0.4) and (abs_val < 0.5):
                scores.append(2)
            elif (cor > 0.95) and (mean_rel > 0.2) and (abs_val < 0.7):
                scores.append(3)
            elif (cor > 0.85) and (mean_rel > 0.1) and (abs_val < 0.85):
                scores.append(4)
            else:
                scores.append(5)
            
        print( metric, list(zip(names,scores)))
            
        # Coverage
        metric = title + "Coverage"
        real, sims = convert_to_numpy([st.Coverage() for st in stats])

        real /= np.sum(real)

        scores = []
        for sim in sims:
            sim /= np.sum(sim)
        
            abs_val = np.sum(np.abs(real-sim))

            if abs_val < 0.1:
                scores.append(1)
            elif abs_val < 0.35:
                scores.append(2)
            elif abs_val < 0.62:
                scores.append(3)
            elif abs_val < 0.9:
                scores.append(4)
            else:
                scores.append(5)
                
        print( metric, list(zip(names,scores)))
            
        # Duplications
        metric = title + "Duplications"
        real, sims = convert_to_numpy([st.FragmentDuplicationNumber() for st in stats])
        
        real = np.trim_zeros(real, 'b')
        
        len_x = len(real)
        max_x = len(real)
        for i in range(2,len(real)):
            if 0 == real[i]: # Find first zero
                max_x = i
                break
            
        if max_x < 3:
            print( "Duplications in File are all below 50: ", real )
            sys.exit(1)
        
        real = real[2:max_x]/np.sum(real)

        scores = []
        for sim in sims:
            sim = sim / np.sum(sim)
            if len_x < len(sim):
                prop = np.sum(sim[len_x:])
            else:
                prop = 0.0
            
            sim = sim[2:max_x]
            with np.errstate(divide='ignore', invalid='ignore'):
                rel = np.where((0 == real) | (0 == sim), sim == real, np.where(sim > real, real/sim, sim/real))
            mean_rel = np.mean(rel)
            
            if mean_rel > 0.75:
                scores.append(1)
            elif mean_rel > 0.25:
                scores.append(2)
            elif mean_rel > 0.05:
                scores.append(3)
            elif mean_rel > 0.01:
                scores.append(4)
            else:
                scores.append(5)
                
            if prop > 0:
                if prop > 1e-4:
                    scores[-1] += 2
                else:
                    scores[-1] += 1
            scores[-1] = min(5, scores[-1])
                
        print( metric, list(zip(names,scores)))
            
        # Fragment length
        metric = title + "FragLen"
        real, sims = convert_to_numpy([st.InsertLengths() for st in stats])

        real /= np.sum(real)
 
        scores = []
        for sim in sims:
            sim /= np.sum(sim)
        
            abs_val = np.sum(np.abs(real-sim))
            
            if abs_val < 0.01:
                scores.append(1)
            elif abs_val < 0.05:
                scores.append(2)
            elif abs_val < 0.2:
                scores.append(3)
            elif abs_val < 0.6:
                scores.append(4)
            else:
                scores.append(5)
                
        print( metric, list(zip(names,scores)))


def usage():
    print( "Usage: python score_reseq.py [OPTIONS] File File2 [File3 ...]" )
    print( "Scores systematic errors, coverage, duplications and fragment length of File2 [File3 ...] compared to File." )
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

    score_reseq(args,title)
    pass

if __name__ == "__main__":
    main(sys.argv[1:])
    pass

