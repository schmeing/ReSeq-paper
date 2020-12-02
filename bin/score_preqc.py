#!/usr/bin/env python3

import csv
import getopt
import numpy as np
import os.path
from scipy.stats import spearmanr
import simplejson as json
import sys

def load_preqc_datafiles(preqc_files):
    # load the data
    data = []
    deserial_fail = []
    for f in preqc_files:
        f = os.path.abspath(f)
        if( os.path.getsize(f) <= 0 ):
            print("Warning: empty file '%s' ... skipping"%f)
            continue
        if( f in [d['file'] for d in data] ):
            print("Warning: duplicate file '%s' ... skipping"%f)
            continue
        try:
            deserial = json.load(open(f, 'r'))
        except ValueError:
            deserial_fail.append(f)
            continue
        data.append(deserial)
        data[-1]['file'] = f
    for failed in deserial_fail:
        print("Warning: failed to de-serialize file:", failed)
    return data

def n50(values):
    """Return the N50 of the list of numbers"""
    values.sort(reverse=True)
    target = sum(values) / 2.
    total = 0
    for v in values:
        total += v
        if total >= target:
            return v
    return 0

def count_matches(real, sim, dist):
    # Fill overhang to both sides
    full = np.zeros(len(real)+2*dist)
    full[dist:-dist] = real
    
    for d in range(1,dist+1):
        full[dist-d] = 2*full[dist] - full[dist+1]
        full[-dist-1+d] = 2*full[-dist-1] - full[-dist-2]
        
    # Calculate minimum and maximum in the range of dist for each value
    cmin = real
    cmax = real
    for d in range(1,dist+1):
        cmin = np.minimum(cmin,full[dist-d:len(full)-dist-d])
        cmax = np.maximum(cmax,full[dist-d:len(full)-dist-d])
        
        cmin = np.minimum(cmin,full[dist+d:len(full)-dist+d])
        cmax = np.maximum(cmax,full[dist+d:len(full)-dist+d])

    # count
    return np.sum( (sim >= cmin) & (sim <= cmax) )    

def score_preqc(preqcFiles, title = ""):
    names = []
    for pf in preqcFiles:
        if 'gz' == pf.split('.')[-1]:
            names.append( pf.split('.',2)[0] )
        else:
            names.append( pf.rsplit('.',1)[0] )

    for i, name in enumerate(names):
        if '/' in name:
            names[i] = name.rsplit('/',1)[1].split('.')[0]
    
    names = names[1:]

    preqc = load_preqc_datafiles(preqcFiles)
    
    if len(title):
        print(title)
        title = title + " "

    if 1 < len(preqc):
        # Quality Scores
        metric = title + "Quality"

        real = np.array(preqc[0]['QualityScores']['mean_quality'])
        sims = [ np.array(sim['QualityScores']['mean_quality']) for sim in preqc[1:] ]
        
        for i in range(len(sims)):
            sims[i].resize(len(real))
        
        scores = []
        for name, sim in zip(names,sims):
            abs_val = np.sum(np.abs(sim-real))
            
            #print( metric, name, abs_val )
            
            if abs_val < 20:
                scores.append(1)
            elif abs_val < 53:
                scores.append(2)
            elif abs_val < 100:
                scores.append(3)
            elif abs_val < 500:
                scores.append(4)
            else:
                scores.append(5)
            
        print( metric, list(zip(names,scores)))
        
        # Error Rates
        metric = title + "Errors"

        real = np.array(preqc[0]['ErrorsPerBase']['error_count'])/np.array(preqc[0]['ErrorsPerBase']['base_count'])
        sims = [ np.array(sim['ErrorsPerBase']['error_count'])/np.array(sim['ErrorsPerBase']['base_count']) for sim in preqc[1:] ]
        
        for i in range(len(sims)):
            sims[i].resize(len(real))
        
        mean_rate = np.mean(real)
        real /= mean_rate
        
        scores = []
        for name, sim in zip(names,sims):
            sim /= mean_rate
            abs_val = np.sum(np.abs(sim-real))

            #print( metric, name, abs_val )
            
            if abs_val < 20:
                scores.append(1)
            elif abs_val < 35:
                scores.append(2)
            elif abs_val < 65:
                scores.append(3)
            elif abs_val < 100:
                scores.append(4)
            else:
                scores.append(5)
            
        print( metric, list(zip(names,scores)))
        
        # Assembly continuity
        metric = title + "N50"

        real = np.array([ n50(asm_sim['walk_lengths']) for asm_sim in preqc[0]['SimulateAssembly'] ]) #, dtype=np.float64
        sims = [ np.array([ n50(asm_sim['walk_lengths']) for asm_sim in sim['SimulateAssembly'] ]) for sim in preqc[1:] ]

        for i in range(len(sims)):
            sims[i].resize(len(real))

        scores = []
        for name, sim in zip(names,sims):
            count1 = count_matches(real, sim, 1)
            count2 = count_matches(real, sim, 2)
            count4 = count_matches(real, sim, 4)
            
            real_sign = real[1:] - real[:-1]
            real_sign = np.where(real_sign < 0, -1, np.where(real_sign > 0, 1, 0))
            
            sim_sign = sim[1:] - sim[:-1]
            sim_sign = np.where(sim_sign < 0, -1, np.where(sim_sign > 0, 1, 0))
            
            sign_diff = np.abs(real_sign - sim_sign)
            max_diff = 0
            cur_diff = 0
            cur_sign = 0
            for sdiff, rsign in zip(sign_diff, real_sign):
                if sdiff == 0 or abs(rsign - cur_sign) == 2 or (cur_diff == 0 and sdiff == 1):
                    cur_diff = 0
                    cur_sign = 0
                else:
                    cur_diff += 1
                    max_diff = max(max_diff, cur_diff)
                    cur_sign = rsign
            
            #print( metric, name, count1, count2, count4, max_diff )
            
            if count1 >= 13 and count2 >= 14 and count4 >= 15:
                scores.append(1)
            elif count1 >= 5 and count2 >= 10 and count4 >= 12:
                scores.append(2)
            elif count2 >= 5 and count4 >= 9:
                scores.append(3)
            elif count4 >= 6:
                scores.append(4)
            else:
                scores.append(5)
            
            if max_diff >= 4:
                scores[-1] = min(5, scores[-1] + 1)
            
        print( metric, list(zip(names,scores)))



def usage():
    print( "Usage: python score_preqc.py [OPTIONS] File File2 [File3 ...]" )
    print( "Scores quality values, error rates and assembly continuity from preqc File2 [File3 ...] compared to File." )
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

    score_preqc(args,title)
    pass

if __name__ == "__main__":
    main(sys.argv[1:])
    pass

