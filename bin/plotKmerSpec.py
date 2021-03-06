#!/usr/bin/python

import csv
import getopt
from itertools import izip
from math import ceil
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PatchCollection
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
from matplotlib.ticker import LogFormatter
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import os
import sys


text_size = 20
tick_size = 16
colour_scheme = ["#D92120","#488BC2","#7FB972","#E6642C","#781C81","#D9AD3C","#BBBBBB","#4065B1"]

line_width = 4

def findIdOfMinimum(y, interval=15, max_id=int(1e6)):
    min_id = 1
    while min_id+1 < len(y) and y[min_id] < y[min_id + 1]: # Find the maximum of the systematic error peak (which might be already at 1)
        min_id += 1
    while min_id+1 < len(y) and y[min_id] > min(y[min_id+1:min_id+2+interval]): # Find the minimum splitting systematic errors and signal peak
        min_id = y.index(min(y[min_id+1:min_id+2+interval]))

    if min_id < len(y) and min_id < max_id:
        return min_id
    else:
        return min(len(y), max_id)-1

def plot(pdf, xtitle, ytitle, names, hists, xlim=False, ylim=False, log=False, legend='upper right', normalize=False, title="", sys_deviation=[], cov_deviation=[], vline=0):
    plt.close()
    
    if vline:
        plt.plot(hists[0][0], hists[0][1], linewidth=1, marker=None, linestyle='-', label='')
        plt.vlines(vline, 0, 1e99, linewidth=3)
    
    for (name, hist, col) in izip(names, hists, colour_scheme):
        x_vals = hist[0]
        y_vals = hist[1]
        if normalize:
            y_sum = sum(y_vals)
            y_vals = [y/y_sum for y in y_vals]
        plt.plot(x_vals, y_vals, linewidth=line_width, marker=None, linestyle='-', label=name, color=col)
        pass

    ax = plt.gca()
    if log:            
        ax.set_yscale('log', nonposy='clip')
        ax.yaxis.set_major_formatter(LogFormatter())
        #ax.ticklabel_format(useMathText=False, style='plain')
    else:
        # Set at what range of exponents they are not plotted in exponential format for example (-3,4): [0.001-1000[
        ax.get_yaxis().get_major_formatter().set_powerlimits((-3,4))

    if xlim:
        ax.set_xlim(right=xlim)
        
    if ylim:
        ax.set_ylim(ylim[0], ylim[1])
    
    plt.xlabel(xtitle, fontsize=text_size)
    plt.ylabel(ytitle, fontsize=text_size)
    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)
    ax.yaxis.offsetText.set_fontsize(tick_size)
    
    if title != "":
        plt.title(title, {'fontsize':text_size})
    
    if legend != 'none':
        plt.legend(loc=legend, shadow=True)
        pass
    
    if len(sys_deviation):
        if legend == 'none':
            x_coord = 0.64
            y_coord = 0.93
        else:
            if log: 
                x_coord = 0.39
                y_coord = 0.93
            else:
                x_coord = 0.64
                y_coord = 0.54
        for i in range(len(sys_deviation)):
            plt.annotate(r'$R={:.2f}$'.format(sys_deviation[i]).replace('e+0', 'e'), xycoords='axes fraction', xy=(x_coord,y_coord-i*0.06), color=colour_scheme[i+2], fontsize=text_size )
            plt.annotate(r'|', xycoords='axes fraction', xy=(x_coord+0.205,y_coord-i*0.06), color=colour_scheme[i+2], fontsize=text_size )
            plt.annotate(r'${:.2f}$'.format(cov_deviation[i]).replace('e+0', 'e'), xycoords='axes fraction', xy=(x_coord+0.23,y_coord-i*0.06), color=colour_scheme[i+2], fontsize=text_size )
        
    plt.tight_layout()
    pdf.savefig()

    return

def plotKmerSpec(specFiles, oFile, k, xlim, title, log, legend):
    names=[]
    for i, sf in enumerate(specFiles):
        names.append( sf.rsplit('.',1)[0] )
    
    if not oFile:
        oFile = names[0] + ".pdf"

    for i, name in enumerate(names):
        if '/' in name:
            names[i] = name.rsplit('/',1)[1].split('.')[0][0:10]

    y_min=sys.maxsize
    y_max=0

    specs = []
    for sf in specFiles:
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
            
            if xlim:
                max_id = sum(np.array(x) <= xlim)
            else:
                max_id = len(x)
                
            if log:
                y_min = min(y_min, min(y[:max_id]))
                y_max = max(y_max, max(y[:max_id]))
            else:
                y_min = min(y_min, min(y[:max_id]))
                y_max = max(y_max, max(y[findIdOfMinimum(y, interval=5,max_id=max_id):max_id]))
            
            specs.append((x,y))

    if specs:
        y_min *= 0.9
        y_max *= 1.1
        
        if xlim:
            x_minimum = specs[0][0][ findIdOfMinimum(specs[0][1], xlim//50) ]
        else:
            x_minimum = specs[0][0][ findIdOfMinimum(specs[0][1]) ]
        max_x = max([max(specs[i][0]) for i in range(len(specs))])
        real = np.zeros(max_x+1)
        real[specs[0][0]] = specs[0][1]
        
        sys_deviation = []
        cov_deviation = []
        for spec in specs[1:len(specs)]:
            sim = np.zeros(max_x+1)
            sim[spec[0]] = spec[1]

            sys_deviation.append( np.sum(np.abs((sim[1:x_minimum+1]-real[1:x_minimum+1])/real[1:x_minimum+1])) )
            cov_deviation.append( np.sum(np.abs(sim[x_minimum+1:]-real[x_minimum+1:])) )
            
        if 1 < len(sys_deviation):
            sys_deviation = [d/sys_deviation[0] for d in sys_deviation[1:]]
            cov_deviation = [d/cov_deviation[0] for d in cov_deviation[1:]]
        
        with PdfPages(oFile) as pdf:
            plt.ioff()

            plot( pdf, '-'.join([k, "mer frequency x"]), ''.join(["# unique ",k,"-mers with frequency x"]), names, specs, xlim=xlim, ylim=[y_min, y_max], log=log, legend=legend, title=title, sys_deviation=sys_deviation, cov_deviation=cov_deviation, vline=x_minimum )
    return

def usage():
    print "Usage: python plotKmerSpec.py [OPTIONS] File"
    print "Plots the DataStats from an boost archive from readar."
    print "  -h, --help            display this help and exit"
    print "  -k, --kvalue          value set for k in axis labels"
    print "  -o, --output          define plotting output file [File with ending pdf]"
    print "  -x, --xlim            maximum x to be plotted"
    pass

def main(argv):
    try:
        optlist, args = getopt.getopt(argv, 'hik:ln:o:x:', ['help','nolegend','kvalue=','linear','output=','name=','xlim='])
        pass
    except getopt.GetoptError:
        print "Unknown option\n"
        usage()
        sys.exit(2)
        pass

    oFile = ''
    k = 'k'
    xlim = False
    log = True
    legend = 'upper right'
    title = ""
    for opt, par in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
            pass
        elif opt in ("-i", "--nolegend"):
            legend = 'none'
            pass
        elif opt in ("-k", "--kvalue"):
            k = par
            pass
        elif opt in ("-l", "--linear"):
            log = False
            pass
        elif opt in ("-n", "--name"):
            title = par
            pass
        elif opt in ("-o", "--output"):
            oFile = par
            pass
        elif opt in ("-x", "--xlim"):
            xlim = int(par)
            pass
        pass

    if 1 > len(args) or len(args) > 8:
        print "Wrong number of arguments. Only one to eight files are supported.\n"
        usage()
        sys.exit(2)
        pass

    plotKmerSpec(args, oFile, k, xlim, title, log, legend)
    pass

if __name__ == "__main__":
    main(sys.argv[1:])
    pass
