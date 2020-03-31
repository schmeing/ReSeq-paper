#!/usr/bin/python

import getopt
import matplotlib.pyplot as pyplot
import numpy
import pysam
from scipy.optimize import leastsq
import sys

def plotInsertSize(bamFile,oFile=False,maxSize=False):
    if False == maxSize:
        maxSize = 1000;
        pass
    if False == oFile:
        oFile = bamFile[:-4]+".pdf"
        pass

    sam = pysam.AlignmentFile(bamFile, "r")
    newStandard = False

    insertSizes = numpy.zeros(maxSize+1, numpy.int32)
    overflow = 0
    underflow = 0
    Nreads = 0
    readLength = 0
    for read in sam.fetch():
        # Can't use template length since definition depends on aligner
        # bwa aligner uses old rightmost of reverse - leftmost of forward (old specification) instead of rightmost of both - leftmost of both (new specification == insert size)
        # bowtie2 does it acording to the new specification
        # New specification: Always use the second read of the pair since the length of that one is needed to calculate the righmost part of the alignment (second read is the one with the higher position)
        # Mapping quality threshold can only be applied to the second read without requiring the sam/bam file to be sorted, bowtie2 assigns mapping quality to pair, but bwa assigns it to the reads individually, there manual sorting is needed to remove contamination from low quality first reads
        # Old specification: Always use reverse read
        if not read.is_unmapped and not read.mate_is_unmapped and 0 != read.template_length and read.mapping_quality > 10 and read.reference_id == read.next_reference_id\
                and (newStandard and (read.next_reference_start < read.reference_start or read.next_reference_start == read.reference_start and read.is_read2)\
                or not newStandard and read.is_reverse):
            curInsertSize = read.reference_end - read.next_reference_start
            if maxSize < curInsertSize:
                overflow += 1
                pass
            elif 0 < curInsertSize:
                insertSizes[curInsertSize] += 1
                pass
            else:
                underflow += 1
                pass
            pass

        Nreads += 1
        readLength += len(read.query_sequence)
        pass

    if overflow:
        print "overflow:", overflow
        pass
    if underflow:
        print "underflow:", underflow
        pass

    # Fit Gaussian
    x = xrange(maxSize+1)
    y = insertSizes.tolist()

    gaussian = lambda p, x: p[2]*numpy.exp(-0.5*((x-p[0])/p[1])**2)
    errfunc  = lambda p, x, y: (y - gaussian(p, x))

    out = leastsq( errfunc, [400.0, 100.0, insertSizes[400]], args=(x,y) )
    fit = out[0]
    print "Gaussian mean insert size:", fit[0]
    print "Gaussian standard deviation:", fit[1]
    print "Mean read length:", float(readLength)/Nreads

    # Plot
    patches = pyplot.hist(x, maxSize+1, (0,maxSize), weights=y, color='red', edgecolor='red')

    pyplot.xlabel("Insert/Fragment size")
    pyplot.ylabel("#Pairs")

    pyplot.plot(x, gaussian(fit, x), 'k-', lw=2, label='fit')
    pyplot.annotate(r'$\mu={:.2f}$'.format(fit[0])+'\n'+r'$\sigma={:.2f}$'.format(fit[1]), xycoords='axes fraction', xy=(0.85,0.9) )

    pyplot.savefig(oFile)
    sam.close()
    pass

def usage():
    print "Usage: python plotInsertSize.py [OPTIONS] File"
    print "Creates a plot of the insert Sizes in File(bam Format)."
    print "  -h, --help            display this help and exit"
    print "  -o, --output FILE     file to which the plot should be written to (File with .bam replaced by .pdf)"
    print "  -r, --maxSize [int]   insert size up to where the plot is plotted (1000)"
    pass

def main(argv):
    try:
        optlist, args = getopt.getopt(argv, 'ho:r:', ['help','output=','maxSize='])
        pass
    except getopt.GetoptError:
	print "Unknown option\n"
        usage()
        sys.exit(2)
        pass

    oFile = False
    maxSize = False
    for opt, par in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
            pass
        elif opt in ("-o", "--output"):
            oFile = par
            pass
        elif opt in ("-r", "--maxSize"):
            try:
                maxSize = int(par)
                pass
            except ValueError:
                print "-r,--maxSize option only accepts integers"
                sys.exit()
                pass
            pass
        pass

    if 1 != len(args):
        print "Wrong number of files. Exactly one file is required.\n"
        usage()
        sys.exit(2)
        pass

    plotInsertSize(args[0],oFile,maxSize)
    pass

if __name__ == "__main__":
    main(sys.argv[1:])
    pass
