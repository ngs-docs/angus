#! /usr/bin/env python
"""
This script calculates the mismatch profile by read location for SAM files,
given the SAM file and the reference genome.

(This was developed for https://peerj.com/preprints/890/.)
"""
import sys
import argparse
import screed
import math


MAX_SEQ_LEN = 5000


def ignore_at(iter):
    for item in iter:
        if item.startswith('@'):
            continue
        yield item


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('genome')
    parser.add_argument('samfile')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    
    args = parser.parse_args()

    # load the sequences from the genome.
    genome_dict = dict([ (record.name, record.sequence) for record in \
                        screed.open(args.genome) ])

    # initialize the array of mismatch counts.
    positions = [0] * MAX_SEQ_LEN
    lengths = []

    # walk through the SAM file, calculating mismatches.
    n = 0
    n_rev = n_fwd = 0
    for samline in ignore_at(open(args.samfile)):
        n += 1
        if n % 100000 == 0:
            print >>sys.stderr, '...', n

        readname, flags, refname, refpos, _, _, _, _, _, seq = \
                  samline.split('\t')[:10]

        if refname == '*' or refpos == '*':
            # (don't count these as skipped)
            continue
        
        refpos = int(refpos)

        ref = genome_dict[refname][refpos-1:refpos+len(seq) - 1]

        # for each read, look where the mismatches are by position.
        errors = []
        for pos, (a, b) in enumerate(zip(ref, seq)):
            if a != b:
                # see http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#sam-output - '16' is the revcomp flag.
                if int(flags) & 16:
                    pos = len(seq) - pos - 1
                    n_rev += 1
                else:
                    n_fwd += 1
                positions[pos] += 1
        lengths.append(len(seq))

    # normalize for length
    lengths.sort()
    max_length = lengths[-1]

    length_count = [0]*max_length
    for j in range(max_length):
        length_count[j] = sum([1 for i in lengths if i >= j])

    # output!
    args.outfile.write('position,error_count,error_fraction\n')
    for n, i in enumerate(positions[:max_length]):
        print >>args.outfile, "%s,%s,%s" % (n, i, float(i) /
                                            float(length_count[n]))

    print >>sys.stderr, "error rate: %.2f%%" % \
          (100.0 * sum(positions) / float(sum(lengths)))
    print >>sys.stderr, 'logratio of fwd to rev: %.2f' % (math.log(n_fwd / float(n_rev), 2))

if __name__ == '__main__':
    main()
