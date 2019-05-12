#!/usr/bin/env python
import collections
import argparse
import sys
import os
import yaml

import chimerafcns

helpstr = """
Calculates a metric (max-abs-diff) that helps to identify chimeric sequences.
In each sequence, the metric looks for positions within V at which the SHM rate is dramatically different to the left vs to the right.
We find the position that maximizes this difference, and call that difference the maximum absolute difference, or max-abs-diff.
See explanation-slides.pdf for details.

This script calculates max-abs-diff for each input sequence.
The distribution of these values can be plotted and compared to the distributions in explanation-slides.pdf.
In addition, this script prints the fraction of input sequences that have "very high" max-abs-diff (above a threshold specified by --cutoff).
These are sequences that are quite likely to be chimeric, and thus the larger the fraction that this represents of your repertoire, the more likely it is that you have an atypically large numbers of chimeric sequences.
A repertoire with zero chimeric sequences will have a very small such fraction, less than around one percent.

Input yaml file should be a list of dicts, where each dict has a uid, naive V sequence, and mature V sequence.
One way to make it would be something like this:

seqfos = [
    {'uid' : 'a',
     'v_naive' : 'CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGAGGTCCCTGAGACTCTCCTGTGCAGCGTCTGGATTCACCTTCAGTAGCTATGGCATGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATGGTATGATGGAAGTAATAAATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAG',
     'v_mature' : 'CGCAGGCAGCTGGTGCAGTCTGAGGGAGGCGTGGTCCAGCCTGGGAGGTCCCTGAGACTCTCCTGTGTAACGTCTGGATTCTTCTTCAGCAGTTATGGCCTGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCATTTATTTGGTCTGATGGAACTAAGAAATACTACACAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAATTTTAAGAGCACACTGTATCTGCAGATGAACAGCCTGAGAGTCGACGACACGGCTAGGTATTATTGTGTGAGGG',
    },
    {'uid' : 'b',
     'v_naive' : 'CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCACAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGCAGTGGTGGTTACTACTGGAGCTGGATCCGCCAGCACCCAGGGAAGGGCCTGGAGTGGATTGGGTACATCTATTACAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTCGAGTTACCATATCAGTAGACACGTCTAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACTGCCGCGGACACGGCCGTGTATTACTGTGCGAGAGA',
     'v_mature' : 'GCAGTACAGCTGCAGCAGTCGGGCCCAGGACTGGTGAAGCCTTCACAGACCCTGTCCCTCAGCTGCACTGTCTCTGGTGACTCCATCAACAATGGTGGTTACTACTGGACCTGGATCCGCCAGCACCCAGGGAAGGGCCTGGAGTGGATTGGGTACATCTATTACAGTGGGCTCACCTACTACAACCCGTCCCTCAGGAGTCGAGTTACCATGTCAGTAGACACGTCTAAAAACCACTTCTCCCTGAGGCTGAGTTTTGTGACTGCCGCGGACACGGCCGTGTATTACTGTGCGAGAGA',
    },
]
with open('example-input.yaml', 'w') as yfile:
    yaml.dump(seqfos, yfile)
"""

class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter, epilog=helpstr)
parser.add_argument('infile', help='input yaml file with naive/mature V region sequence pairs (see example-input.yaml)')
parser.add_argument('--chunk-len', default=75, help='length over which to calculate maximum absolute difference (see explanation slides)')
parser.add_argument('--cutoff', default=0.3, help='max-abs-diff value above which we assume most sequences are chimeric')
args = parser.parse_args()

with open(args.infile) as yfile:
    seqfos = yaml.load(yfile, Loader=yaml.Loader)

chinfo = {}
for sfo in seqfos:
    chinfo[sfo['uid']] = chimerafcns.get_chimera_max_abs_diff(sfo['v_naive'], sfo['v_mature'], chunk_len=args.chunk_len)

n_above_cutoff = len([_ for cfo in chinfo.values() if cfo['max_abs_diff'] > args.cutoff])
chimeric_fraction = n_above_cutoff / float(len(chinfo))
print '  %d / %d = %.3f above chimeric cutoff of %.2f' % (n_above_cutoff, len(chinfo), chimeric_fraction, args.cutoff)
