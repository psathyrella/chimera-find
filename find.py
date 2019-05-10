#!/usr/bin/env python
import collections
import argparse
import sys
import os
import yaml

import chimerafcns

helpstr = """
Input yaml file should be a list of dicts, where each dict has a uid, naive V sequence, and mature V sequence.
One way to make it would be something like this:

seqfos = [
    {'uid' : 'a',
     'naive' : 'CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGAGGTCCCTGAGACTCTCCTGTGCAGCGTCTGGATTCACCTTCAGTAGCTATGGCATGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATGGTATGATGGAAGTAATAAATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAG',
     'mature' : 'CGCAGGCAGCTGGTGCAGTCTGAGGGAGGCGTGGTCCAGCCTGGGAGGTCCCTGAGACTCTCCTGTGTAACGTCTGGATTCTTCTTCAGCAGTTATGGCCTGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCATTTATTTGGTCTGATGGAACTAAGAAATACTACACAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAATTTTAAGAGCACACTGTATCTGCAGATGAACAGCCTGAGAGTCGACGACACGGCTAGGTATTATTGTGTGAGGG',
    },
    {'uid' : 'b',
     'naive' : 'CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCACAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGCAGTGGTGGTTACTACTGGAGCTGGATCCGCCAGCACCCAGGGAAGGGCCTGGAGTGGATTGGGTACATCTATTACAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTCGAGTTACCATATCAGTAGACACGTCTAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACTGCCGCGGACACGGCCGTGTATTACTGTGCGAGAGA',
     'mature' : 'GCAGTACAGCTGCAGCAGTCGGGCCCAGGACTGGTGAAGCCTTCACAGACCCTGTCCCTCAGCTGCACTGTCTCTGGTGACTCCATCAACAATGGTGGTTACTACTGGACCTGGATCCGCCAGCACCCAGGGAAGGGCCTGGAGTGGATTGGGTACATCTATTACAGTGGGCTCACCTACTACAACCCGTCCCTCAGGAGTCGAGTTACCATGTCAGTAGACACGTCTAAAAACCACTTCTCCCTGAGGCTGAGTTTTGTGACTGCCGCGGACACGGCCGTGTATTACTGTGCGAGAGA',
    },
]
with open('example-input.yaml', 'w') as yfile:
    yaml.dump(seqfos, yfile)
"""

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, epilog=helpstr)
parser.add_argument('infile', help='input yaml file with naive/mature V region sequence pairs (see example-input.yaml)')
parser.add_argument('--chunk-len', default=75, type=int)
args = parser.parse_args()

with open(args.infile) as yfile:
    seqfos = yaml.load(yfile, Loader=yaml.Loader)

for sfo in seqfos:
    print chimerafcns.get_chimera_max_abs_diff(sfo['naive'], sfo['mature'], chunk_len=args.chunk_len)
    sys.exit()

# chfo = {uid : {k : v for k, v in zip(('imax', 'max_abs_diff'), utils.get_chimera_max_abs_diff(annotations[uid], iseq=0, chunk_len=args.chunk_len))} for uid in annotations}
# biggest_adiffs = sorted(chfo, key=lambda q: chfo[q]['max_abs_diff'], reverse=True)
# for uid in biggest_adiffs[:5]:
#     print '%-3d  %6.3f' % (chfo[uid]['imax'], chfo[uid]['max_abs_diff'])
#     utils.print_reco_event(annotations[uid])
