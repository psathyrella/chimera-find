#!/usr/bin/env python
import collections
import argparse
import sys
import os
import yaml

import chimerafcns

class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter)
parser.add_argument('infile', help='input yaml file with naive/mature V region sequence pairs (see example-input.yaml)')
parser.add_argument('--chunk-len', default=75, help='length over which to calculate maximum absolute difference (see explanation slides)')
parser.add_argument('--cutoff', default=0.3, help='max-abs-diff value above which we assume most sequences are chimeric')
args = parser.parse_args()

with open(args.infile) as yfile:
    seqfos = yaml.load(yfile, Loader=yaml.Loader)

chinfo = {}
for sfo in seqfos:
    chinfo[sfo['uid']] = chimerafcns.get_chimera_max_abs_diff(sfo['v_naive'], sfo['v_mature'], chunk_len=args.chunk_len)

n_above_cutoff = len([cfo for cfo in chinfo.values() if cfo['max_abs_diff'] > args.cutoff])
chimeric_fraction = n_above_cutoff / float(len(chinfo))
print '  %d / %d = %.3f above chimeric cutoff of %.2f' % (n_above_cutoff, len(chinfo), chimeric_fraction, args.cutoff)
