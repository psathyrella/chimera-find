Calculates a metric (max-abs-diff) that helps to identify chimeric sequences.
In each sequence, we look for positions within V at which the SHM rate is dramatically different to the left vs to the right.
We find the position that maximizes this difference, and call that difference the maximum absolute difference, or max-abs-diff.
See explanation-slides.pdf for details.

The script `find.py` calculates max-abs-diff for each input sequence.
The distribution of these values can be plotted and compared to the distributions in explanation-slides.pdf.
In addition, the script prints the fraction of input sequences that have "very high" max-abs-diff (above a threshold specified by --cutoff).
These are sequences that are quite likely to be chimeric, and thus the larger the fraction that this represents of your repertoire, the more likely it is that you have an atypically large numbers of chimeric sequences.
A repertoire with zero chimeric sequences will typically have a fraction of around one percent.

Two example simulation samples are included in `examples/`, one with no chimeras, and one consisting entirely of chimeras.
The former has a value for the fraction described above of 0.012 (1.2%), while the latter has 0.21 (21%).
These samples correspond to the simulation plots in explanation-slides.pdf.

Input yaml file should be a list of dicts, where each dict has a uid, naive V sequence, and mature V sequence.
One way to make it would be something like the following python code:

```
import yaml

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
```
