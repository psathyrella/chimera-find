""" this is all copied from partis/python/utils.py """

ambiguous_bases = ['N', ]
gap_chars = ['.', '-']

# ----------------------------------------------------------------------------------------
def ambig_frac(seq):
    ambig_seq = filter(ambiguous_bases.__contains__, seq)
    return float(len(ambig_seq)) / len(seq)

# ----------------------------------------------------------------------------------------
def hamming_distance(seq1, seq2, return_mutated_positions=False):
    if len(seq1) != len(seq2):
        raise Exception('unequal length sequences %d %d:\n  %s\n  %s' % (len(seq1), len(seq2), seq1, seq2))
    if len(seq1) == 0:
        if return_len_excluding_ambig:
            return 0, 0
        else:
            return 0

    skip_chars = set(ambiguous_bases + gap_chars)

    distance, len_excluding_ambig = 0, 0
    mutated_positions = []
    for ich in range(len(seq1)):  # already made sure they're the same length
        if seq1[ich] in skip_chars or seq2[ich] in skip_chars:
            continue
        len_excluding_ambig += 1
        if seq1[ich] != seq2[ich]:
            distance += 1
            if return_mutated_positions:
                mutated_positions.append(ich)

    if return_mutated_positions:
        return distance, mutated_positions
    else:
        return distance

# ----------------------------------------------------------------------------------------
def get_chimera_max_abs_diff(naive_v_seq, mature_v_seq, chunk_len=75, max_ambig_frac=0.1, debug=False):
    if ambig_frac(naive_v_seq) > max_ambig_frac or ambig_frac(mature_v_seq) > max_ambig_frac:
        if debug:
            print '  too much ambig %.2f %.2f' % (ambig_frac(naive_v_seq), ambig_frac(mature_v_seq))
        return None, 0.

    # if debug:
    #     color_mutants(naive_v_seq, mature_v_seq, print_result=True)
    #     print ' '.join(['%3d' % s for s in isnps])

    _, isnps = hamming_distance(naive_v_seq, mature_v_seq, return_mutated_positions=True)

    max_abs_diff, imax = 0., None
    for ipos in range(chunk_len, len(mature_v_seq) - chunk_len):
        if debug:
            print ipos

        muts_before = [isn for isn in isnps if isn >= ipos - chunk_len and isn < ipos]
        muts_after = [isn for isn in isnps if isn >= ipos and isn < ipos + chunk_len]
        mfreq_before = len(muts_before) / float(chunk_len)
        mfreq_after = len(muts_after) / float(chunk_len)

        if debug:
            print '    len(%s) / %d = %.3f' % (muts_before, chunk_len, mfreq_before)
            print '    len(%s) / %d = %.3f' % (muts_after, chunk_len, mfreq_after)

        abs_diff = abs(mfreq_before - mfreq_after)
        if imax is None or abs_diff > max_abs_diff:
            max_abs_diff = abs_diff
            imax = ipos

    return {'imax' : imax, 'max_abs_diff': max_abs_diff}  # <imax> is break point
