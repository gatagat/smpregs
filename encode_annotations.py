#!/usr/bin/env python

import os
from pyfasta import Fasta
from region_utils import regions_reader

def create_encoded_fasta(input_filename, output_filename, flatten=True, codes=None):
    codes_pool = 'abcdefghijklmnopqrstuvwxyz'
    if codes is None:
        codes = {}
    else:
        codes = codes.copy()
        codes_pool = sorted(list(set(codes_pool) - set(codes.values())))

    print 'Creating fasta'
    with open(output_filename, 'w') as fw:
        chrom = None
        for region in regions_reader(input_filename):
            if chrom != region.chr:
                if chrom is not None:
                    fw.write('\n')
                fw.write('>' + region.chr + '\n')
                chrom = region.chr
                assert region.start == 0
                old_stop = region.stop
            else:
                assert old_stop == region.start
                old_stop = region.stop
            if not region.name in codes:
                codes[region.name] = codes_pool[0]
                codes_pool = codes_pool[1:]
            fw.write(codes[region.name] * (region.stop - region.start))
    
    print 'Encoded as: %s' % str(codes)

    print 'Flattening'
    if os.path.exists(output_filename + '.gdx'): # cleanup non-inplace idx if any
        os.unlink(output_filename + '.gdx')
    Fasta(output_filename, flatten_inplace=True)


if __name__ == '__main__':
    import sys

    if len(sys.argv) != 4:
        print 'Usage: %s input encoding output' % sys.argv[0]
        print
        print 'input  - BED file with genome annotations'
        print 'encoding - tab separated file with two columns: "annotation  coding character"'
        print 'output - FASTA file where genome annotations are stored encoded by single characters'
        print
        sys.exit(-1)

    input_filename = sys.argv[1]
    encoding_filename = sys.argv[2]
    output_filename = sys.argv[3]
    codes = {}
    with open(encoding_filename) as fr:
        for line in fr:
            toks = line.rstrip('\r\n').split('\t')
            codes[toks[0]] = toks[1]
    create_encoded_fasta(input_filename, output_filename, codes=codes)

