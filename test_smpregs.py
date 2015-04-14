import numpy as np
import os
from tempfile import mktemp
from region_utils import AllowedSpace, RegionAcceptorApproxGC, count_g_and_c, RegionAcceptorApproxHistogram, histogram_intersection, KmerHistogram
from smpregs import sample_regions, get_genome #, _setup_log
from kmers import count_kmers, all_kmers


def create_regions_file(length, cnt, fasta, no_ns=False):
    tempfile = mktemp()
    chroms = fasta.keys()
    with open(tempfile, 'w') as fw:
        for i in range(1, 1+cnt):
            chrom = np.random.choice(chroms)
            pos = np.random.randint(0, len(fasta[chrom])-length-1)
            if no_ns:
                seq = fasta[chrom][pos:pos+length]
                while not np.in1d(list(str(seq).upper()), list('ACGT')).all():
                    print 'N-s encountered, drawing a new region.'
                    chrom = np.random.choice(chroms)
                    pos = np.random.randint(0, len(fasta[chrom])-length-1)
                    seq = fasta[chrom][pos:pos+length]
            fw.write('%s\t%d\t%d\t%s\n' % \
                    (chrom, pos, pos+length, 'reg%03d' % i))
    return tempfile


def region_length(r):
    return r.stop - r.start


def region_sequence(fasta, r):
    return fasta[r.chrom][r.start:r.stop]


def region_gc(fasta, r):
    return count_g_and_c(region_sequence(fasta, r))


def test_sample_regions_simple():
    genome_fasta = get_genome('dm3')
    regions_file = create_regions_file(301, 10, genome_fasta)
    for input_region, random_region in sample_regions(
            regions_file, AllowedSpace(genome_fasta), [], genome_fasta):
        assert region_length(input_region) == region_length(random_region)
        assert input_region.chrom == random_region.chrom
    os.unlink(regions_file)


def test_sample_regions_gc():
    #_setup_log()
    np.random.seed(1234L)
    genome_fasta = get_genome('dm3')
    regions_file = create_regions_file(301, 100, genome_fasta)
    for input_region, random_region in sample_regions(
            regions_file, AllowedSpace(genome_fasta),
            [(RegionAcceptorApproxGC, dict(threshold=10))], genome_fasta):
        assert region_length(input_region) == region_length(random_region)
        input_gc = region_gc(genome_fasta, input_region)
        random_gc = region_gc(genome_fasta, random_region)
        print input_region
        print input_gc, random_gc
        assert input_gc - 10 <= random_gc
        assert input_gc + 10 >= random_gc
    os.unlink(regions_file)


def test_sample_regions_gc_relative():
    #_setup_log()
    np.random.seed(1234L)
    genome_fasta = get_genome('dm3')
    for size in [100, 200, 300]:
        regions_file = create_regions_file(size, 10, genome_fasta)
        for input_region, random_region in sample_regions(
                regions_file, AllowedSpace(genome_fasta),
                [(RegionAcceptorApproxGC, dict(threshold=0.1))], genome_fasta):
            assert region_length(input_region) == region_length(random_region)
            input_gc = region_gc(genome_fasta, input_region)
            random_gc = region_gc(genome_fasta, random_region)
            print input_region
            print input_gc, random_gc, size
            assert input_gc - 0.1 * size <= random_gc
            assert input_gc + 0.1 * size >= random_gc
    os.unlink(regions_file)


def test_sample_regions_kmers():
    #_setup_log()
    np.random.seed(1234L)
    genome_fasta = get_genome('dm3')
    regions_file = create_regions_file(301, 10, genome_fasta, no_ns=True)
    kmer_k = 2
    kmer_keys = all_kmers(kmer_k)
    for input_region, random_region in sample_regions(
            regions_file, AllowedSpace(genome_fasta),
            [(RegionAcceptorApproxHistogram, dict(threshold=20, histogram=KmerHistogram(fasta=genome_fasta, k=kmer_k), features_per_nt=2))], genome_fasta):
        assert region_length(input_region) == region_length(random_region)
        input_kmers = count_kmers(kmer_k, str(genome_fasta[input_region.chrom][input_region.start:input_region.stop]).upper())
        random_kmers = count_kmers(kmer_k, str(genome_fasta[random_region.chrom][random_region.start:random_region.stop]).upper())
        print input_region
        print input_kmers, random_kmers
        errors = 2*(input_region.stop - input_region.start) - histogram_intersection(dict(zip(kmer_keys, input_kmers)), dict(zip(kmer_keys, random_kmers)))
        print errors
        assert errors < 20
    os.unlink(regions_file)
