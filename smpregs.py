#!/usr/bin/env python
#
# Sample regions matching input regions (1:1)
#  - according to GC-content
#  - according to genomic distribution at the peak
#  - accroding to approx. genomic distribution
#  - according to approx. k-mer content (for small k)
#  - according to motif count
#
# Any output matched region is placed on the same chromosome as input region
# and has to be inside the given allowed space.
#
# Exclude overlaps of the generated regions. Exclude input regions.
#

# XXX: probably always we want regions without N's - check it always?

#TODO: now check the random regions selected accroding to peak summit genomic distribution
#    ? does it still look like before (large gd regions in active seqs, and small gd regions in inactive)?

#TODO: allow for passing "-" to read regions (and other inputs) from stdin - save it to tempdir to be able to go through it twice

import logging
import sys
from pyfasta import Fasta
from region_utils import regions_reader, AllowedSpace, generate, \
    RegionAcceptorApproxGC, RegionAcceptorGenomicAnnotation, RegionAcceptorApproxHistogram, GenomicAnnotationsHistogram, KmerHistogram, RegionAcceptorNoNs
from region_utils import get_log
import socket


def _setup_log(level=logging.INFO):
    logging.basicConfig()
    log = logging.getLogger()
    log.handlers = []
    try:
        handler = logging.StreamHandler(stream=sys.stderr)
    except TypeError:
        handler = logging.StreamHandler(strm=sys.stderr)
    formatter = logging.Formatter(
            fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y.%m.%d %I:%M:%S %p'
        )
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(level)


def get_genome(genome_assembly, tempdir):
    """
    Return Fasta object with the required genome.
    """
    logger = get_log('generate')
    if socket.gethostname().split('.')[0].upper() == 'MPBA02':
        fasta_filename = '/Users/kazmar/tmp/%s.fa' % genome_assembly
    else:
        fasta_filename = '/groups/stark/kazmar/data/genomes/%s.fa' % genome_assembly
    logger.debug('Getting genome from %s', fasta_filename)
    return Fasta(fasta_filename)


def parse_filters(filters, genome_fasta, genomic_annotations=None):
    """
    Parse input filters to acceptors.

    Parameters:
    ===========
    filters: list of str
        - string specifications of acceptors
    genomic_annotations: filename
        - Fasta file with encoded genomic annotations.

    Return:
    =======
    accpetors: list of Acceptor objects
    """
    logger = get_log('generate')
    all_acceptor_classes = dict(
            GC=RegionAcceptorApproxGC,
            GAPos=RegionAcceptorGenomicAnnotation,
            GAHist=RegionAcceptorApproxHistogram,
            KMer=RegionAcceptorApproxHistogram)
    acceptors = []
    logger.debug('Parsing filters: %s', str(filters))
    for f in filters:
        try:
            filter_name, filter_args = f.split(':', 1)
            if not filter_name in all_acceptor_classes:
                raise ValueError('Unknown filter name %s.' % filter_name)
            filter_opts = []
            kmer_k = 2
            for a in filter_args.split(','):
                k, v = a.split('=', 1)
                if filter_name == 'KMer' and k == 'k':
                    kmer_k = int(v)
                else:
                    filter_opts += [(k, float(v))]
            if filter_name == 'KMer':
                filter_opts += [('histogram', KmerHistogram(fasta=genome_fasta, k=kmer_k))]
                filter_opts += [('features_per_nt', 2)]
            if filter_name.startswith('GA'):
                if genomic_annotations is None:
                    raise ValueError('Genomic annotations required for filter %s' % filter_name)
                if filter_name == 'GAPos':
                    filter_opts += [('filename', genomic_annotations)]
                elif filter_name == 'GaHist':
                    filter_opts += [('histogram',
                        GenomicAnnotationsHistogram(genomic_annotations))]
                    filter_opts += [('features_per_nt', 1)]
                else:
                    assert False
        except Exception, e:
            logger.exception('Error parsing filter %s' , f)
            raise ValueError('Error parsing filter %s (%s)\n' % (f, str(e)))
        acceptors += [(all_acceptor_classes[filter_name], dict(filter_opts))]
    return acceptors


import contextlib
@contextlib.contextmanager
def output_file_wrapper(filename=None):
    """
    Context manager for either output file or stdout.
    """
    if filename is None:
        yield sys.stdout
    else:
        writer = open(filename, 'w')
        yield writer
        writer.close()

@contextlib.contextmanager
def tempdir():
    from tempfile import mkdtemp
    tmp = mkdtemp()
    yield tmp
    from shutil import rmtree
    rmtree(tmp)


def output_region(stream, region):
    """
    Output region to stream.
    """
    if region.name is None:
        s = '\t'.join([region.chr, str(region.start), str(region.stop)])
    else:
        s = '\t'.join([region.chr, str(region.start), str(region.stop), region.name])
    stream.write(s + '\n')


def sample_regions(regions_file, allowed_space, acceptors, fasta):
    """
    Generator providing random regions that match input regions.

    Output regions are obtained as follows:
    1. sample locations in the genome uniformly
    2. check whether the generated region is inside the allowed space
    3. check whether it fulfils all the acceptors

    Parameters:
    ===========
    regions_file:
    allowed_space: AllowedSpace object
        - Space to which candidate regions are restricted.
    acceptors: list of Acceptor objects.
        - Accept/reject candidate regions.
    fasta: Fasta object
        - Provides access to sequences in regions_file or allowed_space.

    Returns:
    ========
    Yields tuples of input and sampled regions:
    (input_region, matching_random_region)
    """
    logger = get_log('generate')
    for input_region in regions_reader(regions_file):
        acceptor_instances = []
        for acceptor in acceptors:
            acceptor_instances += [acceptor[0](
                    template=input_region,
                    fasta=fasta,
                    **acceptor[1])]
        accepted = False
        while not accepted:
            candidate = generate(input_region, allowed_space)
            accepted = True
            for acceptor in acceptor_instances:
                if not acceptor.accept(candidate):
                    accepted = False
                    logger.info('REJ %s on %s(%s)', candidate, acceptor.__class__.__name__, acceptor.reason)
                    break
        if accepted:
            logger.info('ACC %s', candidate)
            allowed_space.remove(candidate)
            yield input_region, candidate


if __name__ == '__main__':
    import argparse
    import textwrap

    parser = argparse.ArgumentParser(
            description='Sample random regions one-by-one matching given input regions.',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=textwrap.dedent('''\
                Filters
                  One or multiple filter(s) can be specified to enforce properties
                  in which the random regions have to be similar to their matching
                  regions.

                  Allowed filters are: GC, GAPos, GAHist, and KMer.

                  GC:threshold=10 allows at most 10 more/less of GC nucleotides.
                  GAPos:pos=101 enforces equal genomic annotation at position 101 in the sequence.
                  GAHist:threshold=150 allows at most 150 errors when matching histograms of genomic annotations.
                  KMer:k=2,threshold=50 analog. to GAHist but for k-mer sequence content.

                  All filters have to be fulfilled at once (logical AND). Multiple filters of the
                  same kind are allowed.
                '''))
    parser.add_argument('-r', '--regions', dest='regions', required=True,
            action='store', default=None, help='Regions BED file')
    parser.add_argument('-o', '--output', dest='output', required=False,
            action='store', default=None, help='Output BED file')
    parser.add_argument('-i', '--include', dest='include', required=False,
            action='store', default=None, help='Include BED file. Generated  \
            random regions will be fully contained inside one of the regions \
            from this file. If not specified, the whole genome excluding \
            input regions is taken.')
    parser.add_argument('-e', '--exclude', dest='exclude', required=False,
            action='store', default=None, help='Exclude BED file. Generated \
            random regions will _NOT_ overlap any of the regions from this  \
            file.')
    parser.add_argument('-g', '--genome-assembly', dest='genome_assembly',
            required=False, action='store', default='dm3', help='Assembly of \
            the genome [Default: dm3]')
    parser.add_argument('-n', '--genomic-annotations', dest='genomic_annotations',
            required=False, action='store', default=None, help='Genomic \
            annotations FASTA file. Use encode_genomic_annotations.py to create \
            it.')
    parser.add_argument('filters', action='store', nargs='*', help='Filters. \
            See below.')
    parser.add_argument('-v', '--verbose', action='count', default=0)
    opts = parser.parse_args()

    loglevel = max(logging.DEBUG, logging.WARNING - opts.verbose*10)
    _setup_log(level=loglevel)
    logger = get_log('main')
    logger.debug('Logging started at level %d', loglevel)

    with tempdir() as tmp:
        genome_fasta = get_genome(opts.genome_assembly, tmp)
        acceptors = parse_filters(opts.filters, genome_fasta, opts.genomic_annotations)
        acceptors = [(RegionAcceptorNoNs, {})] + acceptors
        allowed_space_opts = {}
        if opts.include is not None:
            allowed_space_opts['include'] = regions_reader(opts.include)
        if opts.exclude is None:
            allowed_space_opts['exclude'] = regions_reader(opts.regions)
        else:
            allowed_space_opts['exclude'] = regions_reader([opts.regions, opts.exclude])
        allowed_space = AllowedSpace(fasta=genome_fasta, **allowed_space_opts)
        with output_file_wrapper(opts.output) as fw:
            for _, region in sample_regions(
                    opts.regions, allowed_space, acceptors, genome_fasta):
                output_region(fw, region)
