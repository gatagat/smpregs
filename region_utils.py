import numpy as np
from collections import namedtuple, Counter
from interval_linked_list import IntervalLinkedList
from pyfasta import Fasta
import logging
from kmers import all_kmers, count_kmers

Region = namedtuple('Region', ['chrom', 'start', 'stop', 'name'])


def get_log(name, level=None):
    logger = logging.getLogger(name)
    if level is not None:
        logger.setLevel(level)
    return logger


def regions_reader(*args):
    """
    Provide generator access to regions in one or more BED files.
    """
    for filename in args:
        with open(filename, 'r') as f:
            for line in f:
                line = line.rstrip('\r\n')
                toks = line.split('\t', 4)
                if len(toks) < 3:
                    RuntimeError('At least 3 columns expected in input %s. Only %d found on line %s.' % \
                            (filename, len(toks), line))
                if len(toks) == 3:
                    yield Region(chrom=toks[0], start=int(toks[1]), stop=int(toks[2]), name=None)
                else:
                    yield Region(chrom=toks[0], start=int(toks[1]), stop=int(toks[2]), name=toks[3])


class RegionAcceptor(object):
    def __init__(self, template=None, fasta=None, **kwargs):
        self.template = template
        self.fasta = fasta

    def accept(self, region):
        raise NotImplemented('Subclasses of RegionAcceptor have to implement this method.')

    @property
    def reason(self):
        if self._reason_args == True:
            return 'accepted'
        else:
            return self._reason_args[0] % self._reason_args[1:]


def count_g_and_c(seq):
    """
    Count the number of C/G in the string.
    """
    #TODO: cythonize
    return (np.in1d(list(seq), ['c', 'g', 'C', 'G'])).sum()


class RegionAcceptorApproxGC(RegionAcceptor):
    """
    Acceptor of regions depending on the GC-content.
    """

    def __init__(self, threshold=10, **kwargs):
        assert threshold >= 0
        super(RegionAcceptorApproxGC, self).__init__(**kwargs)
        seq = self.fasta[self.template.chrom][self.template.start:self.template.stop]
        self.gc = count_g_and_c(seq)
        if threshold <= 1.:
            self.threshold = threshold * len(seq)
        else:
            self.threshold = threshold

    def accept(self, region):
        seq = self.fasta[region.chrom][region.start:region.stop]
        gc = count_g_and_c(seq)
        diff = abs(self.gc - gc)
        if diff <= self.threshold:
            self._reason_args = True
            return True
        else:
            self._reason_args = ('difference %d', diff)
            return False

class RegionAcceptorNoNs(RegionAcceptor):
    """
    Acceptor of regions requiring no unknown (N) nucleotides.
    """

    def __init__(self, **kwargs):
        super(RegionAcceptorNoNs, self).__init__(**kwargs)

    def accept(self, region):
        seq = self.fasta[region.chrom][region.start:region.stop]
        if not 'N' in seq and 'n' not in seq:
            self._reason_args = True
            return True
        else:
            self._reason_args = ('N in sequence', )
            return False


class RegionAcceptorAND(RegionAcceptor):
    """
    Meta region acceptor allowing logical AND of a given set of acceptors.
    """
    def __init__(self, acceptors=None, **kwargs):
        super(RegionAcceptorAND, self).__init__(**kwargs)
        assert len(acceptors) > 0
        self.acceptors = acceptors

    def accept(self, region):
        for a in self.acceptors:
            if not a.accept(region):
                return False
        return True


class GenomicAnnotationsHistogram(object):
    """
    Wraps computation of histograms over genomic annotations.

    Genomic annotations for the whole genome have to be encoded in fasta
    format. Use encode_annotations.py to create it from a BED file.
    """
    def __init__(self, filename):
        """
        Instantiate GenomicAnnotationsHistogram object.

        Parameters:
        ===========
        filename: string
            FASTA file with genomic annotations.
        """
        self.regions_fa = Fasta(filename)

    def __call__(self, region):
        """
        Compute the genomic annotation histogram for a given region.
        """
        return dict(Counter(
            self.regions_fa[region.chrom][region.start:region.stop]))


def histogram_intersection(p, q):
    ret = 0
    for k in set(p.keys()).intersection(q.keys()):
        ret += min(p[k], q[k])
    return ret


class RegionAcceptorApproxHistogram(RegionAcceptor):
    def __init__(self, histogram=None, dissimilarity=None, threshold=None, features_per_nt=None, **kwargs):
        """
        histogram: callable
            Function returning a dict with histogram for a given region.
        dissimilarity: callable
            Function returning a dissimilarity of histograms of template and
            candidate regions. If None, complement of histogram intersection is taken.
        threshold: float
            Upper bound on accepted dissimilarity.
        """
        super(RegionAcceptorApproxHistogram, self).__init__(**kwargs)
        if dissimilarity is None:
            dissimilarity=lambda p, q: \
                features_per_nt*(self.template.stop - self.template.start) \
                - histogram_intersection(p, q)
        self.histogram = histogram
        self.dissimilarity = dissimilarity
        self.threshold = threshold
        self.template_hist = self.histogram(self.template)
        if self.template_hist is None:
            raise ValueError, 'Cannot generate a matching region for %s' % self.template

    def accept(self, region):
        hist = self.histogram(region)
        if hist is None:
            return False
        dis = self.dissimilarity(self.template_hist, hist)
        if dis < self.threshold:
            self._reason_args = True
            return True
        else:
            self._reason_args = ('dissimilarity %d', dis)
            return False


class GenomicAnnotationsAtPosition(object):
    def __init__(self, filename):
        self.regions_fa = Fasta(filename)

    def __call__(self, chrom, position):
        return self.regions_fa[chrom][position]


class RegionAcceptorGenomicAnnotation(RegionAcceptor):
    def __init__(self, filename=None, pos=None, **kwargs):
        """
        filename: string
            File with genomic annotations for the whole genome encoded in fasta format.
            (Use encode_annotations.py to create it.)
        pos: int
            Take into account only a single position (eg. peak summit, 0 == 1st bp)

        Exactly one option of either position or dissimilarity-and-threshold has to be specified.
        """
        super(RegionAcceptorGenomicAnnotation, self).__init__(**kwargs)
        self.ga = GenomicAnnotationsAtPosition(filename)
        self.position = int(pos)
        self.template_ga = self.ga(
                self.template.chrom, self.template.start + self.position)
        assert pos >= 0 and pos < (self.template.stop - self.template.start)

    def accept(self, region):
        ga = self.ga(region.chrom, region.start + self.position)
        if self.template_ga == ga:
            self._reason_args = True
            return True
        else:
            self._reason_args = ('%s != %s', ga, self.template_ga)
            return False


class KmerHistogram(object):
    """
    Wraps computation of histograms of k-mers.
    """
    def __init__(self, fasta, k=2):
        self.fasta = fasta
        self.k = k
        self.keys = all_kmers(self.k)

    def __call__(self, region):
        """
        Compute the k-mer histogram for a given region.
        """
        seq = self.fasta[region.chrom][region.start:region.stop]
        if not np.in1d(list(str(seq).upper()), list('ACGT')).all():
            return None
        return dict(zip(self.keys, count_kmers(self.k, str(seq).upper())))


class RegionAcceptorFeatureCount(RegionAcceptor):
    def __init__(self, filename=None, threshold=None, **kwargs):
        raise NotImplemented
        #super(RegionAcceptorFeatureCount, self).__init__(**kwargs)
        #assert threshold > 0
        #self.threshold = threshold
        #cnt = self.template
        #self.lb = None
        #self.ub = None

    def accept(self, region):
        #cnt = ...
        #return (cnt > self.lb) and (cnt < self.ub)
        raise NotImplemented


class AllowedSpace(object):
    """
    Represent remaining available space where new regions are allowed.
    """

    def __init__(self, fasta=None, include=None, exclude=None):
        """
        fasta - pyfasta.Fasta object
        include - iterable of Region-s
        exclude - iterable of Region-s
        """
        logger = get_log('AllowedSpace')
        if include is None and fasta is None:
            raise ValueError('Either include or fasta have to be specified.')
        self._range = {}
        self._space = {}
        if include is None:
            for k in fasta.keys():
                self._space[k] = IntervalLinkedList([(0, len(fasta[k]))])
        else:
            intervals = {}
            for region in include:
                if region.start >= region.stop:
                    raise ValueError('Region %r is invalid (start >= stop).' % region)
                if region.chrom not in intervals:
                    intervals[region.chrom] = []
                intervals[region.chrom] += [(region.start, region.stop)]
            for k, v in intervals.items():
                pos = v[0][0]
                is_sorted = True
                for start, stop in v:
                    is_sorted = is_sorted and start >= pos
                    pos = stop
                if is_sorted:
                    self._space[k] = IntervalLinkedList(v)
                else:
                    logger.warn('Sorting %d include regions for chromosome %s.' % (len(intervals), k))
                    self._space[k] = IntervalLinkedList(sorted(v))
        if exclude is not None:
            for region in exclude:
                self._space[region.chrom].remove((region.start, region.stop))
        for k in self._space.keys():
            self._update_range(k)


    def remove(self, region):
        """
        Remove region from the allowed space.
        """
        self._space[region.chrom].remove((region.start, region.stop))
        self._update_range(region.chrom)


    def _update_range(self, chrom):
        current = self._space[chrom].next
        start = current.data[0]
        while current is not None:
            stop = current.data[1]
            current = current.next
        self._range[chrom] = (start, stop)


    def contains(self, region):
        """
        Check whether the given region is fully inside this space.
        """
        return region.chrom in self._space and (region.start, region.stop) in self._space[region.chrom]


    def range(self, chrom):
        return self._range[chrom]


def the_random_regions_lair(region, lo, hi):
    """
    Infinite generator of a random regions in the interval [lo, hi).

    Output regions will be on the same chromosome and will have the same
    length as the given region.
    """
    region_len = region.stop - region.start
    while True:
        start = np.random.randint(lo, hi)
        stop = start + region_len
        yield Region(chrom=region.chrom, start=start, stop=stop, name='rnd_' + region.name)


def generate(input_region, allowed_space, max_generate_iter=10000):
    """
    Generate a random region for the given region and in the allowed space.
    """
    logger = get_log('generate')
    i = 0
    lo, hi = allowed_space.range(input_region.chrom)
    for region in the_random_regions_lair(input_region, lo, hi):
        if allowed_space.contains(region):
            logger.debug('GEN %s', region)
            return region
        logger.info('NOT %s', region)
        i += 1
        if i > max_generate_iter:
            raise RuntimeError('Failed to generate a non-overlapping region.')


