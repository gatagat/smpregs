__all__ = ['regions_for', 'generate']

import os
from .smpregs import get_assembly, parse_filters, sample_regions
from .region_utils import AllowedSpace, RegionAcceptorNoNs, generate

def regions_for(regions, include_file=None, ctrl_props=None):
    """
    Standalone control region generator.

    Simple interface to produce control regions from Python code without any disk access.
    """
    if ctrl_props is None:
        ctrl_props = []
    genome_fasta = get_assembly('dm3')
    annotations_file = os.path.expanduser('~/projs/smpregs/data/genomic-annotations-dm3.fa')
    acceptors = parse_filters(ctrl_props, genome_fasta, annotations_file)
    acceptors = [(RegionAcceptorNoNs, {})] + acceptors
    allowed_space = AllowedSpace(fasta=genome_fasta, include=include_file)
    for _, ctrl in sample_regions(regions, allowed_space, acceptors, genome_fasta):
        yield ctrl


