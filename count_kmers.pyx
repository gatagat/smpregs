# Tomas Kazmar, 2012, BSD 2-clause license, see LICENSE

import numpy as np
cimport numpy as np
cimport cython
from libc.stdlib cimport malloc, free
from libc.stdint cimport int32_t

def all_kmers(k, i=None):
    if i is None:
        i = k-1
    letters = 'ACGT'
    if i > 0:
        kmers = []
        for c in letters:
            kmers += [s + c for s in all_kmers(k, i-1)]
        return kmers
    else:
        return [c for c in letters]


def count_kmers(k, seq=None):
    if seq is None:
        return all_kmers(k)
    map_to_index = { 'A': 3, 'C': 2, 'G': 1, 'T': 0 }
    cdef np.ndarray[np.uint32_t, ndim=1] seqi = np.array([map_to_index[c] for c in seq], dtype=np.uint32)
    #cdef np.ndarray[np.uint32_t, ndim=1] seqi = np.array([
    #    3 if c == 'A' else 2 if c == 'C' else 1 if c == 'G' else 0 for c in seq], dtype=np.uint32)
    #cdef np.ndarray[np.uint32_t, ndim=1] map_to_index = np.array([3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.uint32)
    #cdef np.ndarray[np.uint32_t, ndim=1] seqi = np.array([map_to_index[ord(c) - ord('A')] for c in seq], dtype=np.uint32)
    return _count_kmers(k, seqi)

@cython.boundscheck(False)
def _count_kmers(unsigned int k, np.ndarray[np.uint32_t, ndim=1] seqi):
    cdef np.ndarray[np.uint32_t, ndim=1] counts = np.zeros(4**k, dtype=np.uint32)
    cdef np.uint32_t index_mask = 4 ** k - 1
    cdef np.uint32_t index
    cdef np.uint32_t c
    index = 0
    for c in seqi[:-1+<int>k]:     # initialize index
        index = index*4 + c
    for c in seqi[-1+<int>k:]:     # count for the opposite strand
        index = ((index << 2) & index_mask) + c
        counts[index] += 1
    seqi = 3 - seqi                # take complement (back to original)
    index = 0
    for c in seqi[:-<int>k:-1]:    # initialize index for the fwd strand
        index = index*4 + c
    for c in seqi[-<int>k::-1]:    # count
        index = ((index << 2) & index_mask) + c
        counts[index] += 1
    return counts


def create_kmers_lookup(k, seq):
    map_to_index = { 'A': 3, 'C': 2, 'G': 1, 'T': 0 }
    cdef np.ndarray[np.uint32_t, ndim=1] seqi = np.array([map_to_index[c] for c in seq], dtype=np.uint32)
    return _create_kmers_lookup(k, seqi)


@cython.boundscheck(False)
def _create_kmers_lookup(unsigned int k, np.ndarray[np.uint32_t, ndim=1] seqi):
    cdef np.ndarray[np.uint16_t, ndim=2] counts = np.zeros((len(seqi)+2-<int>k, 4**k), dtype=np.uint16)
    cdef np.uint32_t index_mask = 4 ** k - 1
    cdef np.uint32_t index
    cdef np.uint32_t c
    index = 0
    for c in seqi[:-1+<int>k]:     # initialize index
        index = index*4 + c
    i = 1
    for c in seqi[-1+<int>k:]:     # count for the opposite strand
        index = ((index << 2) & index_mask) + c
        counts[i, index] = 1
        i += 1
    seqi = 3 - seqi                # take complement (back to original)
    index = 0
    for c in seqi[:-<int>k:-1]:    # initialize index for the fwd strand
        index = index*4 + c
    i = len(seqi) + 1 - <int>k
    for c in seqi[-<int>k::-1]:    # count
        index = ((index << 2) & index_mask) + c
        counts[i, index] += 1
        i -= 1
    return np.cumsum(counts, axis=0, dtype=counts.dtype)


def count_kmers_lookup(np.ndarray[np.uint16_t, ndim=2] lookup, start, stop):
    k = int(np.round(np.log(lookup.shape[1]) / np.log(4.)))
    return lookup[stop - k + 1, :] - lookup[start, :]
