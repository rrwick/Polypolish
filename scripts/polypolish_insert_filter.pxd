# This file contains Cython definitions for the polypolish_insert_filter.py script. It is only
# used when using Cython to create a compiled version of the script.

# See:
#   https://github.com/rrwick/Polypolish/wiki/Insert-size-filter
#   https://cython.readthedocs.io/en/latest/src/tutorial/pure.html


cdef get_orientation(Alignment alignment_1, Alignment alignment_2)

cdef int get_insert_size(alignment_1, alignment_2)

cdef int filter_sam(str in_filename, str out_filename, alignments, int low, int high, str correct_orientation, int read_num)

cdef bint alignment_pass_qc(Alignment a, this_alignments, pair_alignments, int low, int high, str correct_orientation)

cdef class Alignment:
    cdef public int sam_flags, ref_start, ref_end
    cdef public str read_name, ref_name
    cdef bint is_aligned(self)
    cdef bint is_on_forward_strand(self)
    cdef bint has_flag(self, int flag)

cdef int get_ref_end(int ref_start, str cigar)
