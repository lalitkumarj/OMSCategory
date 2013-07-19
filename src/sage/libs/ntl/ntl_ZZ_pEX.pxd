include "decl.pxi"
include "sage/ext/cdefs.pxi"

from sage.libs.ntl.ntl_ZZ_pEContext cimport ntl_ZZ_pEContext_class

cdef class ntl_ZZ_pEX:
    cdef ZZ_pEX_c x
    cdef ntl_ZZ_pEContext_class c
    #cdef void setitem_from_int(ntl_ZZ_pX self, long i, int value)
    #cdef int getitem_as_int(ntl_ZZ_pX self, long i)
    cdef ntl_ZZ_pEX _new(self)