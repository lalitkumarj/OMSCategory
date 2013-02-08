# cython: profile=True

# This has been placed in "working_copy" simply to not break the sagelib
# For now we can add all cython code into here

#*****************************************************************************
#       Copyright (C) 2012 Robert Pollack <rpollack@math.bu.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.arith import binomial, bernoulli
from sage.rings.padics.padic_fixed_mod_element cimport pAdicFixedModElement
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.finite_rings.integer_mod cimport IntegerMod_gmp
from sage.misc.misc import verbose, cputime
from sage.structure.sequence import Sequence
from sage.rings.padics.factory import Qp, Zp
from sage.misc.misc_c import prod
from sage.functions.other import factorial
from sage.matrix.constructor import Matrix
from sage.matrix.matrix cimport Matrix

include "stdsage.pxi"
include "cdefs.pxi"

cdef extern from "zn_poly/zn_poly.h":
    pass

#from sage.libs.flint.zmod_poly cimport *, zmod_poly_t
#from sage.libs.flint.long_extras cimport *

cdef long overflow = 1 << (4*sizeof(long)-1)
cdef long underflow = -overflow
cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) - 1


###############################################################################


def ps_normalize(f, Integer p, int p_prec):
    """reduces all of the coefficients of the power series modulo p^N"""
    v = Sequence(f)
    cdef Integer P = p ** p_prec
    cdef int a
    cdef int Len = len(v)
    v = [(v[a]%P) for a in range(Len)]
    S = f.parent()
    return S(v)

def ps_normalize_long(f, long p, int p_prec):
    """reduces all of the coefficients of the power series modulo p^N"""
    v = Sequence(f)
    cdef long P = p ** p_prec
    cdef int a
    cdef int Len = len(v)
    v = [(v[a]%P) for a in range(Len)]
    S = f.parent()
    return S(v)

def logp_fcn(Integer p, int p_prec, z):
    """this is the *function* on Z_p^* which sends z to log_p(z) using a power series truncated at p_prec terms"""
    R = Zp(p, 2 * p_prec)
    cdef int m
    z = z / R.teichmuller(z)
    return sum([( (-1) * (1 - z) ** m ) / m for m in range(1, p_prec)])

def logp_fcn_long(long p, int p_prec, z):
    """this is the *function* on Z_p^* which sends z to log_p(z) using a power series truncated at p_prec terms"""
    raise NotImplementedError

#Don't know how to cdef any better than this. Rational coefficients => ???
def logpp(Integer p, int p_prec):
    """returns the (integral) power series for log_p(1+p*z) -- extra p here!"""
    cdef int m
    SS = PolynomialRing(QQ, 'y')
    y = SS.gen()
    return sum([((-1) ** (m - 1)) * ((p * y) ** m) / m for m in range(1, p_prec)])

def logpp_long(long p, int p_prec):
    """returns the (integral) power series for log_p(1+p*z) -- extra p here!"""
    raise NotImplementedError

def logpp_gam(Integer p, int p_prec):
    """returns the (integral) power series log_p(1+p*z)*(1/log_p(1+p)) where the denominator is computed with some accuracy"""
    L = logpp(p, p_prec)
    loggam = ZZ(logp_fcn(p, p_prec * (p ** 2), 1 + p))
    return ps_normalize(L / loggam, p, p_prec)

def logpp_gam_long(long p, int p_prec):
    """returns the (integral) power series log_p(1+p*z)*(1/log_p(1+p)) where the denominator is computed with some accuracy"""
    raise NotImplementedError

#@cached_function
def logpp_binom(int n, Integer p, int p_prec):
    """returns the (integral) power series p^n*(log_p(1+p*z)/log_p(1+p) choose n)"""
    cdef int j
    L = logpp_gam(p, p_prec)
    ans = prod([(L - j) for j in range(n)])
    ans *= (p ** n) / factorial(n)
    return ps_normalize(ans.truncate(p_prec), p, p_prec)

def logpp_binom_long(int n, long p, int p_prec):
    """returns the (integral) power series p^n*(log_p(1+p*z)/log_p(1+p) choose n)"""
    raise NotImplementedError

def automorphy_factor_vector(Integer p, a, c, k, chi, int p_prec, int var_prec, R):
    cdef int n
    S = PolynomialRing(R, 'z')
    z = S.gens()[0]
    w = R.gen()
    aut = S(1)
    for n in range(1, var_prec):
        LB = logpp_binom(n, p, p_prec)
        ta = ZZ(Qp(p, 2 * max(p_prec, var_prec)).teichmuller(a))
        arg = (a / ta - 1) / p + c / (p * ta) * z
        aut += LB(arg) * (w ** n)
    aut *= (ta ** k)
    if not (chi is None):
        aut *= chi(a)    
    return Sequence(aut)

def automorphy_factor_vector_long(long p, a, c, k, chi, int p_prec, int var_prec, R):
    raise NotImplementedError

#@cached_function
def automorphy_factor_matrix(Integer p, a, c, k, chi, int p_prec, int var_prec, R):
    cdef int i,j
    cdef int Len
    cdef int min1, min2
    aut = automorphy_factor_vector(p, a, c, k, chi, p_prec, var_prec, R)
    Len = len(aut)
    M = Matrix(R, p_prec)
    min1 = min(p_prec, Len)
    for i in range(min1):
        min2 = min(p_prec, Len - i)
        for j in range(min2):
            M[j + i, i] = aut[j]
    return M

def automorphy_factor_matrix(long p, a, c, k, chi, int p_prec, int var_prec, R):
    raise NotImplementedError

