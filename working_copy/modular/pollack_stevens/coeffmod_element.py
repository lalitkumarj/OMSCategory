#This corresponds to dist.pyx and hence should be cythoned.
# cython: profile=True

#*****************************************************************************
#       Copyright (C) 2012 Robert Pollack <rpollack@math.bu.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
#from sage.rings.rational_field import QQ
#from sage.rings.polynomial.all import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
#from sage.rings.finite_rings.integer_mod_ring import Zmod
#from sage.rings.arith import binomial, bernoulli
#from sage.modules.free_module_element import vector, zero_vector
#from sage.matrix.matrix cimport Matrix
#from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.all import matrix
#from sage.misc.prandom import random
#from sage.functions.other import floor
#from sage.structure.element cimport RingElement, Element
import operator
##from sage.modular.overconvergent.pollack.S0p import S0
#from sage.rings.padics.padic_generic import pAdicGeneric
#from sage.rings.padics.padic_capped_absolute_element cimport pAdicCappedAbsoluteElement
#from sage.rings.padics.padic_capped_relative_element cimport pAdicCappedRelativeElement
#from sage.rings.padics.padic_fixed_mod_element cimport pAdicFixedModElement
#from sage.rings.integer cimport Integer
#from sage.rings.rational cimport Rational
#from sage.misc.misc import verbose, cputime

#cdef extern from "zn_poly/zn_poly.h":
#    pass
#from sage.libs.flint.zmod_poly cimport *, zmod_poly_t
#from sage.libs.flint.long_extras cimport *

from sage.modular.pollack_stevens.sigma0 import Sigma0

#cdef long overflow = 1 << (4*sizeof(long)-1)
#cdef long underflow = -overflow
#cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) - 1

#include "stdsage.pxi"
#include "cdefs.pxi"

from sage.categories.action import Action
from sage.structure.element import ModuleElement
from sage.modular.pollack_stevens.families_util import automorphy_factor_matrix

class WeightKAction_generic(Action):
    def __init__(self, Dk, character, adjuster, on_left, dettwist):
        self._k = Dk._k
        self._adjuster = adjuster
        self._character = character
        self._dettwist = dettwist
        #self._p = Dk._p
        self._actmat = {}
        self._maxprecs = {}
        if character is None:
            self._Np = ZZ(1) # all of M2Z acts
        else:
            self._Np = character.modulus()
        #Do something about this in a derived class
        #if not self._symk:
        #    self._Np = self._Np.lcm(self._p)
        #if padic:
        #    try:
        #        self._p = Dk._p
        #        self._Np = self._Np.lcm(self._p)
        #    except AttributeError:
        #        pass
        #    Action.__init__(self, Sigma0(self._Np, base_ring=Dk.base_ring(), adjuster=self._adjuster), Dk, on_left, operator.mul)
        #else:
        #    Action.__init__(self, Sigma0(self._Np, base_ring=ZZ, adjuster=self._adjuster), Dk, on_left, operator.mul)
    
    def clear_cache(self):
        r"""
            
        """
        self._actmat = {}
        self._maxprecs = {}
    
    #cpdef acting_matrix(self, g, M):
    def acting_matrix(self, g, M):
        r"""
        

        INPUT:

        - ``g`` -- an instance of
          :class:`sage.matrices.matrix_integer_2x2.Matrix_integer_2x2`
          or of :class:`sage.matrix.matrix_generic_dense.Matrix_generic_dense`

        - ``M`` -- a positive integer giving the precision at which
          ``g`` should act.

        OUTPUT:

        - An `M \times M` matrix so that the action of `g` on a
          distribution with `M` moments is given by a vector-matrix
          multiplication.

        .. NOTE::

            This function caches its results.  To clear the cache use
            :meth:`clear_cache`.
        """
        g = g.matrix()
        if not self._maxprecs.has_key(g):
            A = self._compute_acting_matrix(g, M)
            self._actmat[g] = {M:A}
            self._maxprecs[g] = M
            return A
        else:
            mats = self._actmat[g]
            if mats.has_key(M):
                return mats[M]
            maxprec = self._maxprecs[g]
            if M < maxprec:
                A = mats[maxprec][:M,:M] # submatrix; might want to reduce precisions
                mats[M] = A
                return A
            if M < 2*maxprec:
                maxprec = 2*maxprec
            else:
                maxprec = M
            self._maxprecs[g] = maxprec
            mats[maxprec] = self._compute_acting_matrix(g, maxprec) # could lift from current maxprec
            if M == maxprec:
                return mats[maxprec]
            A = mats[maxprec][:M,:M] # submatrix; might want to reduce precisions
            mats[M] = A
            return A
        
    #cpdef _compute_acting_matrix(self, g, M):
    def _compute_acting_matrix(self, g, M):
        r"""
        

        INPUT:

        - ``g`` -- an instance of
          :class:`sage.matrices.matrix_integer_2x2.Matrix_integer_2x2`

        - ``M`` -- a positive integer giving the precision at which
          ``g`` should act.

        OUTPUT:

        - 

        Forms a large M x M matrix say G such that if v is the vector of
        moments of a distribution mu, then v*G is the vector of moments of
        mu|[a,b;c,d]
        """
        raise NotImplementedError

#This class should eventually be moved into its own file
class WeightKAction_OMS_fam(WeightKAction_generic):
    """
    """
    def __init__(self, Dk, character, adjuster, on_left, dettwist, padic=True):
        #Only difference is that it adds a cache for automorphy factors and
        #ensures there's a p in the level.
        #self._Np = self._Np.lcm(self._p)
        self._autfactors = {}
        WeightKAction_generic.__init__(self, Dk, character, adjuster, on_left, dettwist)
        self._Np = self._Np.lcm(Dk._p)
        Action.__init__(self, Sigma0(self._Np, base_ring=Dk.base_ring().base_ring(), \
                        adjuster=self._adjuster), Dk, on_left, operator.mul)
    
    def clear_cache(self):
        #Only difference is that it clears the cache for automorphy factors.
        self._actmat = {}
        self._maxprecs = {}
        self._autfactors = {}
    
    def _compute_aut_factor_matrix(self, g, M):
        #compute the power series
        D = self.domain()
        p_prec, var_prec = D.precision_cap()
        a, b, c, d = self._adjuster(g._mat)
        #print D.prime(), a, c, self._k, self._character, M, var_prec, D.base_ring()
        return automorphy_factor_matrix(D.prime(), a, c, self._k, \
                                self._character, M, var_prec, D.base_ring())
    
    #TODO: much of this should be moved to a separte utility file
    #cpdef _compute_acting_matrix(self, g, M):
    def _compute_acting_matrix(self, g, M):
        r"""
        

        INPUT:

        - ``g`` -- an instance of
          :class:`sage.matrices.matrix_integer_2x2.Matrix_integer_2x2`
          or :class:`sage.matrix.matrix_generic_dense.Matrix_generic_dense`

        - ``M`` -- a positive integer giving the precision at which
          ``g`` should act.

        OUTPUT:

        - 
        """
        #tim = verbose("Starting")
        a, b, c, d = self._adjuster(g)
        # if g.parent().base_ring().is_exact():
        #     self._check_mat(a, b, c, d)
        k = self._k
        if g.parent().base_ring() is ZZ:
            if self._symk:
                base_ring = QQ
            else:
                base_ring = Zmod(self._p**M)
        else:
            base_ring = self.domain().base_ring()
        #cdef Matrix B = matrix(base_ring,M,M)
        B = matrix(base_ring,M,M) #
        if M == 0:
            return B.change_ring(self.codomain().base_ring())
        R = PowerSeriesRing(base_ring, 'y', default_prec = M)
        y = R.gen()
        #tim = verbose("Checked, made R",tim)
        # special case for small precision, large weight
        scale = (b+d*y)/(a+c*y)
        t = (a+c*y)**k # will already have precision M
        #cdef long row, col #
        #tim = verbose("Made matrix",tim)
        for col in range(M):
            for row in range(M):
                #B.set_unsafe(row, col, t[row])
                B[row, col] = t[row]
            t *= scale
        #verbose("Finished loop",tim)
        # the changering here is annoying, but otherwise we have to change ring each time we multiply
        B = B.change_ring(self.codomain().base_ring())
        if self._character is not None:
            B *= self._character(a)
        if self._dettwist is not None:
            B *= (a*d - b*c)**(self._dettwist)
        return B
    
    def get_action_matrices(self, g, M):
        #g = M2Z(g)
        #g.set_immutable()
        if not self._maxprecs.has_key(g):
            AF = self._compute_aut_factor_matrix(g, M)
            self._autfactors[g] = {M : AF}
        else:
            auts = self._autfactors[g]
            if auts.has_key(M):
                AF = auts[M]
            else:
                maxprec = self._maxprecs[g]
                if M < maxprec:
                    AF = auts[maxprec][:M,:M]
                    auts[M] = AF
                else:
                    if M < 2 * maxprec:
                        maxprec = 2 * maxprec
                    else:
                        maxprec = M
                    auts[maxprec] = self._compute_aut_factor_matrix(g, maxprec)
                    if M == maxprec:
                        AF = auts[maxprec]
                    else:
                        AF = auts[maxprec][:M,:M]
                        auts[M] = AF
        
        A = self.acting_matrix(g, M)
        return [AF, A]
    
    def _call_(self, v, g):
        AF, A = self.get_action_matrices(g, len(v._moments))
        w = (v._moments * AF) * A
        Elem = v.parent().Element
        return Elem(w, v.parent(), ordp=v.ordp, check=False)

#This class should eventually be moved into its own file
class WeightKAction_OMS(WeightKAction_generic):
    """
    """
    def __init__(self, Dk, character, adjuster, on_left, dettwist, padic=True):
        #ensures there's a p in the level.
        #self._Np = self._Np.lcm(self._p)
        self._autfactors = {}
        WeightKAction_generic.__init__(self, Dk, character, adjuster, on_left, dettwist)
        self._Np = self._Np.lcm(Dk._p)
        Action.__init__(self, Sigma0(self._Np, base_ring=Dk.base_ring().base_ring(), \
                        adjuster=self._adjuster), Dk, on_left, operator.mul)
    
    def clear_cache(self):
        self._actmat = {}
        self._maxprecs = {}
    
    #TODO: much of this should be moved to a separte utility file
    #cpdef _compute_acting_matrix(self, g, M):
    def _compute_acting_matrix(self, g, M):
        r"""
        

        INPUT:

        - ``g`` -- an instance of
          :class:`sage.matrices.matrix_integer_2x2.Matrix_integer_2x2`
          or :class:`sage.matrix.matrix_generic_dense.Matrix_generic_dense`

        - ``M`` -- a positive integer giving the precision at which
          ``g`` should act.

        OUTPUT:

        - 
        """
        #tim = verbose("Starting")
        a, b, c, d = self._adjuster(g)
        # if g.parent().base_ring().is_exact():
        #     self._check_mat(a, b, c, d)
        k = self._k
        if g.parent().base_ring() is ZZ:
            if self._symk:
                base_ring = QQ
            else:
                base_ring = Zmod(self._p**M)
        else:
            base_ring = self.domain().base_ring()
        #cdef Matrix B = matrix(base_ring,M,M)
        B = matrix(base_ring,M,M) #
        if M == 0:
            return B.change_ring(self.codomain().base_ring())
        R = PowerSeriesRing(base_ring, 'y', default_prec = M)
        y = R.gen()
        #tim = verbose("Checked, made R",tim)
        # special case for small precision, large weight
        scale = (b+d*y)/(a+c*y)
        t = (a+c*y)**k # will already have precision M
        #cdef long row, col #
        #tim = verbose("Made matrix",tim)
        for col in range(M):
            for row in range(M):
                #B.set_unsafe(row, col, t[row])
                B[row, col] = t[row]
            t *= scale
        #verbose("Finished loop",tim)
        # the changering here is annoying, but otherwise we have to change ring each time we multiply
        B = B.change_ring(self.codomain().base_ring())
        if self._character is not None:
            B *= self._character(ZZ(a))
        if self._dettwist is not None:
            B *= (a*d - b*c)**(self._dettwist)
        return B
    
    def acting_matrix(self, g, M):
        g = g.matrix()
        if not self._maxprecs.has_key(g):
            A = self._compute_acting_matrix(g, M)
            self._actmat[g] = {M:A}
            self._maxprecs[g] = M
            return A
        else:
            mats = self._actmat[g]
            if mats.has_key(M):
                return mats[M]
            maxprec = self._maxprecs[g]
            if M < maxprec:
                A = mats[maxprec][:M,:M] # submatrix; might want to reduce precisions
                mats[M] = A
                return A
            if M < 2*maxprec:
                maxprec = 2*maxprec
            else:
                maxprec = M
            self._maxprecs[g] = maxprec
            mats[maxprec] = self._compute_acting_matrix(g, maxprec) # could lift from current maxprec
            if M == maxprec:
                return mats[maxprec]
            A = mats[maxprec][:M,:M] # submatrix; might want to reduce precisions
            mats[M] = A
            return A
    
    def _call_(self, v, g):
        r"""
            EXAMPLES::
            
                sage: D = OverconvergentDistributions(0,7,base=Qp(7,10))
                sage: mu = D([1,1,1])
                sage: mu
                (1 + O(7^3), 1 + O(7^2), 1 + O(7))
                sage: mu = mu / 7
                sage: mu
                7^-1 * (1 + O(7^3), 1 + O(7^2), 1 + O(7))
                sage: S0 = D.action().actor()
                sage: A = S0([1,1,0,1])
                sage: mu * A
                7^-1 * (1 + O(7^3), 2 + O(7^2), 4 + O(7))
                sage: A = S0([1,0,0,1])
                sage: mu * A
                7^-1 * (1 + O(7^3), 1 + O(7^2), 1 + O(7))
        """
        #Could change this to look like corresponding function for families?
        A = self.acting_matrix(g, len(v._moments))
        ans = v.parent()(v._moments * A)
        #print "v =", v
        #print "ans = ", v._moments * A
        ans.ordp += v.ordp   #may be redundant
        return ans

class CoefficientModuleElement_generic(ModuleElement):
    pass
