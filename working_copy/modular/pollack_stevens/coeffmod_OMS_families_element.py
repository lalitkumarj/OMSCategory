# This should be cythoned once it's done.

from copy import copy   #Not necessary in cython version
from sage.rings.infinity import Infinity
from sage.modular.pollack_stevens.coeffmod_element import CoefficientModuleElement_generic
from sage.modular.pollack_stevens.coeffmod_OMS_element import CoeffMod_OMS_element
from sage.rings.integer_ring import ZZ
from sage.functions.other import ceil

class CoeffMod_OMS_Families_element(CoefficientModuleElement_generic):
    # Implementation currently ignores ordp
    def __init__(self, moments, parent, ordp=0, check=True, var_prec=None):
        CoefficientModuleElement_generic.__init__(self, parent)
        #TODO: finish this
        p = parent.prime()
        if var_prec == None:
            var_prec = parent._prec_cap[1]
        self._var_prec = var_prec
        if check:
            if isinstance(moments, CoeffMod_OMS_Families_element):
                ordp = moments.ordp
                moments = moments._moments.change_ring(parent.base_ring())
            if isinstance(moments, CoeffMod_OMS_element):
                ordp = moments.ordp
                moments = self.parent().approx_module(p_prec=len(moments._moments), var_prec=self.parent().precision_cap()[1])(moments._moments)
            elif hasattr(moments, '__len__'):
                #Need to modify if uniformiser is not p
                #Deal with var_prec
                R = self.parent().base_ring()
                moments = [R(x) for x in moments]
                ordp = min(map(lambda x : _padic_val_of_pow_series(x, p=p), moments))
                if ordp == Infinity:
                    ordp = parent._prec_cap[0]
                    moments = parent.approx_module(0)([])
                else:
                    M = len(moments)
                    moments = parent.approx_module(M)(moments)
                    if ordp != 0:
                        moments *= p ** (-ordp)
            elif moments == 0:
                ordp = parent._prec_cap[0]
                moments = parent.approx_module(parent._prec_cap[0], parent._prec_cap[1])(moments)
            else:
                moments = parent.approx_module(1, parent.precision_cap()[1])([moments])
                ordp = moments[0].valuation()
                if ordp == Infinity:
                    ordp = parent._prec_cap[0]
                    moments = parent.approx_module(0)([])
                elif ordp != 0:
                    moments *= p ** (-ordp)
                
        self._moments = moments
        #if var_prec is None:
        #    self._var_prec = parent._prec_cap[1]
        #else:
        #    self._var_prec = var_prec
        self.ordp = ordp   #eventually change this maybe
        self._cp = (p-2) / (p-1)
    
    #def _relprec(self):
    #    return len(self._moments)
    def _repr_(self):
        #RH: adapted from coeffmod_OMS_element.py
        self.normalize()
        valstr = ""
        if self.ordp == 1:
            valstr = "%s * "%(self.parent().prime())
        elif self.ordp != 0:
            valstr = "%s^%s * "%(self.parent().prime(), self.ordp)
        if len(self._moments) == 1:
            return valstr + repr(self._moments[0])
        else:
            return valstr + repr(self._moments)
    
    def character(self):
        """
        Returns the tame character of self.
        """
        return self._character
    
    def _add_(self, right):
        #RH: copied from coeffmod_OMS_element.py
        return self._addsub(right, False)
    
    def _sub_(self, right):
        #RH: copied from coeffmod_OMS_element.py
        return self._addsub(right, True)
    
    def _addsub(self, right, negate):
        #RH: adapted from coeffmod_OMS_element.py
        ans = self.parent()(0)
        saprec = self.precision_absolute()
        raprec = right.precision_absolute()
        #print "In _addsub"
        #print saprec, raprec
        aprec = min(saprec[0], raprec[0])
        #print "aprec", aprec
        ans.ordp = min(self.ordp, right.ordp)
        #print self.ordp, right.ordp
        #print "ans.ordp", ans.ordp
        rprec = aprec - ans.ordp
        #print rprec
        #print self.parent().precision_cap()
        var_prec = min(saprec[1], raprec[1])
        V = ans.parent().approx_module(rprec, var_prec)
        R = V.base_ring()
        smoments = self._moments
        rmoments = right._moments
        if smoments.parent() is not V:
            smoments = V(smoments.list(copy=False)[:rprec] + ([R(0)] * (rprec - len(smoments)) if rprec > len(smoments) else []))
        if rmoments.parent() is not V:
            rmoments = V(rmoments.list(copy=False)[:rprec] + ([R(0)] * (rprec - len(rmoments)) if rprec > len(rmoments) else []))
        # We multiply by the relative power of p
        if self.ordp > right.ordp:
            smoments *= self.parent().prime()**(self.ordp - right.ordp)
        elif self.ordp < right.ordp:
            rmoments *= self.parent().prime()**(right.ordp - self.ordp)
        if negate:
            rmoments = -rmoments
        ans._moments = smoments + rmoments
        return ans
    
    def _lmul_(self, right):
        #RH: adapted from coeffmod_OMS_element.py
        """
        Scalar multiplication self*right.
        """
        ans = self.parent()(0)
        if right.is_zero():
            return ans
        p = self.parent().prime()
        if right.is_zero():
            ans._moments = self.parent().approx_module(0, 0)([])
            ans.ordp = self.parent().precision_cap()[0]
        else:
            v, u = _padic_val_unit_of_pow_series(right, p)
            ans._moments = self._moments * u
            ans.ordp = self.ordp + v
        return ans
    
    def _rmul_(self, left):
        #RH: copied from coeffmod_OMS_element.py
        """
        Scalar multiplication left*self.
        """
        return self._lmul_(left)
    
    def __cmp__(left, right):
        #RH: adapted from coeffmod_OMS_element.py
        left = copy(left)
        right = copy(right)
        left.normalize()
        right.normalize()
        lrprec = left.precision_relative()[0]
        rrprec = right.precision_relative()[0]
        if lrprec == 0:
            if rrprec == 0:
                return 0
            else:
                return -1
        elif rrprec == 0:
            return 1
        rprec = min(lrprec, rrprec)
        #if lrprec > rprec:
        #    left.reduce_precision(rprec)
        #elif rrprec > rprec:
        #    right.reduce_precision(rprec)
        #left.normalize()
        #right.normalize()
        p = left.parent().prime()
        if left.ordp > right.ordp:
            shift = p ** (left.ordp - right.ordp)
            for i in range(rprec):
                L = (shift * left._unscaled_moment(i)).padded_list()
                R = (right._unscaled_moment(i)).padded_list()
                c = cmp(L, R)
                if c: return c
        elif left.ordp < right.ordp:
            shift = p ** (right.ordp - left.ordp)
            for i in range(rprec):
                L = (left._unscaled_moment(i)).padded_list()
                R = (shift * right._unscaled_moment(i)).padded_list()
                c = cmp(L, R)
                if c: return c
        else:
            for i in range(rprec):
                L = (left._unscaled_moment(i)).padded_list()
                R = (right._unscaled_moment(i)).padded_list()
                c = cmp(L, R)
                if c: return c
        return 0

#        parent_prec = self.parent()._prec_cap
#        var_prec = min(self._var_prec, parent_prec[1], other._var_prec)
#        p_prec = min(len(self._moments), parent_prec[0], len(other._moments))
#        for i in range(p_prec):
#            selflist = _add_big_ohs_list(self._moments[i], [((p_prec - i) * self._cp).ceil(), var_prec])
#            otherlist = _add_big_ohs_list(other._moments[i], [((p_prec - i) * self._cp).ceil(), var_prec])
#            c = cmp(selflist, otherlist)
#            if c:
#                return c
#        return 0
    
#    def __nonzero__(self):
#        """
#        Checks that self is non-zero up to precision ...
#        """
#        parent_prec = self.parent()._prec_cap
#        var_prec = min(self._var_prec, parent_prec[1])
#        p_prec = min(len(self._moments), parent_prec[0])
#        for i in range(p_prec):
#            selflist = _add_big_ohs_list(self._moments[i], [((p_prec - i) * self._cp).ceil(), var_prec])
#            for c in selflist:
#                if c != 0:
#                    return True
#        return False
    
    def _coerce_map_from_(self, other):
        if isinstance(other, CoeffMod_OMS_Families_element) \
            and other._k  == self._k \
            and self._character == other._character \
            and self.base_ring().has_coerce_map_from(other.base_ring()):
            return True
        elif isinstance(other, CoeffMod_OMS_element):
            #kdiff = other._k - self._k
            #self._character = other._character
            #and self.base_ring().has_coerce_map_from(other.base_ring()):
            return False
        else:
            return False
    
    def is_zero(self, prec=None):
        if prec is None:
            prec = self.precision_relative()
        else:
            from sage.modular.pollack_stevens.coeffmod_OMS_families_space import _prec_cap_parser
            prec = _prec_cap_parser(prec)
        for i in range(min(prec[0], len(self._moments))):
            selflist = _add_big_ohs_list(self._moments[i], [((prec[0]-i)*self._cp).ceil(), prec[1]])
            for c in selflist:
                if c != 0:
                    return False
        return True
    
    def precision_relative(self):
        #RH: copied from coeffmod_OMS_element.py
        return [ZZ(len(self._moments)), self._var_prec]
    
    def precision_absolute(self):
        #RH: copied from coeffmod_OMS_element.py
        return [ZZ(len(self._moments) + self.ordp), self._var_prec]
    
    def valuation(self):
        #RH: adapted from coeffmod_OMS_element.py
        p = self.parent().prime()
        n = self.precision_relative()[0]
        return self.ordp + min([n] + [_padic_val_of_pow_series(self._unscaled_moment(a), p) for a in range(n) if not self._unscaled_moment(a).is_zero()])
    
    def normalize(self):
        #RH: adapted from coeffmod_OMS_element.py
        #Not tested
        V = self._moments.parent()
        R = V.base_ring()
        n = self.precision_relative()[0]
        p = self.parent().prime()
        for i in range(n):
            cutoff = ceil((n-i) * self._cp)
            f = self._moments[i].list()
            #f = [f[j].add_bigoh(cutoff) if f[j] != 0 else f[j] for j in range(min(len(f), self._var_prec))]
            f = [f[j].add_bigoh(cutoff) for j in range(min(len(f), self._var_prec))]
            self._moments[i] = R(f, self._var_prec)
        shift = self.valuation() - self.ordp
        if shift != 0:
            V = self.parent().approx_module(n-shift)
            self.ordp += shift
            p_to_shift = p**shift
            new_moments = []
            for i in range(n-shift):
                cutoff = ceil((n-shift-i) * self._cp)
                f = self._moments[i].list()
                f = [(f[j] // p_to_shift).add_bigoh(cutoff) if f[j] != 0 else f[j] for j in range(min(len(f), self._var_prec))]     #should shift prec in else?
                new_moments.append(R(f, self._var_prec))
            self._moments = V(new_moments)
        return self
    
    def moment(self, n):
        #RH: "copied" from dist.pyx
        r"""
        Returns the `n`-th moment.
        """
        if self.ordp == 0:
            return self._unscaled_moment(n)
        else:
            return self.parent().prime()**(self.ordp) * self._unscaled_moment(n)
    
    def _unscaled_moment(self, n):
        #RH: "copied" from coeffmod_OMS_element.py
        r"""
        Returns the `n`-th moment, unscaled by the overall power of p stored in self.ordp.
        """
        return self._moments[n]
    
    def reduce_precision(self, new_prec):
        rprec = self.precision_relative()
        try:
            p_prec, var_prec = new_prec
        except TypeError:   #if new_prec isn't iterable (wrong size is ValueError)
            p_prec = new_prec
            var_prec = rprec[1]
        if p_prec > rprec[0] or var_prec > rprec[1]:
            raise ValuError("Precisions specified must be less than current precisions.")
        moments = self._moments[:p_prec]
        ordp = self.ordp
        return CoeffMod_OMS_Families_element(moments, self.parent(), ordp=ordp, check=False, var_prec=var_prec)
    
    def solve_diff_eqn(self):
        #Do something about ordp
        if self.is_zero():
            return self.parent()(0)
        if self._moments[0] != 0:
            raise ValueError("Family of distribution must have total measure 0 to be in image of difference operator.")
        M = len(self._moments)
        if M == 1:
            return self.parent().zero()
        from sage.modular.pollack_stevens.coeffmod_OMS_space import OverconvergentDistributions
        R = self.base_ring().base_ring()
        DD = self.parent()
        D = OverconvergentDistributions(0, base=R, prec_cap=M, character=DD._character, adjuster=DD._adjuster, act_on_left=DD.action().is_left(), dettwist=DD._dettwist)
        V = D.approx_module(M)
        Elem = D.Element
        v = V([R.zero(), R.one()] + [R.zero()]*(M-2))
        mu = Elem(v, D, ordp=0, check=False)
        #Normalize mu?
        mus = self.moment(1) * mu.solve_diff_eqn().lift(DD)
        for j in range(2, M):
            mu._moments[j] = R.one()
            mu._moments[j-1] = R.zero()
            mus += self.moment(j) * mu.solve_diff_eqn().lift(DD)
        return mus

def _padic_val_of_pow_series(f, p=None):
    r"""
        Given a power series ``f`` return its ``p``-adic valuation, i.e. the
        minimum ``p``-adic valuation of its coefficients
    """
    if f == 0:
        return Infinity
    return min([coeff.valuation(p) for coeff in f if not coeff.is_zero()])

def _padic_val_unit_of_pow_series(f, p=None):
    r"""
        Given a power series ``f``, return a tuple `(v, u)`, where `v` is the
        ``p``-adic valuation of ``f`` and `u` is the ``p``-adic unit part of
        ``f``.
    """
    if f == 0:
        return (Infinity, self.parent()(0,0))
    if p is None:
        p = self.parent().base_ring().prime()
    v = _padic_val_of_pow_series(f, p)
    u = f.parent()([(coeff / (p ** v)) for coeff in f])
    return (v, u)

def _add_big_ohs_list(f, prec_cap):
    r"""
    Returns a (padded) list of length (at most) ``prec_cap``[1] of the coefficients
    up to `p`-adic precision ``prec_cap``[0]. The input is checked. 
    
    INPUT:
    
        - ``f`` -- a power series over a `p`-adic ring
        - ``prec_cap`` -- a pair [``p_prec``, ``var_prec``, where ``p_prec`` is
        the desired `p`-adic precision of the coefficients of f and ``var_prec``
        is the desired variable-adic precision
    
    OUTPUT:
    
        - a list of coefficients of ``f`` of precision up to (at most)
        `O(p^{\text{p_prec}})`. The length of this list is the minimum of the
        precision of ``f`` and of ``var_prec``. The padded part should be exact zeroes.
    
    EXAMPLES::
    
        sage: from sage.modular.pollack_stevens.coeffmod_OMS_families_element import _add_big_ohs_list
        sage: R = PowerSeriesRing(Zp(3), 'w')
        sage: f = R([0,-1,2,0,-3])
        sage: _add_big_ohs_list(f, [10,7])
        [O(3^10), 2 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + O(3^10), 2 + O(3^10), O(3^10), 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^11 + 2*3^12 + 2*3^13 + 2*3^14 + 2*3^15 + 2*3^16 + 2*3^17 + 2*3^18 + 2*3^19 + 2*3^20 + O(3^21), 0, 0]
        sage: _add_big_ohs_list(f, [2,4])
        [O(3^2), 2 + 2*3 + O(3^2), 2 + O(3^2), O(3^2)]
    """
    p_prec, var_prec = prec_cap
    flist = f.padded_list(var_prec)
    for i in range(0, min(f.degree(), var_prec)):
        flist[i] = flist[i].add_bigoh(p_prec)
    return flist