# This should be cythoned once it's done.

from sage.modular.pollack_stevens.coeffmod_element import CoefficientModuleElement_generic
from sage.modular.pollack_stevens.coeffmod_OMS_element import CoeffMod_OMS_element
from sage.rings.integer_ring import ZZ

class CoeffMod_OMS_Families_element(CoefficientModuleElement_generic):
    # Implementation currently ignores ordp
    def __init__(self, moments, parent, ordp=0, check=True, var_prec=None):
        CoefficientModuleElement_generic.__init__(self, parent)
        #TODO: finish this
        if var_prec == None:
            var_prec = parent._prec_cap[1]
        self._var_prec = var_prec
        if check:
            if isinstance(moments, CoeffMod_OMS_Families_element):
                moments = moments._moments.change_ring(parent.base_ring())
            #if isinstance(moments, CoeffMod_OMS_element):
            #    moments = self.parent().approx_module(p_prec=len(moments._moments), var_prec=self.parent().precision_cap()[1])(moments._moments)
            elif hasattr(moments, '__len__'):
                M = len(moments)
                moments = parent.approx_module(M)(moments)
            elif moments == 0:
                moments = parent.approx_module(parent._prec_cap[0], parent._prec_cap[1])(moments)
            else:
                moments = parent.approx_module(parent._prec_cap[0], parent.precision_cap()[1])([moments]*parent._prec_cap[0])
        self._moments = moments
        #if var_prec is None:
        #    self._var_prec = parent._prec_cap[1]
        #else:
        #    self._var_prec = var_prec
        self.ordp = 0   #eventually change this maybe
        p = parent._p
        self._cp = (p-2) / (p-1)
    
    #def _relprec(self):
    #    return len(self._moments)
    def _repr_(self):
        return repr(self._moments)
    
    def _add_(self, right):
        return self._addsub(right, False)
    
    def _sub_(self, right):
        return self._addsub(right, True)
    
    def _addsub(self, right, negate):
        new_moments = []
        new_prec = [min(len(self._moments), len(right._moments)), min(self._var_prec, right._var_prec)]
        R = self.parent().base_ring()
        for i in range(new_prec[0]):
            selflist = _add_big_ohs_list(self._moments[i], new_prec)
            rightlist = _add_big_ohs_list(right._moments[i], new_prec)
            for j in range(len(selflist)):
                if negate:
                    selflist[j] -= rightlist[j]
                else:
                    selflist[j] += rightlist[j]
            new_moments.append(R(selflist))
        return self.parent()(new_moments)
    
    def _lmul_(self, right):
        """
        Scalar multiplication self*right.
        """
        return self.parent()(self._moments*right)
    
    def _rmul_(self, left):
        """
        Scalar multiplication left*self.
        """
        return self._lmul_(left)
    
    def __eq__(self, other):
        parent_prec = self.parent()._prec_cap
        var_prec = min(self._var_prec, parent_prec[1], other._var_prec)
        p_prec = min(len(self._moments), parent_prec[0], len(other._moments))
        for i in range(p_prec):
            selflist = _add_big_ohs_list(self._moments[i], [((p_prec - i) * self._cp).ceil(), var_prec])
            otherlist = _add_big_ohs_list(other._moments[i], [((p_prec - i) * self._cp).ceil(), var_prec])
            for j in range(len(selflist)):
                if selflist[j] != otherlist[j]:
                    return False
        return True
    
    def __nonzero__(self):
        """
        Checks that self is non-zero up to precision ...
        """
        parent_prec = self.parent()._prec_cap
        var_prec = min(self._var_prec, parent_prec[1])
        p_prec = min(len(self._moments), parent_prec[0])
        for i in range(p_prec):
            selflist = _add_big_ohs_list(self._moments[i], [((p_prec - i) * self._cp).ceil(), var_prec])
            for c in selflist:
                if c != 0:
                    return True
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
        return [ZZ(len(self._moments)), self._var_prec]
    
    def precision_absolute(self):
        return [ZZ(len(self._moments)) + self.ordp, self._var_prec]
    
    def normalize(self):
        #customized
        p = self.parent()._p
        moments = self._moments
        M = len(moments)
        R = self.parent().base_ring()
        #new_moments = [m.add_bigoh(self._var_prec) for m in self._moments]
        new_moments = []
        for i in range(M):
            cutoff = ceil((M-i) * self._cp)
            f = moments[i].list()
            f = [f[j].add_bigoh(cutoff) for j in range(self._var_prec)]
            moments[i] = R(f, self._var_prec)
        return self
    
    def reduce_precision(self, new_prec):
        pass    #TODO
    
    def solve_diff_eqn(self):
        M = len(self._moments)
        D0 = self.parent().specialize(0)
        mus = self.parent().zero()
        v = [1] + [0]*(M-1)
        for j in range(M):
            mu = D0(v)
            nus = self.parent()(mu.solve_diff_eqn())
            mus += nus
            v.pop()
            v.insert(0, 0)
        return

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