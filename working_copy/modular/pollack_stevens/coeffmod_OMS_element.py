# This should be cythoned once it's done.

from sage.rings.arith import binomial, bernoulli
from sage.modular.pollack_stevens.coeffmod_element import CoefficientModuleElement_generic
from sage.rings.integer_ring import ZZ

class CoeffMod_OMS_element(CoefficientModuleElement_generic):
    # Implementation currently ignores ordp
    def __init__(self, moments, parent, ordp=0, check=True):
        CoefficientModuleElement_generic.__init__(self, parent)
        #TODO: finish this
        if check:
            if isinstance(moments, CoeffMod_OMS_element):
                moments = moments._moments.change_ring(parent.base_ring())
            elif hasattr(moments, '__len__'):
                M = len(moments)
                moments = parent.approx_module(M)(moments)
            elif moments == 0:
                moments = parent.approx_module(parent._prec_cap)(moments)
            else:
                moments = parent.approx_module(parent._prec_cap)([moments]*parent._prec_cap)
        self._moments = moments
        #if var_prec is None:
        #    self._var_prec = parent._prec_cap[1]
        #else:
        #    self._var_prec = var_prec
        self.ordp = 0   #eventually change this maybe
    
    #def _relprec(self):
    #    return len(self._moments)
    def _repr_(self):
        return repr(self._moments)
    
    def _add_(self, right):
        return self._addsub(right, False)
    
    def _sub_(self, right):
        return self._addsub(right, True)
    
    def _addsub(self, right, negate):
        self_moments = self._moments
        right_moments = right._moments
        new_prec = min(len(self._moments), len(right._moments))
        if negate:
            new_moments = [self_moments[i] - right_moments[i] for i in range(new_prec)] 
        else:
            new_moments = [self_moments[i] + right_moments[i] for i in range(new_prec)]
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
        #print other, type(other)
        prec = min(len(self._moments), parent_prec, len(other._moments))    #this is a problem for checking bla == 0
        for i in range(prec):
            if self._moments[i].add_bigoh(prec - i) != other._moments[i].add_bigoh(prec - i):
                return False
        return True
    
    def __nonzero__(self):
        """
        Checks that self is non-zero up to precision ...
        """
        parent_prec = self.parent()._prec_cap
        prec = min(len(self._moments), parent_prec)
        for i in range(prec):
            if self._moments[i].add_bigoh(prec - i) != 0:
                return True
        return False
    
    def is_zero(self, prec=None):
        if prec is None:
            prec = self.precision_relative()
        for i in range(min(prec, len(self._moments))):
            if self._moments[i].add_bigoh(prec - i) != 0:
                return False
        return True
    
    def precision_relative(self):
        return ZZ(len(self._moments))
    
    def precision_absolute(self):
        return ZZ(len(self._moments)) + self.ordp
    
    def normalize(self):
        #customized
        M = len(moments)
        for i in range(M):
            self._moments[i] = self._moments[i].add_bigoh(M-i)
        return self
    
    def moment(self, n):
        r"""
        Returns the `n`-th moment.

        INPUT:

        - ``n`` -- an integer or slice, to be passed on to moments.

        OUTPUT:

        - the `n`-th moment, or a list of moments in the case that `n`
          is a slice.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        if self.ordp == 0:
            return self._unscaled_moment(n)
        else:
            return self.parent().prime()**(self.ordp) * self._unscaled_moment(n)
    
    def _unscaled_moment(self, n):
        r"""
        Returns the `n`-th moment, unscaled by the overall power of p stored in self.ordp.
        """
        return self._moments[n]
    
    def reduce_precision(self, new_prec):
        pass    #TODO
    
    def lift(self, var_prec=None, variable_name='w'):
        V = self.parent().lift(var_prec, variable_name)
        #k = V._k
        #p = V._p
        #prec = V.precision_cap()
        #R = V.base_ring()
        return V([self.moment(i) for i in range(len(self._moments))])
    
    def solve_diff_eqn(self):
        r"""
        Solves the difference equation.

        See Theorem 4.5 and Lemma 4.4 of [PS].

        OUTPUT:

        - a distribution v so that self = v | Delta, where Delta = [1, 1; 0, 1] - 1.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        # assert self._moments[0][0]==0, "not total measure zero"
        # print "result accurate modulo p^",self.moment(0).valuation(self.p)
        #v=[0 for j in range(0,i)]+[binomial(j,i)*bernoulli(j-i) for j in range(i,M)]
        M = self.precision_relative()
        R = self.parent().base_ring()
        K = R.fraction_field()
        V = self._moments.parent()
        v = [K(0) for i in range(M)]
        bern = [bernoulli(i) for i in range(0,M,2)]
        minhalf = ~K(-2)
        for m in range(1,M):
            scalar = K(self.moment(m) / m)
            # bernoulli(1) = -1/2; the only nonzero odd bernoulli number
            v[m] += m * minhalf * scalar
            for j in range(m-1,M,2):
                v[j] += binomial(j,m-1) * bern[(j-m+1)//2] * scalar
        p = self.parent().prime()
        if p == 0:
            if R.is_field():
                ans = self.parent()(v)
                ans.ordp = 0    #Is this redundant?
            else:
                newparent = self.parent().change_ring(K)
                ans = newparent(v)
        else:
            ans = self.parent().zero()
            ans.ordp = min(a.valuation(p) for a in v)
            if ans.ordp < 0:
                scalar = K(p) ** (-ans.ordp)
                ans._moments = V([R(a * scalar) for a in v])
            elif ans.ordp > 0:
                scalar = K(p) ** ans.ordp
                ans._moments = V([R(a // scalar) for a in v])
            else:
                ans._moments = V([R(a) for a in v])
        return ans