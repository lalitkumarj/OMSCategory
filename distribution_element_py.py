cdef class DistributionElementPy(DistributionElementBase):
    r"""
    A distribution is stored as a vector whose `j`-th entry is the `j`-th moment of the distribution.

    The `j`-th entry is stored modulo `p^(N-j)` where `N` is the total number of moments.
    (This is the accuracy that is maintained after acting by `\Gamma_0(p)`.)

    INPUTS:

    - ``moments`` -- the list of moments.  If ``check == False`` it
      must be a vector in the appropriate approximation module.

    - ``parent`` -- a :class:`distributions.Distributions_class` or
      :class:`distributions.Symk_class` instance

    - ``check`` -- (default: True) boolean, whether to validate input

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.distributions import Distributions
    """
    def __init__(self,moments, parent, check=True):
        """
        Initialization.

        TESTS::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
        """
        Dist.__init__(self,parent)
        if check:
            base = parent.base_ring()
            try:
                M = len(moments)
            except TypeError:
                M = 1
                moments = [moments]
            moments = parent.approx_module(M)(moments)
        self.moments = moments

    def __reduce__(self):
        return (self.__class__,(self.moments,self.parent(),False))

    cdef Dist_vector _new_c(self):
        r"""
        Creates an empty distribution.

        OUTPUT:

        - A distribution with no moments.  The moments are then filled
          in by the calling function.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
        """
        cdef Dist_vector ans = PY_NEW(Dist_vector)
        ans._parent = self._parent
        return ans

    def _repr_(self):

        r"""
        String representation.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
        """
        r"""
        Displays the moments of the distribution
        """
        self.normalize()
        if len(self.moments) == 1:
            return repr(self.moments[0])
        else:
            return repr(self.moments)

    def _rational_(self):
        """
        Convert to a rational number.

        EXAMPLES::

            sage: D = Symk(0); d = D(4/3); d
            4/3
            sage: QQ(d)
            4/3

        We get a TypeError if there is more than 1 moment::

            sage: D = Symk(1); d = D([1,2]); d
            (1, 2)
            sage: QQ(d)
            Traceback (most recent call last):
            ...
            TypeError: k must be 0
        """
        if len(self.moments) == 1:
            return QQ(self.moments[0])
        raise TypeError, "k must be 0"

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
        r"""
        Returns the `n`-th moment
        """
        return self.moments[n]

    cpdef ModuleElement _add_(self, ModuleElement _right):
        r"""
        Sum of two distributions.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        cdef Dist_vector ans = self._new_c()
        cdef Dist_vector right = _right
        if len(self.moments) == len(right.moments):
            ans.moments = self.moments + right.moments
        elif len(self.moments) < len(right.moments):
            ans.moments = self.moments + right.moments[:len(self.moments)]
        else:
            ans.moments = self.moments[:len(right.moments)] + right.moments
        return ans

    cpdef ModuleElement _sub_(self, ModuleElement _right):
        r"""
        Difference of two distributions.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        cdef Dist_vector ans = self._new_c()
        cdef Dist_vector right = _right
        if len(self.moments) == len(right.moments):
            ans.moments = self.moments - right.moments
        elif len(self.moments) < len(right.moments):
            ans.moments = self.moments - right.moments[:len(self.moments)]
        else:
            ans.moments = self.moments[:len(right.moments)] - right.moments
        return ans

    cpdef ModuleElement _lmul_(self, RingElement right):
        r"""
        Scalar product of a distribution with a ring element that coerces into the base ring.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        cdef Dist_vector ans = self._new_c()
        ans.moments = self.moments * right
        return ans

    def precision_absolute(self):
        r"""
        Returns the precision of this distribution.

        The precision is just the number of moments stored, which is
        also k+1 in the case of Sym^k(R).  For overconvergent
        distributions, the precision is the integer `m` so that the
        sequence of moments is known modulo `Fil^m`.

        OUTPUT:

        - An integer giving the number of moments.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        r"""
        Returns the number of moments of the distribution
        """
        return len(self.moments)

    cdef int _cmp_c_impl(left, Element right) except -2:
        r"""
        Comparison.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        return cmp(left.moments, right.moments)

    def zero(self):
        r"""
        

        OUTPUT:

        - 

        """
        cdef Dist_vector ans = self._new_c()
        ans.moments = 0 * self.moments # could make this faster
        return ans

    cpdef normalize(self):
        r"""
        Normalize by reducing modulo `Fil^N`, where `N` is the number of moments.

        If the parent is Symk, then normalize has no effect.  If the
        parent is a space of distributions, then normalize reduces the
        `i`-th moment modulo `p^{N-i}`.

        OUTPUT:

        - this distribtion, after normalizing.

        .. WARNING::

            This function modifies the distribution in place as well as returning it.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        p = self.parent()._p
        if not self.parent().is_symk(): # non-classical
            if self.valuation() < 0:
                verbose("Negative valuation!")
                verbose("%s"%(self.moments))
            #assert self.valuation() >= 0, "moments not integral in normalization"
            V = self.moments.parent()
            R = V.base_ring()
            n = self.precision_absolute()
            if isinstance(R, pAdicGeneric):
                self.moments = V([self.moment(i).add_bigoh(n-i) for i in range(n)])
            else:
                self.moments = V([self.moment(i)%(p**(n-i)) for i in range(n)])
        return self

    def reduce_precision(self, M):
        r"""
        Only hold on to `M` moments.

        INPUT:

        - ``M`` -- a positive integer less than the precision of this
          distribution.

        OUTPUT:

        - a new distribution with `M` moments equal to the first `M`
          moments of this distribution.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        assert M<=self.precision_absolute(),"not enough moments"

        cdef Dist_vector ans = self._new_c()
        ans.moments = self.moments[:M]
        return ans

    def solve_diff_eqn(self):
        r"""
        Solves the difference equation.

        See Theorem 4.5 and Lemma 4.4 of [PS].

        OUTPUT:

        - a distribution v so that self = v | Delta, where Delta = [1, 1; 0, 1] - 1.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        # assert self.moments[0][0]==0, "not total measure zero"
        # print "result accurate modulo p^",self.moment(0).valuation(self.p)
        #v=[0 for j in range(0,i)]+[binomial(j,i)*bernoulli(j-i) for j in range(i,M)]
        M = self.precision_absolute()
        K = self.parent().base_ring().fraction_field()
        V = self.moments.parent()
        v = [K(0) for i in range(M)]
        for m in range(1,M):
            scalar = K(self.moment(m) / m) # division can take us out of the base ring.
            for j in range(m-1,M):
                v[j] += binomial(j,m-1) * bernoulli(j-m+1) * scalar
        cdef Dist_vector ans = self._new_c()
        ans.moments = V(v)
        return ans

    #def lift(self):
    #    r"""
    #    Increases the number of moments by `1`.
    #    """
    #    n = len(self.moments)
    #    if n >= self.parent()._prec_cap:
    #        raise ValueError("Cannot lift above precision cap")
    #    cdef Dist_vector ans = self._new_c()
    #    ans.moments = self.parent().approx_module(n+1)(list(self.moments) + [0])
    #    return ans
