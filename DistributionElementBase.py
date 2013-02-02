cdef class DistributionElementBase(ModuleElement):
    r"""
        The main p-adic distribution class, implemented as per the paper
        'Overconvergent Modular Symbols and p-adic L-functions' by Pollack
        & Stevens
    """

    def scale(self,left):
        r"""
        Scales the moments of the distribution by `left`

        INPUT:

        - ``left`` -- scalar

        OUTPUT:

        - Scales the moments by `left`

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
            sage: D = Distributions(5, 7, 15)
            sage: v = D([1,2,3,4,5]); v
            (1 + O(7^5), 2 + O(7^4), 3 + O(7^3), 4 + O(7^2), 5 + O(7))
            sage: v.scale(2)
            (2 + O(7^5), 4 + O(7^4), 6 + O(7^3), 1 + 7 + O(7^2), 3 + O(7))
        """
        if isinstance(self, Dist_long) and isinstance(left, (Integer, pAdicCappedRelativeElement, pAdicCappedAbsoluteElement, pAdicFixedModElement)):
            return self._lmul_(left)
        R = left.parent()
        base = self.parent().base_ring()
        if base is R:
            return self._lmul_(left)
        elif base.has_coerce_map_from(R):
            return self._lmul_(base(left))
        else:
            from sage.categories.pushout import pushout
            new_base = pushout(base, R)
            V = self.parent().change_ring(new_base)
            scalar = new_base(left)
            return V([scalar * new_base(self.moment(i)) for i in range(self.precision_absolute())])

    def is_zero(self, p=None, M=None):
        r"""
        Returns True if all of the moments are either identically zero or zero modulo p^M

        INPUT:

        - ``p`` -- prime

        - ``M`` -- precision

        OUTPUT:

        - True/False

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
            sage: D = Distributions(5, 7, 15)
            sage: v = D([1,2,3,4,5]); v
            (1 + O(7^5), 2 + O(7^4), 3 + O(7^3), 4 + O(7^2), 5 + O(7))
            sage: v.is_zero()
            False
            sage: v = D(5*[0])
            sage: v.is_zero()
            True
        """
        n = self.precision_absolute()
        if M is None:
            return all([self.moment(a).is_zero() for a in range(n)])
        else:
            return all([self.moment(a).valuation(p) >= M for a in range(n)])

    def find_scalar(self, other, p, M = None, check=True):
        r"""
        Returns an ``alpha`` with ``other = self * alpha``, or raises a ValueError.

        It will also raise a ValueError if this distribution is zero.

        INPUT:

        - ``other`` -- another distribution

        - ``p`` -- an integral prime (only used if the parent is not a Symk)

        - ``M`` -- (default: None) an integer, the relative precision
          to which the scalar must be determined

        - ``check`` -- (default: True) boolean, whether to validate
          that ``other`` is actually a multiple of this element.

        OUTPUT:

        - A scalar ``alpha`` with ``other = self * alpha``.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
            sage: D = Distributions(5, 7, 15)
            sage: v = D([1,2,3,4,5])
            sage: w = D([3,6,9,12,15])
            sage: v.find_scalar(w,p=7)
            3 + O(7^5)

            sage: u = D([1,4,9,16,25])
            sage: v.find_scalar(u,p=7)
            Traceback (most recent call last):
            ...
            ValueError: not a scalar multiple

        """
        i = 0
        n = self.precision_absolute()
        if n != other.precision_absolute():
            raise ValueError("other should have the same number of moments")
        verbose("n = %s"%n)
        verbose("moment 0")
        a = self.moment(i)
        verbose("a = %s"%(a))
        padic = isinstance(a.parent(), pAdicGeneric)
        if self.parent().is_symk():
            while a == 0:
                if other.moment(i) != 0:
                    raise ValueError("not a scalar multiple")
                i += 1
                verbose("moment %s"%i)
                try:
                    a = self.moment(i)
                except IndexError:
                    raise ValueError("self is zero")
            alpha = other.moment(i) / a
            if check:
                i += 1
                while i < n:
                    verbose("comparing moment %s"%i)
                    if self.moment(i) != alpha * other.moment(i):
                        raise ValueError("not a scalar multiple")
                    i += 1
        else:
            p = self.parent().prime()
            v = a.valuation(p)
            while v >= n - i:
                i += 1
                verbose("p moment %s"%i)
                try:
                    a = self.moment(i)
                except IndexError:
                    raise ValueError("self is zero")
                v = a.valuation(p)
            relprec = n - i - v
            verbose("p=%s, n-i=%s\nself.moment=%s, other.moment=%s"%(p, n-i, a, other.moment(i)),level=2)
            if padic:
                alpha = (other.moment(i) / a).add_bigoh(n-i)
            else:
                alpha = (other.moment(i) / a) % p**(n-i)
            verbose("alpha = %s"%(alpha))
            while i < n-1:
                i += 1
                verbose("comparing p moment %s"%i)
                a = self.moment(i)
                if check:
                    verbose("self.moment=%s, other.moment=%s"%(a, other.moment(i)))
                    if (padic and other.moment(i) != alpha * a) or \
                       (not padic and other.moment() % p**(n-i) != alpha * a % p**(n-i)):
                        raise ValueError("not a scalar multiple")
                v = a.valuation(p)
                if n - i - v > relprec:
                    verbose("Reseting alpha: relprec=%s, n-i=%s, v=%s"%(relprec, n-i, v))
                    relprec = n - i - v
                    if padic:
                        alpha = (other.moment(i) / a).add_bigoh(n-i)
                    else:
                        alpha = (other.moment(i) / a) % p**(n-i)
                    verbose("alpha=%s"%(alpha))
            if relprec < M:
                raise ValueError("result not determined to high enough precision")
        try:
            return self.parent().base_ring()(alpha)
        except ValueError:
            return alpha

    cpdef ModuleElement _rmul_(self, RingElement _left):
        """
        Scalar multiplication.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
        """
        return self._lmul_(_left)

    def diagonal_valuation(self, p=None):
        """
        Returns the largest `m` so that this distribution lies in `Fil^m`.

        INPUT:

        - ``p`` -- (default: None) a positive integral prime

        OUTPUT:

        - the largest integer `m` so that `p^m` divides the `0`-th
          moment, `p^{m-1}` divides the first moment, etc.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
            sage: D = Distributions(8, 7, 15)
            sage: v = D([7^(5-i) for i in range(1,5)])
            sage: v
            (O(7^4), O(7^3), O(7^2), O(7))
            sage: v.diagonal_valuation(7)
            4
        """
        if p is None:
            p = self.parent()._p
        n = self.precision_absolute()
        return min([n] + [a + self.moment(a).valuation(p) for a in range(n)])

    def valuation(self, p=None):
        """
        Returns the minimum valuation of any moment.

        INPUT:

        - ``p`` -- (default: None) a positive integral prime

        OUTPUT:

        - 

        .. WARNING::

            Since only finitely many moments are computed, this
            valuation may be larger than the actual valuation of this
            distribution.  Moreover, since distributions are
            normalized so that the top moment has precision 1, this valuation may be smaller than the actual valuation (for example, if the actual valuation is 2)

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
            sage: D = Distributions(8, 7, 15)
            sage: v = D([7^(5-i) for i in range(1,5)])
            sage: v
            (O(7^4), O(7^3), O(7^2), O(7))
            sage: v.valuation(7)
            1
        """
        r"""
        Returns the highest power of `p` which divides all moments of the distribution
        """
        if p is None:
            p = self.parent()._p
        n = self.precision_absolute()
        return min([self.moment(a).valuation(p) for a in range(n)])

    def specialize(self, new_base_ring=None):
        """
        Returns the image of this overconvergent distribution under
        the canonical projection from distributions of weight k to
        Sym^k.

        INPUT:

        - ``new_base_ring`` -- (default: None) a ring giving the
          desired base ring of the result.

        OUTPUT:

        - An element of Sym^k(K), where K is the specified base ring.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
        """
        k=self.parent()._k
        if k < 0:
            raise ValueError("negative weight")
        if self.precision_absolute() < k+1:
            raise ValueError("not enough moments")
        V = self.parent().specialize(new_base_ring)
        new_base_ring = V.base_ring()
        return V([binomial(k, j) * (-1)**j * new_base_ring(self.moment(j)) for j in range(k+1)])

    def lift(self, p=None, M=None, new_base_ring=None):
        r"""
        Lifts a distribution or element of Sym^k to an overconvergent distribution.

        INPUT:

        - ``p`` -- (default: None) a positive integral prime.  If None
          then p must be available in the parent.

        - ``M`` -- (default: None) a positive integer giving the
          desired number of moments.

        - ``new_base_ring`` -- (default: None) a ring giving the
          desired base ring of the result.

        OUTPUT:

        - An overconvergent distribution with `M` moments whose image
          under the specialization map is this element.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
        """
        V = self.parent().lift(p, M, new_base_ring)
        k = V._k
        p = V.prime()
        M = V.precision_cap()
        R = V.base_ring()
        moments = [R(self.moment(j) * (-1)**j / binomial(k, j)) for j in range(k+1)]
        zero = R(0)
        moments.extend([zero] * (M - k - 1))
        mu = V(moments)
        #val = mu.valuation()
        #if val < 0:
        #    # This seems unnatural
        #    print "scaling by %s^%s to keep things integral"%(p, -val)
        #    mu *= p**(-val)
        return mu

    def _is_malformed(self):
        n = self.precision_absolute()
        for i in range(n):
            if self.moment(i).precision_absolute() < n - i:
                return True
        return False

    def act_right(self,gamma):
        r"""
        The image of this element under the right action by a
        `2 \times 2` matrix.

        INPUT:

        - ``gamma`` -- the matrix by which to act

        OUTPUT:

        - ``self | gamma``

        .. NOTE::

            You may also just use multiplication ``self * gamma``.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
        """
        return self.parent()._act(self, gamma)
