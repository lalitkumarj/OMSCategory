from sage.categories.all import Category

class MSCoefficientModule(Category):
    def __repr__(self):
        print "Category of Corfficient Modules"

    class ParentMethods:
        def _repr_(self):
            pass

        def prime():
            """
       Return prime `p` such that this is a space of `p`-adic distributions.

       In case this space is Symk of a non-padic field, this makes no
       sense, and we raise a ValueError.

       OUTPUT:

       - bool

       EXAMPLES::

           sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
           sage: D = Distributions(0, 7); D
           Space of 7-adic distributions with k=0 action and precision cap 20
           sage: D.prime()
           7
           sage: D = Symk(4, base=GF(7)); D
           Sym^4 (Finite Field of size 7)^2
           sage: D.prime()
           Traceback (most recent call last):
           ...
           ValueError: not a space of p-adic distributions

       But Symk of a `p`-adic field does work::

           sage: D = Symk(4, base=Qp(7)); D
           Sym^4 Q_7^2
           sage: D.prime()
           7
           sage: D.is_symk()
           True
       """
            if self._p is None:
                raise ValueError, "not a space of p-adic distributions"

        def random_element():
            pass

        def weight(self):
            """
            Return the weight of this distribution space.  The standard
            caveat applies, namely that the weight of `Sym^k` is
            defined to be `k`, not `k+2`.

            OUTPUT:

            - nonnegative integer

            EXAMPLES::

                sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
                sage: D = Distributions(0, 7); D
                Space of 7-adic distributions with k=0 action and precision cap 20
                sage: D.weight()
                0
                sage: Distributions(389, 7).weight()
                389
            """
            return self._k

        def zero_element(self, M=None):
            """
            Return zero element in the M-th approximating module.

        INPUT:

        - `M` -- None (default), or a nonnegative integer, less than or equal to the precision cap

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Distributions(0, 7, 4); D
            Space of 7-adic distributions with k=0 action and precision cap 4
            sage: D.zero_element()
            (0, 0, 0, 0)
            sage: D.zero_element(0)
            ()
            sage: D.zero_element(1)
            0
            sage: D.zero_element(2)
            (0, 0)
            sage: D.zero_element(3)
            (0, 0, 0)
            sage: D.zero_element(4)
            (0, 0, 0, 0)
            sage: D.zero_element(5)
            Traceback (most recent call last):
            ...
            ValueError: M must be less than or equal to the precision cap
            """
            return self(self.approx_module(M)(0))

        def precision_cap(self):
            """
            Return the precision cap on distributions.

            EXAMPLES::

                sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
                sage: D = Distributions(0, 7, 10); D
                Space of 7-adic distributions with k=0 action and precision cap 10
                sage: D.precision_cap()
                10
                sage: D = Symk(389, base=QQ); D
                Sym^389 Q^2
                sage: D.precision_cap()
                390
            """
            return self._prec_cap

        def approx_module(self, M=None):
            """
        Return the M-th approximation module, or if M is not specified,
        return the largest approximation module.

        INPUT::

        - `M` -- None or nonnegative integer that is at most the precision cap

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
            sage: D = Distributions(0, 5, 10)
            sage: D.approx_module()
            Ambient free module of rank 10 over the principal ideal domain 5-adic Ring with capped absolute precision 10
            sage: D.approx_module(1)
            Ambient free module of rank 1 over the principal ideal domain 5-adic Ring with capped absolute precision 10
            sage: D.approx_module(0)
            Ambient free module of rank 0 over the principal ideal domain 5-adic Ring with capped absolute precision 10

        Note that M must be at most the precision cap, and must be nonnegative::

            sage: D.approx_module(11)
            Traceback (most recent call last):
            ...
            ValueError: M must be less than or equal to the precision cap
            sage: D.approx_module(-1)
            Traceback (most recent call last):
            ...
            ValueError: rank (=-1) must be nonnegative
            """
            if M is None:
                M = self._prec_cap
            elif M > self._prec_cap:
                raise ValueError("M must be less than or equal to the precision cap")
            return self.base_ring()**M

        def clear_cache(self):
            """
           Clear some caches that are created only for speed purposes.

           EXAMPLES::

               sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
               sage: D = Distributions(0, 7, 10)
               sage: D.clear_cache()
         """
            self.approx_module.clear_cache()
            self._act.clear_cache()


        def basis(self, M=None):
            """
    Return a basis for this space of distributions.

    EXAMPLES::

    sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
    sage: D = Distributions(0, 7, 4); D
    Space of 7-adic distributions with k=0 action and precision cap 4
    sage: D.basis()
    [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]

    sage: D = Symk(3, base=QQ); D
    Sym^3 Q^2
    sage: D.basis()
    [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
            """
            V = self.approx_module(M)
            return [self(v) for v in V.basis()]


class ElementMethods:
        def normalize():
            pass
        def valuation(Self, p=None);
        def scale
        def is_zero
        def _rmul_
        def lift():
            pass
        def find_scalar()
        def _is_malformed()
        def act_right()

