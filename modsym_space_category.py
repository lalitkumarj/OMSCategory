from sage.categories.all import Category

class ModsymSpace(Category):
    def super_categories(self):
        return []

    def __repr__(self):
        return "Category of Modular Symbol Spaces"


    class ParentMethods:
        def coefficient_module():
            pass
        
        def prime(self):
            if self.p is None:
                raise ValueError, "not a space of p-adic distributions"


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
            return (self.M_max,self.d_max)

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

        def manin_relations():
            pass

        def basis(self, M=None):
            pass
        def char(self):
            return self.coefficient_module.char()
        def weight(self):
            pass

    class ElementMethods:
        def __init__(self,map_data,parent,construct = False):
            #ModuleElement.__init__(self, parent)
            if construct:
                self._map = map_data
            else:
                self._map = ManinMap(parent._coefficients, parent._source, map_data)
            self.parent = parent
        def _repr_(self):
            r"""
            Returns the print representation of the symbol.

            EXAMPLES::

                sage: E = EllipticCurve('11a')
                sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
                sage: phi = ps_modsym_from_elliptic_curve(E)
                sage: phi._repr_()
                'Modular symbol with values in Sym^0 Q^2'
            """
            return "Modular symbol with values in %s"%(self.parent().coefficient_module())
            
        def dict(self):
            r"""
            Returns dictionary on the modular symbol self, where keys are generators and values are the corresponding values of self on generators

            EXAMPLES::

                sage: E = EllipticCurve('11a')
                sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
                sage: phi = ps_modsym_from_elliptic_curve(E)
                sage: phi.dict()
                {[1 0]
                [0 1]: -1/5, [ 0 -1]
                [ 1  3]: 3/2, [-1 -1]
                [ 3  2]: -1/2}
            """
            return self._map._dict
        def values(self):
            r"""
            Returns the values of the symbol self on our chosen generators (generators are listed in self.dict().keys())

            EXAMPLES::

                 sage: E = EllipticCurve('11a')
                 sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
                 sage: phi = ps_modsym_from_elliptic_curve(E)
                 sage: phi.values()
                 [-1/5, 3/2, -1/2]
                 sage: phi.dict().keys()
                 [
                 [1 0]  [ 0 -1]  [-1 -1]
                 [0 1], [ 1  3], [ 3  2]
                 ]
                 sage: phi.values() == phi.dict().values()
                 True
            """
            return list(self._map._dict.values())

        def _normalize(self):
            """
            Normalizes all of the values of the symbol self

            EXAMPLES::

                sage: E = EllipticCurve('11a')
                sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
                sage: phi = ps_modsym_from_elliptic_curve(E)
                sage: phi._normalize()
                Modular symbol with values in Sym^0 Q^2
                sage: phi._normalize().values()
                [-1/5, 3/2, -1/2]
            """
            for val in self._map:
                val.normalize()
            return self

        def __cmp__(self, other):
            """
            Checks if self == other

            EXAMPLES::

                sage: E = EllipticCurve('11a')
                sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
                sage: phi = ps_modsym_from_elliptic_curve(E)
                sage: phi == phi
                True
                sage: phi == 2*phi
                False
                sage: psi = ps_modsym_from_elliptic_curve(EllipticCurve('37a'))
                sage: psi == phi
                False
            """
            gens = self.parent._source.gens()
            for g in gens:
                c = cmp(self._map[g], other._map[g])
                if c: return c
            return 0

        def _add_(self, right):
            """
            Returns self + right

            EXAMPLES::

                sage: E = EllipticCurve('11a')
                sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
                sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
                [-1/5, 3/2, -1/2]
                sage: phi + phi
                Modular symbol with values in Sym^0 Q^2
                sage: (phi + phi).values()
                [-2/5, 3, -1]
            """
            return self.__class__(self._map + right._map, self.parent, construct=True)

        def _lmul_(self, right):
            """
            Returns self * right

            EXAMPLES::

                sage: E = EllipticCurve('11a')
                sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
                sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
                [-1/5, 3/2, -1/2]
                sage: 2*phi
                Modular symbol with values in Sym^0 Q^2
                sage: (2*phi).values()
                [-2/5, 3, -1]
            """
            return self.__class__(self._map * right, self.parent, construct=True)

        def _rmul_(self, right):
            """
            Returns self * right

            EXAMPLES::

                sage: E = EllipticCurve('11a')
                sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
                sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
                [-1/5, 3/2, -1/2]
                sage: phi*2
                Modular symbol with values in Sym^0 Q^2
                sage: (phi*2).values()
                [-2/5, 3, -1]
            """
            return self.__class__(self._map * right, self.parent, construct=True)

        def _sub_(self, right):
            """
            Returns self - right

            EXAMPLES::

                sage: E = EllipticCurve('11a')
                sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
                sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
                [-1/5, 3/2, -1/2]
                sage: phi - phi
                Modular symbol with values in Sym^0 Q^2
                sage: (phi - phi).values()
                [0, 0, 0]
            """
            return self.__class__(self._map - right._map, self.parent, construct=True)
    
        def _get_prime(self, p=None, alpha = None, allow_none=False):
            """
            Combines a prime specified by the user with the prime from the parent.

            INPUT:

            - ``p`` -- an integer or None (default None); if specified
              needs to match the prime of the parent.

            - ``alpha`` -- an element or None (default None); if p-adic
              can contribute a prime.

            - ``allow_none`` -- boolean (default False); whether to allow
              no prime to be specified.

            OUTPUT:

            - a prime or None.  If ``allow_none`` is False then a
              ValueError will be raised rather than returning None if no
              prime can be determined.

            EXAMPLES::

                sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
                sage: D = Distributions(0, 5, 10);  M = PSModularSymbols(Gamma0(2), coefficients=D)
                sage: f = M(1); f._get_prime()
                5
                sage: f._get_prime(5)
                5
                sage: f._get_prime(7)
                Traceback (most recent call last):
                ...
                ValueError: inconsistent prime
                sage: f._get_prime(alpha=Qp(5)(1))
                5
                sage: D = Symk(0);  M = PSModularSymbols(Gamma0(2), coefficients=D)
                sage: f = M(1); f._get_prime(allow_none=True) is None
                True
                sage: f._get_prime(alpha=Qp(7)(1))
                7
                sage: f._get_prime(7,alpha=Qp(7)(1))
                7
                sage: f._get_prime()
                Traceback (most recent call last):
                ...
                ValueError: you must specify a prime
            """
            pp = self.parent.prime()
            ppp = ((alpha is not None) and hasattr(alpha.parent(),'prime') and alpha.parent().prime()) or None
            p = ZZ(p) or pp or ppp
            if not p:
                if not allow_none:
                    raise ValueError("you must specify a prime")
            elif (pp and p != pp) or (ppp and p != ppp):
                raise ValueError("inconsistent prime")
            return p

        def plus_part(self):
            r"""
            Returns the plus part of self -- i.e. self + self | [1,0,0,-1].

            Note that we haven't divided by 2.  Is this a problem?

            OUTPUT:

            - self + self | [1,0,0,-1]

            EXAMPLES::

                sage: E = EllipticCurve('11a')
                sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
                sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
                [-1/5, 3/2, -1/2]
                sage: (phi.plus_part()+phi.minus_part()) == 2 * phi
                True
            """
            return self * minusproj + self

        def minus_part(self):
            r"""
            Returns the minus part of self -- i.e. self - self | [1,0,0,-1]

            Note that we haven't divided by 2.  Is this a problem?

            OUTPUT:

            - self - self | [1,0,0,-1]

            EXAMPLES::

                sage: E = EllipticCurve('11a')
                sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
                sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
                [-1/5, 3/2, -1/2]
                sage: (phi.plus_part()+phi.minus_part()) == phi * 2
                True
            """
            return self - self * minusproj

        def hecke(self, ell, algorithm="prep"):
            r"""
            Returns self | `T_{\ell}` by making use of the precomputations in
            self.prep_hecke()

            INPUT:

            - ``ell`` -- a prime

            - ``algorithm`` -- a string, either 'prep' (default) or
              'naive'

            OUTPUT:

            - The image of this element under the hecke operator
              `T_{\ell}`

            ALGORITHMS:

            - If ``algorithm == 'prep'``, precomputes a list of matrices
              that only depend on the level, then uses them to speed up
              the action.

            - If ``algorithm == 'naive'``, just acts by the matrices
              defining the Hecke operator.  That is, it computes
              sum_a self | [1,a,0,ell] + self | [ell,0,0,1],
              the last term occurring only if the level is prime to ell.

            EXAMPLES::

                sage: E = EllipticCurve('11a')
                sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
                sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
                [-1/5, 3/2, -1/2]
                sage: phi.hecke(2) == phi * E.ap(2)
                True
                sage: phi.hecke(3) == phi * E.ap(3)
                True
                sage: phi.hecke(5) == phi * E.ap(5)
                True
                sage: phi.hecke(101) == phi * E.ap(101)
                True

                sage: all([phi.hecke(p, algorithm='naive') == phi * E.ap(p) for p in [2,3,5,101]])
                True
            """
            return self.__class__(self._map.hecke(ell, algorithm), self.parent, construct=True)

        def valuation(self, p):
            r"""
            Returns the valuation of self at `p`.

            Here the valuation is the minimum of the valuations of the values of self.

            INPUT:

            - ``p`` - prime

            OUTPUT:

            - The valuation of self at `p`

            EXAMPLES::

               sage: E = EllipticCurve('11a')
               sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
               sage: phi = ps_modsym_from_elliptic_curve(E)
               sage: phi.values()
               [-1/5, 3/2, -1/2]
               sage: phi.valuation(2)
               -1
               sage: phi.valuation(3)
               0
               sage: phi.valuation(5)
               -1
               sage: phi.valuation(7)
               0
            """
            return min([val.valuation(p) for val in self._map])

        def diagonal_valuation(self, p):
            """
            Retuns the minimum of the diagonal valuation on the values of self

            INPUT:

            - ``p`` -- a positive integral prime

            EXAMPLES::

                sage: E = EllipticCurve('11a')
                sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
                sage: phi = ps_modsym_from_elliptic_curve(E)
                sage: phi.values()
                [-1/5, 3/2, -1/2]
                sage: phi.diagonal_valuation(2)
                -1
                sage: phi.diagonal_valuation(3)
                0
                sage: phi.diagonal_valuation(5)
                -1
                sage: phi.diagonal_valuation(7)
                0
            """
            return min([val.diagonal_valuation(p) for val in self._map])

