from sage.categories.all import Category

class ModsymSpace(Category):
    def __repr__(self):
        print "Category of Modular Symbol Spaces"

    class ParentMethods:
        def _repr_(self):
            pass

        def coefficient_module():
            pass
        
        def prime():
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

        def manin_relations():
            pass

        def basis(self, M=None):
            

    class ElementMethods:
        def __init__(self,map_data,parent,construct = False):
            #ModuleElement.__init__(self, parent)
            if construct:
                self._map = map_data
            else:
                self._map = ManinMap(parent._coefficients, parent._source, map_data)
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
        def normalize():
            pass
        def valuation(self, p=None):
            pass
        def scale():
            pass
        def is_zero():
            pass
        def _rmul_():
            pass
        def lift():
            pass
        def find_scalar()
        def _is_malformed()
        def act_right()

