import operator

from sage.misc.cachefunc import cached_method
from sage.misc.misc import verbose
from sage.categories.action import Action
from sage.structure.element import ModuleElement

from sage.modular.pollack_stevens.manin_map import ManinMap

minusproj = [1,0,0,-1]

class ModSymAction(Action):
    def __init__(self, actor, MSspace):
        Action.__init__(self, actor, MSspace, False, operator.mul)

    def _call_(self, symb, g):
        return symb.__class__(symb._map * g, symb.parent(), construct=True)

class ModularSymbolElement_generic(ModuleElement):
    def __init__(self, map_data, parent, construct = False):
        ModuleElement.__init__(self, parent)
        if construct:
            self._map = map_data
        else:
            self._map = ManinMap(parent._coefficients, parent._source, map_data)
    
    def _repr_(self):
        return "Modular symbol of level %s with values in %s"%(self.parent().level(), self.parent().coefficient_module())
    
    def dict(self):
        D = {}
        for g in self.parent().source().gens():
            D[g] = self._map[g]
        return D
    
    def weight(self):
        return self.parent().weight()
    
    def values(self):
        return [self._map[g] for g in self.parent().source().gens()]
    
    def _normalize(self):
        for val in self._map:
            val.normalize()
        return self
    
    def __cmp__(self, other):
        gens = self.parent().source().gens()
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
            Modular symbol of level 11 with values in Sym^0 Q^2
            sage: (phi + phi).values()
            [-2/5, 3, -1]
        """
        return self.__class__(self._map + right._map, self.parent(), construct=True)

    def _lmul_(self, right):
        """
        Returns self * right

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: 2*phi
            Modular symbol of level 11 with values in Sym^0 Q^2
            sage: (2*phi).values()
            [-2/5, 3, -1]
        """
        return self.__class__(self._map * right, self.parent(), construct=True)

    def _rmul_(self, right):
        """
        Returns self * right

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi*2
            Modular symbol of level 11 with values in Sym^0 Q^2
            sage: (phi*2).values()
            [-2/5, 3, -1]
        """
        return self.__class__(self._map * right, self.parent(), construct=True)

    def _sub_(self, right):
        """
        Returns self - right

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi - phi
            Modular symbol of level 11 with values in Sym^0 Q^2
            sage: (phi - phi).values()
            [0, 0, 0]
        """
        return self.__class__(self._map - right._map, self.parent(), construct=True)
    
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
        S0N = Sigma0(self.parent().level())
        return self + self * S0N(minusproj)

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
        S0N = Sigma0(self.parent().level())
        return self - self * S0N(minusproj)
    
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
        return self.__class__(self._map.hecke(ell, algorithm), self.parent(), construct=True)
    
    def valuation(self, p):
        if p is None:
            raise ValueError("Must specify p.")
        return min([val.valuation(q) for val in self._map])
    
    def diagonal_valuation(self, p):
        if p is None:
            raise ValueError("Must specify p.")
        return min([val.diagonal_valuation(p) for val in self._map])
    
    @cached_method
    def is_Tq_eigensymbol(self,q,p=None,M=None):
        r"""
        Determines if self is an eigenvector for `T_q` modulo `p^M`

        INPUT:

        - ``q`` -- prime of the Hecke operator

        - ``p`` -- prime we are working modulo

        - ``M`` -- degree of accuracy of approximation

        OUTPUT:

        - True/False

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi_ord = phi.p_stabilize(p = 3, ap = E.ap(3), M = 10, ordinary = True)
            sage: phi_ord.is_Tq_eigensymbol(2,3,10)
            True
            sage: phi_ord.is_Tq_eigensymbol(2,3,100)
            False
            sage: phi_ord.is_Tq_eigensymbol(2,3,1000)
            False
            sage: phi_ord.is_Tq_eigensymbol(3,3,10)
            True
            sage: phi_ord.is_Tq_eigensymbol(3,3,100)
            False
        """
        try:
            aq = self.Tq_eigenvalue(q, p, M)
            return True
        except ValueError:
            return False

    # what happens if a cached method raises an error?  Is it recomputed each time?
    @cached_method
    def Tq_eigenvalue(self, q, p=None, M=None, check=True):
        r"""
        Eigenvalue of `T_q` modulo `p^M`

        INPUT:

        - ``q`` -- prime of the Hecke operator

        - ``p`` -- prime we are working modulo (default: None)

        - ``M`` -- degree of accuracy of approximation (default: None)

        - ``check`` -- check that `self` is an eigensymbol

        OUTPUT:

        - Constant `c` such that `self|T_q - c * self` has valuation greater than
          or equal to `M` (if it exists), otherwise raises ValueError

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi_ord = phi.p_stabilize(p = 3, ap = E.ap(3), M = 10, ordinary = True)
            sage: phi_ord.Tq_eigenvalue(2,3,10) + 2
            O(3^10)

            sage: phi_ord.Tq_eigenvalue(3,3,10)
            2 + 3^2 + 2*3^3 + 2*3^4 + 2*3^6 + 3^8 + 2*3^9 + O(3^10)
            sage: phi_ord.Tq_eigenvalue(3,3,100)
            Traceback (most recent call last):
            ...
            ValueError: not a scalar multiple
        """
        qhecke = self.hecke(q)
        gens = self.parent().source().gens()
        if p is None:
            p = self.parent().prime()   #CHANGE THIS
        i = 0
        g = gens[i]
        verbose("Computing eigenvalue")
        while self._map[g].is_zero(p, M):
            if not qhecke._map[g].is_zero(p, M):
                raise ValueError("not a scalar multiple")
            i += 1
            try:
                g = gens[i]
            except IndexError:
                raise ValueError("self is zero")
        aq = self._map[g].find_scalar(qhecke._map[g], p, M, check)
        verbose("Found eigenvalues of %s"%(aq))
        if check:
            verbose("Checking that this is actually an eigensymbol")
            if p is None or M is None:
                for g in gens[1:]:
                    if qhecke._map[g] != aq * self._map[g]:
                        raise ValueError("not a scalar multiple")
            elif (qhecke - aq * self).valuation(p) < M:
                raise ValueError("not a scalar multiple")
        return aq
    
    def is_ordinary(self, p=None):
        raise NotImplementedError
    
    def reduce_precision(self, M):
        r"""
        Only holds on to `M` moments of each value of self
        """
        return self.__class__(self._map.reduce_precision(M), self.parent(), construct=True)

    def precision_absolute(self):
        r"""
        Returns the number of moments of each value of self
        """
        return min([a.precision_absolute() for a in self._map])
    
    def normalize(self):
        for val in self._map:
            val.normalize()
        return self
    
    #def __call__(self):
    #    pass
    
    