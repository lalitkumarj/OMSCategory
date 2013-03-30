from sage.misc.cachefunc import cached_method
from sage.misc.misc import verbose
from sage.modular.pollack_stevens.modsym_element import ModularSymbolElement_generic

class ModSym_OMS_element(ModularSymbolElement_generic):                        
    def valuation(self, p=None):
        if p is None:
            p = self.parent().prime()
        elif p != self.parent().prime():
            raise ValueError("Specified prime(=%s) must match prime of base ring(=%s)"%(p, self.parent().prime()))
        return min([val.valuation() for val in self._map])
    
    def diagonal_valuation(self, p):
        if p is None:
            p = self.parent().prime()
        elif p != self.parent().prime():
            raise ValueError("Specified prime(=%s) must match prime of base ring(=%s)"%(p, self.parent().prime()))
        return min([val.diagonal_valuation() for val in self._map])
    
    def precision_relative(self):
        return min([mu.precision_relative() for mu in self.values()])

    def list_of_total_measures(self):
        r"""
        Returns the list of total measures of the OMS evaluated at each of our generators

        INPUT:

        - None
        
        OUTPUT:

        - List of p-adic numberes
        """
        z = self.parent().base().zero()
        return [mu.moment(0) if mu.precision_relative() != 0 else z for mu in self.values()]
    
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
        """
#        EXAMPLES::
#
#            sage: E = EllipticCurve('11a')
#            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
#            sage: phi = ps_modsym_from_elliptic_curve(E)
#            sage: phi.values()
#            [-1/5, 3/2, -1/2]
#            sage: phi_ord = phi.p_stabilize(p = 3, ap = E.ap(3), M = 10, ordinary = True)
#            sage: phi_ord.is_Tq_eigensymbol(2,3,10)
#            True
#            sage: phi_ord.is_Tq_eigensymbol(2,3,100)
#            False
#            sage: phi_ord.is_Tq_eigensymbol(2,3,1000)
#            False
#            sage: phi_ord.is_Tq_eigensymbol(3,3,10)
#            True
#            sage: phi_ord.is_Tq_eigensymbol(3,3,100)
#            False
#        """
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
        """
#        EXAMPLES::
#
#            sage: E = EllipticCurve('11a')
#            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
#            sage: phi = ps_modsym_from_elliptic_curve(E)
#            sage: phi.values()
#            [-1/5, 3/2, -1/2]
#            sage: phi_ord = phi.p_stabilize(p = 3, ap = E.ap(3), M = 10, ordinary = True)
#            sage: phi_ord.Tq_eigenvalue(2,3,10) + 2
#            O(3^10)
#
#            sage: phi_ord.Tq_eigenvalue(3,3,10)
#            2 + 3^2 + 2*3^3 + 2*3^4 + 2*3^6 + 3^8 + 2*3^9 + O(3^10)
#            sage: phi_ord.Tq_eigenvalue(3,3,100)
#            Traceback (most recent call last):
#            ...
#            ValueError: not a scalar multiple
#        """
        qhecke = self.hecke(q)
        gens = self.parent().source().gens()
        if p is None:
            p = self.parent().prime()   #CHANGE THIS
        i = 0
        g = gens[i]
        verbose("Computing eigenvalue")
        while self._map[g].is_zero(M):
            if not qhecke._map[g].is_zero(M):
                raise ValueError("not a scalar multiple")
            i += 1
            try:
                g = gens[i]
            except IndexError:
                raise ValueError("self is zero")
        aq = self._map[g].find_scalar(qhecke._map[g], M, check)
        verbose("Found eigenvalues of %s"%(aq))
        if check:
            verbose("Checking that this is actually an eigensymbol")
            if p is None or M is None:
                for g in gens[1:]:
                    if qhecke._map[g] != aq * self._map[g]:
                        raise ValueError("not a scalar multiple")
            elif (qhecke - aq * self).valuation() < M[0]:
                raise ValueError("not a scalar multiple")
        return aq
