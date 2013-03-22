from sage.misc.cachefunc import cached_method
from sage.misc.misc import verbose
from sage.modular.pollack_stevens.modsym_element import ModularSymbolElement_generic

class ModSym_OMS_element(ModularSymbolElement_generic):
    def is_zero(self):
        z = True
        j = 0
        while z and (j<len(self.values())-1):
            if (self.values()[j].is_zero()):
                j += 1
            else:
                z = False
        return z
                        
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

    def list_of_total_measures(self):
        r"""
        Returns the list of total measures of the OMS evaluated at each of our generators

        INPUT:

        - None
        
        OUTPUT:

        - List of p-adic numberes
        """
        return [mu.moment(0) for mu in self.values()]
    
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
            elif (qhecke - aq * self).valuation(p) < M[0]:
                raise ValueError("not a scalar multiple")
        return aq

    def _consistency_check(self):
        """
            Check that the map really does satisfy the Manin relations loop (for debugging).
            The two and three torsion relations are checked and it is checked that the symbol
            adds up correctly around the fundamental domain

            EXAMPLES::

                sage: D = OverconvergentDistributions(0, 7, base=Qp(7,5))
                sage: MS = OverconvergentModularSymbols(14, coefficients=D);
                sage: MR = MS.source()
                sage: V = D.approx_module()
                sage: Phi_dict = {MR.gens()[0]:7 * V((5 + 6*7 + 7^2 + 4*7^3 + 5*7^4 + O(7^5), 6 + 5*7 + 2*7^2 + 3*7^3 + O(7^4), 4 + 5*7 + 7^2 + O(7^3), 2 + 7 + O(7^2), 4 + O(7))), MR.gens()[1]:7 * V((4 + 2*7 + 4*7^2 + 5*7^3 + O(7^5), 5 + 7^2 + 4*7^3 + O(7^4), 1 + 7 + 5*7^2 + O(7^3), 2 + O(7^2), 4 + O(7))), MR.gens()[2]:7 * V((3 + 6*7 + 6*7^2 + 5*7^3 + 7^4 + O(7^5), 6 + 5*7 + 5*7^3 + O(7^4), 3 + 6*7^2 + O(7^3), 6 + 2*7 + O(7^2), O(7))), MR.gens()[3]:7 * V((5 + 3*7 + 4*7^2 + 7^3 + 3*7^4 + O(7^5), 2 + 4*7^2 + 2*7^3 + O(7^4), 1 + 4*7 + 2*7^2 + O(7^3), 6*7 + O(7^2), 6 + O(7))), MR.gens()[4]:7 * V((3 + 2*7^2 + 3*7^3 + 3*7^4 + O(7^5), 5*7 + 4*7^2 + 2*7^3 + O(7^4), 6 + 4*7 + 2*7^2 + O(7^3), 2 + 3*7 + O(7^2), O(7)))}
                sage: Phi = MS(Phi_dict)
                sage: Phi._consistency_check()
                This modular symbol satisfies the manin relations
        """
#            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
#            sage: E = EllipticCurve('37a1')
#            sage: phi = ps_modsym_from_elliptic_curve(E)
#            sage: phi._consistency_check()
#            This modular symbol satisfies the manin relations
        f = self._map
        MR = self._map._manin
        ## Test two torsion relations
        for g in MR.reps_with_two_torsion():
            gamg = MR.two_torsion_matrix(g)
            if not (f[g]*gamg + f[g]).is_zero():
                raise ValueError("Two torsion relation failed with",g)

        ## Test three torsion relations
        for g in MR.reps_with_three_torsion():
            gamg = MR.three_torsion_matrix(g)
            if not (f[g]*(gamg**2) + f[g]*gamg + f[g]).is_zero():
                raise ValueError("Three torsion relation failed with",g)

        ## Test that the symbol adds to 0 around the boundary of the fundamental domain
        t = self.parent().coefficient_module().zero_element()
        for g in MR.gens()[1:]:
            if (not g in MR.reps_with_two_torsion()) and (not g in MR.reps_with_three_torsion()):
                t += f[g] * MR.gammas[g] - f[g]
            else:
                if g in MR.reps_with_two_torsion():
                    t -= f[g] 
                else:
                    t -= f[g]
                    
        id = MR.gens()[0]
        if f[id]*MR.gammas[id] - f[id] != -t:
            print -t
            print f[id]*MR.gammas[id] - f[id]
            raise ValueError("Does not add up correctly around loop")

        print "This modular symbol satisfies the manin relations"

            
