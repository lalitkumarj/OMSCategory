from sage.misc.cachefunc import cached_method
from sage.misc.misc import verbose
from sage.modular.pollack_stevens.modsym_element import ModularSymbolElement_generic
from sage.modular.pollack_stevens.coeffmod_OMS_families_element import _add_big_ohs_list#, _padic_val_of_pow_series
from sage.modular.pollack_stevens.padic_Lfunction import padic_Lfunction_two_variable

class ModSym_OMS_Families_element(ModularSymbolElement_generic):
    """
    TESTS::
    
        sage: D = FamiliesOfOverconvergentDistributions(0, prec_cap=[5,3], base_coeffs=ZpCA(11))
        sage: MS = FamiliesOfOMS(11, coefficients=D)
        sage: Phi = MS.random_element()
        sage: TestSuite(Phi).run()
    """
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
        Returns the list of total measures of the OMS families evaluated at each of our generators.

        INPUT:

        - None
        
        OUTPUT:

        - List of p-adic numbers
        """
        z = self.parent().base().zero()
        p = self.parent().prime()
        return [mu.moment(0) if mu.precision_relative()[0] != 0 else z for mu in self.values()]

    def list_of_total_measures_at_fixed_weight(self,k=None):
        r"""
        Returns the list of total measures of the OMS specialized to weight k and evaluated at each of our generators.
        If no weight is given, just use the weight from the disc the family lives on.

        INPUT:

        - k -- (optional) integer
        
        OUTPUT:

        - List of p-adic numbers
        """
        if k == None:
            k = self.weight()
        z = self.parent().base().zero()
        p = self.parent().prime()
        return [mu.moment(0).polynomial().substitute(w=((1+p)**k-1)/p) if mu.precision_relative()[0] != 0 else z for mu in self.values()]
    
    #@cached_method
    def is_Tq_eigensymbol(self,q,p=None,M=None):
        r"""
        Determines if self is an eigenvector for `T_q` to precision ``M``.
        """
        
        try:
            aq = self.Tq_eigenvalue(q, p, M)
            return True
        except ValueError:
            return False
    
    #@cached_method
    def Tq_eigenvalue(self, q, p=None, M=None, check=True):
        #TODO: deal with var_prec
        if self == 0:
            raise ValueError("0 cannot be an eigenvector")
        qhecke = self.hecke(q)
        if True:
            rel = self.parent().linear_relation([self], qhecke)
            if rel[-1] is None:
                raise ValueError("not a scalar multiple")
            return -rel[0][0] / rel[1]
        gens = self.parent().source().gens()
        if p is None:
            p = self.parent().prime()   #CHANGE THIS
        i = 0
        g = gens[i]
        verbose("Computing eigenvalue to precision %s"%(M))
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
        R = self.parent().base_ring()   #needed for hack below
        Rbase = R.base_ring()
        RK = R.change_ring(Rbase.fraction_field())
        if check:
            verbose("Checking that this is actually an eigensymbol")
            if p is None or M is None:
                for g in gens[1:]:
                    #Hack
                    diff = qhecke._map[g] - R(aq) * self._map[g]
                    if diff != 0:
                    #if qhecke._map[g] != aq * self._map[g]:
                        val = diff.valuation()
                        if val < 1:
                            raise ValueError("not a scalar multiple. difference is %s"%(qhecke._map[g] - R(aq) * self._map[g]))
                        verbose("Resetting aq (lowering p-adic precision to %s)"%(val))
                        aq = RK(_add_big_ohs_list(aq, [val, aq.prec()]), aq.prec())
            #Hack
            elif (M is not None and qhecke - R(aq) * self).valuation(p) < M[0]:            
            #elif (M is not None and qhecke - aq * self).valuation(p) < M[0]:
                raise ValueError("not a scalar multiple")
        return aq
    
    def padic_lfunction(self, var='T', prec=None):
        return padic_Lfunction_two_variable(self, var=var, prec=prec)
