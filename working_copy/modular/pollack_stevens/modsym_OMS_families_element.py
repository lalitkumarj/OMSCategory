from sage.misc.cachefunc import cached_method
from sage.misc.misc import verbose
from sage.modular.pollack_stevens.modsym_element import ModularSymbolElement_generic
#from sage.modular.pollack_stevens.coeffmod_OMS_families_element import _padic_val_of_pow_series

class ModSym_OMS_Families_element(ModularSymbolElement_generic):
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
    
    @cached_method
    def is_Tq_eigensymbol(self,q,p=None,M=None):
        r"""
        Determines if self is an eigenvector for `T_q` to precision ``M``.
        """
        
        try:
            aq = self.Tq_eigenvalue(q, p, M)
            return True
        except ValueError:
            return False
    
    @cached_method
    def Tq_eigenvalue(self, q, p=None, M=None, check=True):
        #TODO: deal with var_prec
        qhecke = self.hecke(q)
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
        if check:
            verbose("Checking that this is actually an eigensymbol")
            if p is None or M is None:
                for g in gens[1:]:
                    #Hack
                    if qhecke._map[g] != R(aq) * self._map[g]:
                    #if qhecke._map[g] != aq * self._map[g]:
                        raise ValueError("not a scalar multiple")
            #Hack
            elif (M is not None and qhecke - R(aq) * self).valuation(p) < M[0]:            
            #elif (M is not None and qhecke - aq * self).valuation(p) < M[0]:
                raise ValueError("not a scalar multiple")
        return aq