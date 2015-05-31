from sage.structure.sage_object import SageObject
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.arith import binomial
from sage.modular.pollack_stevens.manin_map import M2Z
from sage.modular.pollack_stevens.families_util import logp_binom
from sage.all import cached_method

class padic_Lfunction(SageObject):
    pass

class padic_Lfunction_one_variable(padic_Lfunction):
    def __init__(self, Phi):
        pass

class padic_Lfunction_two_variable(padic_Lfunction):
    def __init__(self, Phis, var='T', prec=None):
        #TODO: prec: Default it to None would make us compute it.
        self._Phis = Phis    #should we create a copy of Phis, in case Phis changes? probably
        self._coefficient_ring = Phis.base_ring()
        self._base_ring = PowerSeriesRing(self._coefficient_ring, var)    #What precision?
    
    def base_ring(self):
        return self._base_ring
    
    def coefficient_ring(self):
        return self._coefficient_ring
    
    def variables(self):
        #returns (T, w)
        return (self._base_ring.gens()[0], self._coefficient_ring.gens()[0])
    
    @cached_method
    def _basic_integral(self, a, j):
        r"""
        Computes the integral
        
            .. MATH::
               
               \int_{a+p\ZZ_p}(z-\omega(a))^jd\mu.
        """
        #is this the negative of what we want?
        #if Phis is fixed for this p-adic L-function, we should make this method cached
        p = self._Phis.parent().prime()
        Da = M2Z([1,a,0,p])
        onDa = self._Phis(Da)
        aminusat = a - self._Phis.parent().base_ring().base_ring().teichmuller(a)
        #aminusat = a - self._coefficient_ring.base_ring().teichmuller(a)
        try:
            ap = self._ap
        except AttributeError:
            self._ap = self._Phis.Tq_eigenvalue(p) #catch exception if not eigensymbol
            ap = self._ap
        #print "j =", j, "a = ", a
        ans = onDa.moment(0) * (aminusat ** j)
        #print "\tr =", 0, " ans =", ans
        for r in range(1, j+1):
            ans += binomial(j, r) * (aminusat ** (j - r)) * (p ** r) * onDa.moment(r)
            #print "\tr =", r, " ans =", ans
        #print " "
        return (~ap) * ans
    
    @cached_method
    def _compute_nth_coeff(self, n, twist=None):
        r"""
        Computes the coefficient of T^n.
        """
        #TODO: Check that n is not too big
        #TODO implement twist
        p = self._Phis.parent().prime()
        prec = self._Phis.precision_absolute()[0] #Not quite right, probably
        #print "@@@@n =", n, "prec =", prec
        cjns = list(logp_binom(n, p, prec+1))
        teich = self._Phis.parent().base_ring().base_ring().teichmuller
        return sum([cjns[j] * sum([((~teich(a)) ** j) * self._basic_integral(a, j) for a in range(1,p)]) for j in range(min(prec,len(cjns)))])
    
    def coefficient(self, index, twist=None):
        r"""
        index should be some sort of pair (i,j) corresponding to T^i w^j. Maybe if
        index is just one number, i, it can return the coefficient of T^i to biggest
        possible precision in w.
        
        Should also make it so that one can pass some monomial in the variables.
        """
        pass
    
    def power_series(self, twist=None, prec=None):
        r"""
        returns a power series in base_ring, up to given precision, or max prec if prec is None
        """
        p = self._Phis.parent().prime()
        prec = self._Phis.precision_absolute()[0] #Not quite right, probably
        return self._base_ring([self._compute_nth_coeff(n, twist) for n in range(prec)])
