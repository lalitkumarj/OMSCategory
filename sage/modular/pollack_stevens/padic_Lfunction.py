from sage.structure.sage_object import SageObject
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.arith import binomial
from sage.modular.pollack_stevens.manin_map import M2Z
from sage.modular.pollack_stevens.families_util import logp_binom
from sage.rings.integer_ring import ZZ
from sage.rings.big_oh import O
from sage.modular.dirichlet import kronecker_character
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
        self._prec = prec
    
    def base_ring(self):
        return self._base_ring
    
    def coefficient_ring(self):
        return self._coefficient_ring
    
    def variables(self):
        #returns (T, w)
        return (self._base_ring.gens()[0], self._coefficient_ring.gens()[0])
    
    def _max_coeff(self):
        Phis = self._Phis
        p = Phis.parent().prime()
        p_prec, var_prec = Phis.precision_absolute()
        max_j = Phis.parent().coefficient_module().length_of_moments(p_prec)
        n = 0
        while True:
            if max_j - (n / (p-1)).floor() - min(max_j, n) - (max_j / p).floor() <= 0:
                return n - 1
            n += 1
    
    def _on_Da(self, a, twist):
        r"""
        An internal method used by ``self._basic_integral``. The parameter ``twist`` is assumed to be either ``None`` or a primitive quaratic Dirichlet character.
        """
        p = self._Phis.parent().prime()
        if twist is None:
            return self._Phis(M2Z([1,a,0,p]))
        D = twist.level()
        DD = self._Phis.parent().coefficient_module()
        S0 = DD.action().actor()
        ans = DD.zero()
        for b in range(1, D + 1):
            if D.gcd(b) == 1:
                ans += twist(b) * (self._Phis(M2Z([1, D * a + b * p, 0, D * p])) * S0([1, b / D, 0, 1]))
        return ans
    
    @cached_method
    def _basic_integral(self, a, j, twist=None):
        r"""
        Computes the integral
        
            .. MATH::
               
               \int_{a+p\ZZ_p}(z-\omega(a))^jd\mu_\chi.
        
        If ``twist`` is ``None``, `\\chi` is the trivial character. Otherwise, ``twist`` can be a primitive quadratic character of conductor prime to `p`.
        """
        #is this the negative of what we want?
        #if Phis is fixed for this p-adic L-function, we should make this method cached
        p = self._Phis.parent().prime()
        if twist is None:
            pass
        elif twist in ZZ:
            twist = kronecker_character(twist)
            if twist.is_trivial():
                twist = None
            else:
                D = twist.level()
                assert(D.gcd(p) == 1)
        else:
            if twist.is_trivial():
                twist = None
            else:
                assert((twist**2).is_trivial())
                twist = twist.primitive_character()
                D = twist.level()
                assert(D.gcd(p) == 1)
        
        onDa = self._on_Da(a, twist)#self._Phis(Da)
        aminusat = a - self._Phis.parent().base_ring().base_ring().teichmuller(a)
        #aminusat = a - self._coefficient_ring.base_ring().teichmuller(a)
        try:
            ap = self._ap
        except AttributeError:
            self._ap = self._Phis.Tq_eigenvalue(p) #catch exception if not eigensymbol
            ap = self._ap
        if not twist is None:
            ap *= twist(p)
        if j == 0:
            return (~ap) * onDa.moment(0)
        if a == 1:
            #aminusat is 0, so only the j=r term is non-zero
            return (~ap) * (p ** j) * onDa.moment(j)
        #print "j =", j, "a = ", a
        ans = onDa.moment(0) * (aminusat ** j)
        #ans = onDa.moment(0)
        #print "\tr =", 0, " ans =", ans
        for r in range(1, j+1):
            if r == j:
                ans += binomial(j, r) * (p ** r) * onDa.moment(r)
            else:
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
        Phis = self._Phis
        p = Phis.parent().prime()
        if n == 0:
            return sum([self._basic_integral(a, 0, twist) for a in range(1, p)])
        p_prec, var_prec = Phis.precision_absolute()
        max_j = Phis.parent().coefficient_module().length_of_moments(p_prec)
        ans_prec = max_j - (n / (p-1)).floor() - min(max_j, n) - (max_j / p).floor()
        if ans_prec == 0:
            return self._coefficient_ring(0)
        #prec = self._Phis.parent()#precision_absolute()[0] #Not quite right, probably
        #print "@@@@n =", n, "prec =", prec
        cjns = list(logp_binom(n, p, max_j+1))
        #print cjns
        teich = Phis.parent().base_ring().base_ring().teichmuller
        #Next line should work but loses precision!!!
        ans = sum([cjns[j] * sum([((~teich(a)) ** j) * self._basic_integral(a, j, twist) for a in range(1,p)]) for j in range(1, min(max_j, len(cjns)))])
        #Instead do this messed up thing
        w = ans.parent().gen()
        #ans = 0*w
        #for j in range(1,min(max_j, len(cjns))):
        #    ans_term = [0*w] * var_prec
        #    for a in range(1,p):
        #        term = (((~teich(a)) ** j) * self._basic_integral(a, j, twist)).list()
        #        for i in range(min(var_prec, len(term))):
        #            ans_term[i] += term[i]
        #    ans += cjns[j] * sum([ans_term[i] * w**i for i in range(var_prec)])
        #print ans_prec
        ans_prec = O(p**ans_prec)
        #print ans_prec
        #print ans
        for i in range(ans.degree() + 1):
            ans += ans_prec * w**i
        return ans
    
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
        if prec is None:
            prec = self._max_coeff()#Phis.precision_absolute()[0] #Not quite right, probably
        else:
            pass #do some checks on inputted prec
        return self._base_ring([self._compute_nth_coeff(n, twist) for n in range(prec + 1)])
