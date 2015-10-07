import sage.rings.ring as ring

from sage.structure.factory import UniqueFactory
from sage.rings.integer_ring import ZZ
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.padics.factory import ZpCA, Qp
from sage.misc.cachefunc import cached_method

from sage.modular.pollack_stevens.sigma0 import _default_adjuster
from sage.modular.pollack_stevens.coeffmod_space import CoefficientModule_generic
from sage.modular.pollack_stevens.coeffmod_element import WeightKAction_OMS_fam
from sage.modular.pollack_stevens.coeffmod_OMS_element import CoeffMod_OMS_element
from sage.modular.pollack_stevens.coeffmod_OMS_families_element import CoeffMod_OMS_Families_element

def _prec_cap_parser(prec_cap):
    r"""
    The precision of a family is an ordered pair (p_prec, var_prec) where p_prec
    is the `p`-adic precision cap and var_prec is the precision in the variable
    of the power series ring. This function takes an input ``prec_cap`` and outputs
    an ordered pair. 
    
    EXAMPLES::
    
        sage: from sage.modular.pollack_stevens.coeffmod_OMS_families_space import _prec_cap_parser
        sage: _prec_cap_parser(5)
        [5, None]
        sage: _prec_cap_parser([20, 10])
        [20, 10]
        sage: _prec_cap_parser([7])
        [7, None]
        sage: _prec_cap_parser(None)
        Traceback (most recent call last):
        ...
        ValueError: The case prec_cap is None should be dealt with prior to calling this function.
        sage: _prec_cap_parser([5, 0])
        [5, 0]
    """
    if prec_cap is None:
        raise ValueError("The case prec_cap is None should be dealt with prior to calling this function.")
    if not isinstance(prec_cap, (list, tuple)):
        prec = ZZ(prec_cap)
        if prec < 0:
            ValueError("prec_cap should be at least 0.")
        return [prec, None]
    elif len(prec_cap) == 1:
        prec = ZZ(prec_cap[0])
        if prec < 0:
            ValueError("prec_cap should be at least 0.")
        return [prec, None]
    elif len(prec_cap) > 2:
        raise ValueError("prec_cap should not have length > 2.")
    else:
#        if p_prec < 1 or var_prec < 1:
#            raise ValueError("Precisions should be at least 1.")
        return [ZZ(prec_cap[0]), ZZ(prec_cap[1])]

#def _factor_Dir_char(chi, p):
#    r"""
#    Given a Dirichlet character `\chi` of conductor exactly divisible by the prime
#    `p`, returns a pair 
#    
#    `(\chi_p, \chi^\prime)` such that `\chi=\chi_p\chi^\prime`,
#    `\chi_p` has level dividing `p`, and `\chi^\prime` has conductor prime to `p`
#    (but
#    
#    EXAMPLES::
#    
#        sage: chi = DirichletGroup(3*5*49)
#    """
#    #FINISH THIS
#    if p == 2:
#        raise NotImplementedError
#    DG = chi.parent()
#    unit_gens = DG.unit_gens()
#    for u in unit_gens:
#        if u % p != 1:
#            break
#    return

class CoeffMod_OMS_Families_factory(UniqueFactory):
    """
    EXAMPLES::
    
        sage: D = FamiliesOfOverconvergentDistributions(0, prec_cap=[5,3], base_coeffs=ZpCA(7)); D    # indirect doctest
        Families of overconvergent distributions on the disc 0 over Power Series Ring in w over 7-adic Ring with capped absolute precision 20
        sage: D = FamiliesOfOverconvergentDistributions(2, prec_cap = [8 ,4], base_coeffs=ZpCA(11, 4))
        Traceback (most recent call last):
        ...
        ValueError: Precision cap on coefficients of base ring must be at least the p-adic precision cap of this space.
    """
    def create_key(self, k, p=None, prec_cap=None, base=None, base_coeffs=None, \
                     character=None, adjuster=None, act_on_left=False, \
                     dettwist=None, variable_name = 'w'):
        if base is None:
            if base_coeffs is None:
                if p is None:
                    raise ValueError("Must specify a prime, a base ring, or coefficients.")
                if prec_cap is None:
                    raise ValueError("Must specify a precision cap or a base ring.")
                prec_cap = _prec_cap_parser(prec_cap)
                base_coeffs = Zp(p, prec = prec_cap[0])
            elif not isinstance(base_coeffs, ring.Ring):
                raise TypeError("base_coeffs must be a ring if base is None")
            elif prec_cap is None:
                raise ValueError("Must specify a precision cap or a base ring.")
            else:
                prec_cap = _prec_cap_parser(prec_cap)
                if base_coeffs.precision_cap() < prec_cap[0]:
                    raise ValueError("Precision cap on coefficients of base ring must be at least the p-adic precision cap of this space.")
            base = PowerSeriesRing(base_coeffs, name=variable_name, default_prec=prec_cap[1])
        elif prec_cap is None:
            prec_cap = [ZZ(base.base_ring().precision_cap()), ZZ(base.default_prec())]
        else:
            if base.base_ring().precision_cap() < prec_cap[0]:
                raise ValueError("Precision cap on coefficients of base ring must be at least the p-adic precision cap of this space.")
            if base.default_prec() < prec_cap[1]:
                raise ValueError("Default precision on the variable of base ring must be at least the w-adic precision cap of this space.")
        base_coeffs = None
        p = base.base_ring().prime()
        k_shift = 0
        #if character is not None:
        #    #Should we ensure character is primitive?
        #    cond_val = character.conductor().valuation(p)
        #    if cond_val > 1:
        #        raise ValueError("Level must not be divisible by p^2.")
        #    elif cond_val == 1:
        #        pass
        #    else:
        #        pass
        k = ZZ((k + k_shift) % (p-1))
        #if prec_cap is None:
        #    prec_cap = [base.base_ring().precision_cap, base.default_prec()]
        #else:
        #    prec_cap = list(prec_cap)
        #prec_cap = [base.base_ring().precision_cap(), base.default_prec()]
        if adjuster is None:
            adjuster = _default_adjuster()
        if dettwist is not None:
            dettwist = ZZ(dettwist)
            if dettwist == 0: 
                dettwist = None
        return (k, p, tuple(prec_cap), base, character, adjuster, act_on_left, dettwist)
    
    def create_object(self, version, key):
        k, p, prec_cap, base, character, adjuster, act_on_left, dettwist = key
        return CoeffMod_OMS_Families_space(k, p=p, prec_cap=prec_cap, base=base, character=character, \
                 adjuster=adjuster, act_on_left=act_on_left, dettwist=dettwist)

FamiliesOfOverconvergentDistributions = CoeffMod_OMS_Families_factory('FamiliesOfOverconvergentDistributions')

class CoeffMod_OMS_Families_space(CoefficientModule_generic):
    r"""
    EXAMPLES::
    
        sage: from sage.modular.pollack_stevens.coeffmod_OMS_families_space import FamiliesOfOverconvergentDistributions
        sage: D = FamiliesOfOverconvergentDistributions(0, base=PowerSeriesRing(Zp(3, 10), 'w', default_prec=5)); D
        Families of overconvergent distributions on the disc 0 over Power Series Ring in w over 3-adic Ring with capped relative precision 10
        sage: D.prime()
        3
        sage: D.precision_cap()
        (10, 5)
    
    TEST::
    
        sage: TestSuite(D).run()
    """
    
    def __init__(self, k, p=None, prec_cap=[20, 10], base=None, base_coeffs=None, \
                 character=None, adjuster=None, act_on_left=False, \
                 dettwist=None, action_class = WeightKAction_OMS_fam, \
                 variable_name = 'w'):
        #self._prec_cap = prec_cap
        if base is None:
            if base_coeffs is None:
                base_coeffs = Zp(p, prec = prec_cap[0])
            elif not isinstance(base_coeffs, ring.Ring):
                raise TypeError("base_coeffs must be a ring if base is None")
            base = PowerSeriesRing(base_coeffs, name=variable_name)
        #elif not isinstance(base, ring.Ring):
        #    raise TypeError("base must be a ring")
        self._p = base.base_ring().prime()
        #self._prec_cap = [base.base_ring().precision_cap(), base.default_prec()]
        self._prec_cap = tuple(prec_cap)
        k = k % (self._p - 1)
        #self._cp = (self._p-2) / (self._p-1)
        CoefficientModule_generic.__init__(self, k, base=base, \
                 character=character, adjuster=adjuster, act_on_left=act_on_left, \
                 dettwist=dettwist, action_class=action_class, \
                 element_class=CoeffMod_OMS_Families_element, padic=True)
    
    #def __reduce__(self):
    #    return FamiliesOfOverconvergentDistributions.reduce_data(self)
    
    def _coerce_map_from_(self, other):
        if isinstance(other, CoeffMod_OMS_Families_element) \
            and other._k  == self._k \
            and self._character == other._character \
            and self.base_ring().has_coerce_map_from(other.base_ring()):
            #Should we check data related to action?
            return True
        elif isinstance(other, CoeffMod_OMS_element):
            #ADD THIS?
            #kdiff = other._k - self._k
            #self._character = other._character
            #and self.base_ring().has_coerce_map_from(other.base_ring()):
            return False
        else:
            return False
    
    def _repr_(self):
        s = "Families of overconvergent distributions on the disc %s " % (self._k)
        if self._character is not None:
            s += "twisted by %s " % (self._character)
        if self._dettwist is not None:
            s += "and det^%s " % (self._dettwist)
        s += "over %s" % (self.base_ring())
        return s
        
    def weight(self):
        r"""
        Returns the non-negative integer less than or equal to `p-1` that encodes
        the disc of weight space on which this family lives.
        """
        return self._k
    
    def prime(self):
        return self._p
    
    @cached_method
    def filtration_precisions(self, M=None):
#        r"""
#        Returns a tuple whose `i`th entry is `M-\lfloor i(p-2)/(p-1)\rfloor`,
#        the `p`-adic precision of the `i`th moment of an element of this space.
#        If ``M`` is None, uses the precision cap of this space.
#        """
        if M is None:
            M = self.precision_cap()[0]
        elif M < 0:
            raise ValueError("M(=%s) must be at least 0."%(M))
        if M == 0:
            return ()
        if M == 1:
            return (1, 1)
        pm1 = self._p - 1
        a = M
        i = 0
        precs = []
        while True:
            precs.append(a)
            if i % pm1 != 0:
                a -= 1
            if a == 0:
                break
            i += 1
        return tuple(precs)
        #return tuple((M * self._cp).ceil() - (i * self._cp).floor() for i in range(M))
        #return tuple(((M-i) * self._cp).ceil() for i in range(M))
    
    @cached_method
    def length_of_moments(self, M=None):
        if M is None:
            M = self.precision_cap()[0]
        elif M < 0:
            raise ValueError("M(=%s) must be at least 0."%(M))
        if M == 0:
            return ZZ(0)
        if M == 1:
            return ZZ(2)
        M = ZZ(M)
        return M + (M/(self._p - 2)).ceil()
    
    @cached_method
    def length_reverse_lookup(self, L):
        r"""
        Given a number of moments ``L``, return ``M`` such that a precision of
        ``M`` will have the largest number of moments less than or equal to ``L``.
        """
        if L < 0:
            raise ValueError("L(=%s) must be at least 0."%(L))
        if L <= 1:
            return ZZ(0)
        L = ZZ(L)
        return ZZ(L - 1) - ((L - 1) / (self._p - 1)).floor()
    
    def precision_cap(self):
        return self._prec_cap
    
    def change_ring(self, new_base_ring):
        if new_base_ring == self.base_ring():
            return self
        return FamiliesOfOverconvergentDistributions(self._k, prec_cap=self._prec_cap, base=new_base_ring, character=self._character, adjuster=self._adjuster, act_on_left=self.action().is_left(), dettwist=self._dettwist, variable_name = self.base_ring().variable_name())
    
    def change_precision(self, new_prec):
        """
            Returns a FamiliesOfOMS coefficient module with same input data as self, but with precision cap ``new_prec``
        """
        #print new_prec
        if new_prec == self._prec_cap:
            return self
        base_coeffs = self.base_ring().base_ring()
        if new_prec[0] > base_coeffs.precision_cap():
            #THERE'S NO WAY TO EXTEND PRECISION ON BASE RING!!! This is a crappy hack:
            if base_coeffs.is_field():
                base_coeffs = Qp(self.prime(), new_prec[0])
            else:
                base_coeffs = ZpCA(self.prime(), new_prec[0])
        return FamiliesOfOverconvergentDistributions(self._k, prec_cap = new_prec, base_coeffs=base_coeffs, character=self._character, adjuster=self._adjuster, act_on_left=self.action().is_left(), dettwist=self._dettwist, variable_name = self.base_ring().variable_name())
        
    @cached_method
    def approx_module(self, p_prec=None, var_prec=None):
        if p_prec is None:
            p_prec = self._prec_cap[0]
        elif p_prec < 0 or p_prec > self._prec_cap[0]:
            raise ValueError("p_prec must be between 0 and %s"%(self._prec_cap[0]))
        if var_prec is None:
            var_prec = self._prec_cap[1]
        elif var_prec < 0 or var_prec > self._prec_cap[1]:
            raise ValueError("var_prec must be between 0 and %s"%(self._prec_cap[1]))
        # Should/can we do something about the variable's precision?
        return self.base_ring() ** self.length_of_moments(p_prec)
    
    def random_element(self, prec=None):
        if prec == None:
            prec = self._prec_cap
        else:
            prec = _prec_cap_parser(prec)
        V = self.approx_module(prec[0], prec[1])
        R = self.base_ring().base_ring()
        if R.is_field():
            V = V.change_ring((self.base_ring().change_ring(R.integer_ring())))
        return self(V.random_element()) #Make this better
    
    def clear_cache(self):
        self.approx_module.clear_cache()
        self._act.clear_cache()
    
    def _an_element_(self):
        if self._prec_cap[0] > 1:
            return self([2,1])
        else:
            return self([1])
    
    def specialize(self, k, return_map=False):
        from sage.modular.pollack_stevens.coeffmod_OMS_space import OverconvergentDistributions
        D = OverconvergentDistributions(k, base=self.base_ring().base_ring(), prec_cap=self._prec_cap[0], character=self._character, adjuster=self._adjuster, act_on_left=self.action().is_left(), dettwist=self._dettwist)
        if not return_map:
            return D
        else:
            return D, None
