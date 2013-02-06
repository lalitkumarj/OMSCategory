import sage.rings.ring as ring

from sage.structure.factory import UniqueFactory
from sage.rings.integer_ring import ZZ
from sage.rings.power_series_ring import PowerSeriesRing
from sage.misc.cachefunc import cached_method

from sage.modular.pollack_stevens.sigma0 import _default_adjuster
from sage.modular.pollack_stevens.coeffmod_space import CoefficientModule_generic
from sage.modular.pollack_stevens.coeffmod_element import WeightKAction_OMS_fam
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
        [5, 5]
        sage: _prec_cap_parser([20, 10])
        [20, 10]
        sage: _prec_cap_parser([7])
        [7, 7]
        sage: _prec_cap_parser(None)
        ValueError: The case prec_cap is None should be dealt with prior to calling this function.
        sage: _prec_cap_parser([5, 0])
        ValueError: Precisions should be at least 1.
    """
    if prec_cap is None:
        raise ValueError("The case prec_cap is None should be dealt with prior to calling this function.")
    if not isinstance(prec_cap, (list, tuple)):
        prec = ZZ(prec_cap)
        if prec < 1:
            ValueError("prec_cap should be at least 1.")
        return [prec, prec]
    elif len(prec_cap) == 1:
        prec = ZZ(prec_cap[0])
        if prec < 1:
            ValueError("prec_cap should be at least 1.")
        return [prec, prec]
    elif len(prec_cap) > 2:
        raise ValueError("prec_cap should not have length > 2.")
    else:
        p_prec = ZZ(prec_cap[0])
        var_prec = ZZ(prec_cap[1])
        if p_prec < 1 or var_prec < 1:
            raise ValueError("Precisions should be at least 1.")
        return [p_prec, var_prec]

class CoeffMod_OMS_Families_factory(UniqueFactory):
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
            base = PowerSeriesRing(base_coeffs, name=variable_name, default_prec=prec_cap[1])
        base_coeffs = None
        p = base.base_ring().prime()
        k = ZZ(k % (p-1))
        
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
        return (k, p, base, character, adjuster, act_on_left, dettwist)
    
    def create_object(self, version, key):
        k, p, base, character, adjuster, act_on_left, dettwist = key
        return CoeffMod_OMS_Families_space(k, p=p, base=base, character=character, \
                 adjuster=adjuster, act_on_left=act_on_left, dettwist=dettwist)

FamiliesOfOMS = CoeffMod_OMS_Families_factory('CoeffMod_OMS_Families_space')

class CoeffMod_OMS_Families_space(CoefficientModule_generic):
    r"""
    EXAMPLES::
    
        sage: from sage.modular.pollack_stevens.coeffmod_OMS_families_space import FamiliesOfOMS
        sage: D = FamiliesOfOMS(0, base=PowerSeriesRing(Zp(3, 10), 'w', default_prec=5)); D
        Families of overconvergent distributions on the disc 0 over Power Series Ring in w over 3-adic Ring with capped relative precision 10
        sage: D.prime()
        3
        sage: D.precision_cap()
        [10, 5]
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
        self._prec_cap = [base.base_ring().precision_cap(), base.default_prec()]
        k = k % (self._p - 1)
        CoefficientModule_generic.__init__(self, k, base=base, \
                 character=character, adjuster=adjuster, act_on_left=act_on_left, \
                 dettwist=dettwist, action_class=action_class, \
                 element_class=CoeffMod_OMS_Families_element, padic=True)
    
    def _repr_(self):
        s = "Families of overconvergent distributions on the disc %s " % (self._k)
        if self._character is not None:
            s += "twisted by %s " % (self._character)
        if self._dettwist is not None:
            s += "and det^%s " % (self._dettwist)
        s += "over %s" % (self.base_ring())
        return s
    
    def prime(self):
        return self._p
    
    def precision_cap(self):
        return self._prec_cap
    
    def change_ring(self, new_base_ring):
        if new_base_ring == self.base_ring():
            return self
        return CoeffMod_OMS_Families_space(self._k, prec_cap=self._prec_cap, base=new_base_ring, character=self._character, adjuster=self._adjuster, act_on_left=self.action().is_left(), dettwist=self._dettwist, action_class = self.action().__class__, variable_name = self.base_ring().variable_name())
        
    @cached_method
    def approx_module(self, p_prec=None, var_prec=None):
        if p_prec is None:
            p_prec = self._prec_cap[0]
        else:
            if var_prec is None:
                var_prec = p_prec
            p_prec, var_prec = _prec_cap_parser([p_prec, var_prec])
        if p_prec > self._prec_cap[0]:
            raise ValueError("p_prec must be less than or equal to the p-adic precision cap")
        #Should/can we do something about the variable's precision
        return self.base_ring()**p_prec
    
    def random_element(self, prec):
        pass #TODO
    
    def clear_cache(self):
        self.approx_module.clear_cache()
        self._act.clear_cache()
    
    def _an_element(self):
        pass #TODO
    
    def specialize(self, k):
        pass #TODO