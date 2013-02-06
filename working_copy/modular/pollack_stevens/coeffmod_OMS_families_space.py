import sage.rings.ring as ring

from sage.structure.factory import UniqueFactory
from sage.rings.integer_ring import ZZ
from sage.rings.power_series_ring import PowerSeriesRing
from sage.misc.cachefunc import cached_method

from sage.modular.pollack_stevens.sigma0 import _default_adjuster
from sage.modular.pollack_stevens.coeffmod_space import CoefficientModule_generic
from sage.modular.pollack_stevens.coeffmod_element import WeightKAction_OMS_fam

class CoeffMod_OMS_Families_factory(UniqueFactory):
    def create_key(self, k, p=None, prec_cap=None, base=None, base_coeffs=None, \
                     character=None, adjuster=None, act_on_left=False, \
                     dettwist=None, variable_name = 'w'):
        k = ZZ(k)
        if base is None:
            if base_coeffs is None:
                if p is None:
                    raise ValueError("Must specify a prime, a base ring, or coefficients.")
                if prec_cap is None:
                    raise ValueError("Must specify a precision cap or a base ring.")
                if not isinstance(prec_cap, (list, tuple)):
                    prec_cap = [prec_cap, prec_cap]
                base_coeffs = Zp(p, prec = prec_cap[0])
            elif not isinstance(base_coeffs, ring.Ring):
                raise TypeError("base_coeffs must be a ring if base is None")
            base = PowerSeriesRing(base_coeffs, name=variable_name, default_prec=prec_cap[1])
        base_coeffs = None
        
        #if prec_cap is None:
        #    prec_cap = [base.base_ring().precision_cap, base.default_prec()]
        #else:
        #    prec_cap = list(prec_cap)
        prec_cap = [base.base_ring().precision_cap, base.default_prec()]
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
    
        sage: D = FamiliesOfOMS(0, base=PowerSeriesRing(Zp(3, 10), default_prec=5)); D
        Families of overconvergent distributions of weight 0 over ??
        sage: D.prime()
        5
        sage: D.precision_cap()
        [10, 5]
    """
    
    def __init__(self, k, p=None, prec_cap=[20, 10], base=None, base_coeffs=None, \
                 character=None, adjuster=None, act_on_left=False, \
                 dettwist=None, action_class = WeightKAction_OMS_fam, \
                 variable_name = 'w'):
        #TODO: deal input of prec_cap
        self._prec_cap = prec_cap
        
        if base is None:
            if base_coeffs is None:
                base_coeffs = Zp(p, prec = prec_cap[0])
            elif not isinstance(base_coeffs, ring.Ring):
                raise TypeError("base_coeffs must be a ring if base is None")
            base = PowerSeriesRing(base_coeffs, name=variable_name)
        #elif not isinstance(base, ring.Ring):
        #    raise TypeError("base must be a ring")
        self._p = base.base_ring().prime()
        CoefficientModule_generic.__init__(self, k, base=base, \
                 character=character, adjuster=adjuster, act_on_left=act_on_left, \
                 dettwist=dettwist, action_class=action_class, padic=True)
    
    def _repr_(self):
        s = "Families of overconvergent distributions on the disc %s " % (self._k)
        if self._character is not None:
            s += "twisted by %s " % (self._character)
        if self._dettwist is not None:
            s += "and det^%s " % (self._dettwist)
        s += "over %" % (self.base_ring())
        return s
    
    def prime(self):
        return self._p
    
    def precision_cap(self):
        return self._prec_cap
    
    @cached_method
    def approx_module(self, prec):
        pass #TODO
    
    def random_element(self, prec):
        pass #TODO
    
    def clear_cache(self):
        self.approx_module.clear_cache()
        self._act.clear_cache()
    
    def _an_element(self):
        pass #TODO
    
    def specialize(self, k):
        pass #TODO
