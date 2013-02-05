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
                 adjuster=adjust, act_on_left=act_on_left, dettwist=dettwist)

FamiliesOfOMS = CoeffMod_OMS_Families_factory('CoeffMod_OMS_Families_space')

class CoeffMod_OMS_Families_space(CoefficientModule_generic):
    def __init__(self, k, p=None, prec_cap=[20, 10], base=None, base_coeffs=None, \
                 character=None, adjuster=None, act_on_left=False, \
                 dettwist=None, action_class = WeightKAction_fam, \
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
