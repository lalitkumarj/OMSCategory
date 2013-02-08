import sage.rings.ring as ring

from sage.structure.factory import UniqueFactory
from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method

from sage.modular.pollack_stevens.sigma0 import _default_adjuster
from sage.modular.pollack_stevens.coeffmod_space import CoefficientModule_generic
from sage.modular.pollack_stevens.coeffmod_element import WeightKAction_OMS
from sage.modular.pollack_stevens.coeffmod_OMS_element import CoeffMod_OMS_element

class CoeffMod_OMS_factory(UniqueFactory):
    def create_key(self, k, p=None, prec_cap=None, base=None, \
                     character=None, adjuster=None, act_on_left=False, \
                     dettwist=None):
        k = ZZ(k)
#        if p is None:
#            try:
#                p = base.prime()
#            except AttributeError:
#                raise ValueError("You must specify a prime")
#        else:
#            p = ZZ(p)
        
        if base is None:
            if p is None:
                raise ValueError("Must specify a prime or a base ring.")
            if prec_cap is None:
                base = ZpCA(p)
            else:
                base = ZpCA(p, prec_cap)
        prec_cap = base.precision_cap()
        p = base.prime()
        if adjuster is None:
            adjuster = _default_adjuster()
        if dettwist is not None:
            dettwist = ZZ(dettwist)
            if dettwist == 0: 
                dettwist = None
        return (k, p, prec_cap, base, character, adjuster, act_on_left, dettwist)
    
    def create_object(self, version, key):
        k, p, prec_cap, base, character, adjuster, act_on_left, dettwist = key
        return CoeffMod_OMS_space(k, p=p, prec_cap=prec_cap, base=base, character=character, \
                 adjuster=adjuster, act_on_left=act_on_left, dettwist=dettwist)

OverconvergentDistributions = CoeffMod_OMS_factory('CoeffMod_OMS_Families_space')

class CoeffMod_OMS_space(CoefficientModule_generic):
    r"""
    EXAMPLES::
    #CHANGE THIS
        sage: from sage.modular.pollack_stevens.coeffmod_OMS_families_space import FamiliesOfOMS
        sage: D = FamiliesOfOMS(0, base=PowerSeriesRing(Zp(3, 10), 'w', default_prec=5)); D
        Families of overconvergent distributions on the disc 0 over Power Series Ring in w over 3-adic Ring with capped relative precision 10
        sage: D.prime()
        3
        sage: D.precision_cap()
        [10, 5]
    """
    
    def __init__(self, k, p=None, prec_cap=None, base=None, \
                 character=None, adjuster=None, act_on_left=False, \
                 dettwist=None, action_class = WeightKAction_OMS):
        #self._prec_cap = prec_cap
        
        self._p = base.prime()
        self._prec_cap = base.precision_cap()
        CoefficientModule_generic.__init__(self, k, base=base, \
                 character=character, adjuster=adjuster, act_on_left=act_on_left, \
                 dettwist=dettwist, action_class=action_class, \
                 element_class=CoeffMod_OMS_element, padic=True)
    
    def prime(self):
        return self._p
    
    def precision_cap(self):
        return self._prec_cap
    
    def change_ring(self, new_base_ring):
        if new_base_ring == self.base_ring():
            return self
        return CoeffMod_OMS_space(self._k, base=new_base_ring, character=self._character, adjuster=self._adjuster, act_on_left=self.action().is_left(), dettwist=self._dettwist, action_class = self.action().__class__)
    
    @cached_method
    def approx_module(self, prec=None):
        if prec is None:
            prec = self._prec_cap
        elif prec > self._prec_cap:
            raise ValueError("p_prec must be less than or equal to the p-adic precision cap")
        return self.base_ring()**prec
    
    def random_element(self, prec=None):
        if prec == None:
            prec = self._prec_cap
        return self(self.approx_module(prec).random_element())
    
    def clear_cache(self):
        self.approx_module.clear_cache()
        self._act.clear_cache()
    
    def _an_element_(self):
        if self._prec_cap > 1:
            return self([2,1])
        else:
            return self([1])
    
    def lift(self, var_prec=None, variable_name='w'):
        from sage.modular.pollack_stevens.coeffmod_OMS_families_space import FamiliesOfOMS
        if var_prec is None:
            raise ValueError("Must specify var_prec: the precision on the variable.")
        return FamiliesOfOMS(self._k, prec_cap=[self._prec_cap,var_prec], base_coeffs=self.base_ring(), character=self._character, adjuster=self._adjuster, act_on_left=self.action().is_left(), dettwist=self._dettwist, variable_name=variable_name)