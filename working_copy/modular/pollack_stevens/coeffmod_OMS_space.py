import sage.rings.ring as ring

from sage.misc.prandom import randint
from sage.structure.factory import UniqueFactory
from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method
from sage.rings.padics.factory import ZpCA,Qp

from sage.modular.pollack_stevens.sigma0 import _default_adjuster
from sage.modular.pollack_stevens.coeffmod_space import CoefficientModule_generic
from sage.modular.pollack_stevens.coeffmod_element import WeightKAction_OMS
from sage.modular.pollack_stevens.coeffmod_OMS_element import CoeffMod_OMS_element

class CoeffMod_OMS_factory(UniqueFactory):
    """
    EXAMPLES::
    
        sage: from sage.modular.pollack_stevens.coeffmod_OMS_space import OverconvergentDistributions
        sage: D = OverconvergentDistributions(20, 3, 10); D    # indirect doctest
        Space of 3-adic distributions with k=20 action and precision cap 10
        sage: TestSuite(OverconvergentDistributions).run()
    """
    def create_key(self, k, p=None, prec_cap=None, base=None, \
                     character=None, adjuster=None, act_on_left=False, \
                     dettwist=None):
        k = ZZ(k)
        if base is None:
            if p is None:
                raise ValueError("Must specify a prime or a base ring.")
            if prec_cap is None:
                base = ZpCA(p)
            else:
                base = ZpCA(p, prec_cap)
        if prec_cap is None:
            prec_cap = base.precision_cap()
        elif prec_cap > base.precision_cap():
            raise ValueError("Insufficient precision in base ring (%s < %s)."%(base.precision_cap(), prec_cap))
        if p is None:
            p = base.prime()
        elif p != base.prime():
            raise ValueError("Prime p(=%s) must equal the prime of the base ring(=%s)"%(p, base.prime()))
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

OverconvergentDistributions = CoeffMod_OMS_factory('OverconvergentDistributions')

class CoeffMod_OMS_space(CoefficientModule_generic):
    r"""
    EXAMPLES::
    
        sage: from sage.modular.pollack_stevens.coeffmod_OMS_space import OverconvergentDistributions
        sage: D = OverconvergentDistributions(0, 5, 10)
    
    TEST::
    
        sage: TestSuite(D).run()
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
    
    #def __reduce__(self):
    #    return OverconvergentDistributions.reduce_data(self)
    
    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.modular.pollack_stevens.coeffmod_OMS_space import OverconvergentDistributions
            sage: OverconvergentDistributions(0, 5, 10)._repr_()
            'Space of 5-adic distributions with k=0 action and precision cap 10'
            sage: OverconvergentDistributions(0, 5, 10)
            Space of 5-adic distributions with k=0 action and precision cap 10

        Examples with twists::

            sage: OverconvergentDistributions(0,3,4)                               
            Space of 3-adic distributions with k=0 action and precision cap 4
            sage: OverconvergentDistributions(0,3,4,dettwist=-1)
            Space of 3-adic distributions with k=0 action and precision cap 4 twistted by det^-1
            sage: OverconvergentDistributions(0,3,4,character=DirichletGroup(3).0)
            Space of 3-adic distributions with k=0 action and precision cap 4 twistted by (Dirichlet character modulo 3 of conductor 3 mapping 2 |--> -1)
            sage: OverconvergentDistributions(0,3,4,character=DirichletGroup(3).0,dettwist=-1)
            Space of 3-adic distributions with k=0 action and precision cap 4 twistted by det^-1 * (Dirichlet character modulo 3 of conductor 3 mapping 2 |--> -1)
        """
        s = "Space of %s-adic distributions with k=%s action and precision cap %s"%(self._p, self._k, self._prec_cap)
        twiststuff = []
        if self._dettwist is not None:
           twiststuff.append("det^%s" % self._dettwist)
        if self._character is not None:
            twiststuff.append("(%s)" % self._character)
        if twiststuff:
            s += " twistted by " + " * ".join(twiststuff)
        return s
    
    def prime(self):
        return self._p
    
    def precision_cap(self):
        return self._prec_cap
    
    def change_ring(self, new_base_ring):
        if new_base_ring == self.base_ring():
            return self
        return OverconvergentDistributions(self._k, prec_cap = self._prec_cap, base=new_base_ring, character=self._character, adjuster=self._adjuster, act_on_left=self.action().is_left(), dettwist=self._dettwist)
    
    def change_precision(self, new_prec):
        """
            Returns an OMS coefficient module with same input data as self, but with precision cap ``new_prec``
        """
        if new_prec == self._prec_cap:
            return self
        base = self.base_ring()
        if new_prec > base.precision_cap():
            #THERE'S NO WAY TO EXTEND PRECISION ON BASE RING!!! This is a crappy hack:
            if self.base_ring().is_field():
                base = Qp(self.prime(), new_prec)
            else:
                base = ZpCA(self.prime(), new_prec)
        return OverconvergentDistributions(self._k, prec_cap = new_prec, base=base, character=self._character, adjuster=self._adjuster, act_on_left=self.action().is_left(), dettwist=self._dettwist)
    
    @cached_method
    def approx_module(self, prec=None):
        #print "prec_in, self._prec_cap", prec, self._prec_cap
        if prec is None:
            prec = self._prec_cap
        elif prec > self._prec_cap:
            raise ValueError("p_prec must be less than or equal to the p-adic precision cap")
        return self.base_ring()**prec
    
    def random_element(self, prec=None, total_measure_zero=False):
        #RH: copied from RP's change
        if prec == None:
            prec = self._prec_cap
        if prec == 0:
            V = self.approx_module(0)
            return CoeffMod_OMS_element(V([]), self, ordp=randint(1, self._prec_cap), check=False)
        R = self.base_ring().integer_ring()
        return self((R**prec).random_element())
    
    def clear_cache(self):
        self.approx_module.clear_cache()
        self._act.clear_cache()
    
    def _an_element_(self):
        if self._prec_cap > 1:
            return self([2,1])
        else:
            return self([1])
    
    def lift(self, var_prec=None, variable_name='w'):
        #RH: should be good
        from sage.modular.pollack_stevens.coeffmod_OMS_families_space import FamiliesOfOverconvergentDistributions
        if var_prec is None:
            #raise ValueError("Must specify var_prec: the precision on the variable.")
            var_prec = self._prec_cap
        return FamiliesOfOverconvergentDistributions(self._k, prec_cap=[self._prec_cap,var_prec], base_coeffs=self.base_ring(), character=self._character, adjuster=self._adjuster, act_on_left=self.action().is_left(), dettwist=self._dettwist, variable_name=variable_name)
