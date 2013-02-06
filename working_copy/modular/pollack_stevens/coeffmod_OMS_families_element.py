# This should be cythoned once it's done.

from sage.modular.pollack_stevens.coeffmod_element import CoefficientModuleElement_generic
from sage.rings.integer_ring import ZZ

class CoeffMod_OMS_Families_element(CoefficientModuleElement_generic):
    # Implementation currently ignores ordp
    def __init__(self, moments, parent, ordp=0, check=True):
        CoefficientModuleElement_generic.__init__(self, parent)
        #TODO: finish this
        if check:
            if isinstance(moments, CoeffMod_OMS_Families_element):
                moment = moments._moments.change_ring(parent.base_ring())
            elif hasattr(moments, '__len__'):
                M = len(moments)
                moments = parent.approx_module(M)(moments)
            elif moments == 0:
                moments = parent.approx_module(parent._prec_cap[0], parent._prec_cap[1])(moments)
            else:
                moments = parent.approx_module(1, parent.precision_cap()[1])([moments])
        self._moments = moments
        #if var_prec is None:
        #    self._var_prec = parent._prec_cap[1]
        #else:
        #    self._var_prec = var_prec
        p = parent._p
        self._cp = (p-2) / (p-1)
    
    #def _relprec(self):
    #    return len(self._moments)
    def _repr_(self):
        
        return repr(self._moments)
    
    def precision_relative(self):
        return [ZZ(len(self._moments)), self._var_prec]
    
    def precision_absolute(self):
        return [ZZ(len(self._moments)) + self.ordp, self._var_prec]
    
    def normalize(self):
        #customized
        p = self.parent()._p
        moments = self._moments
        M = len(moments)
        R = self.parent().base_ring()
        #new_moments = [m.add_bigoh(self._var_prec) for m in self._moments]
        new_moments = []
        for i in range(M):
            cutoff = ceil((M-i) * self._cp)
            f = moments[i].list()
            f = [f[j].add_bigoh(cutoff) for j in range(self._var_prec)]
            moments[i] = R(f, self._var_prec)
        return self
    
    def reduce_precision(self, new_prec):
        pass    #TODO
    
    def solve_diff_eqn(self):
        pass    #TODO important!