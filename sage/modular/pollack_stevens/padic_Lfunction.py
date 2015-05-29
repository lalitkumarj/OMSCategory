from sage.structure.sage_object import SageObject
from sage.rings.power_series_ring import PowerSeriesRing

class padic_Lfunction(SageObject):
    pass

class padic_Lfunction_one_variable(padic_Lfunction):
    def __init__(self, Phi):
        pass

class padic_Lfunction_two_variable(padic_Lfunction):
    def __init__(self, Phis, var='T', prec=None):
        #TODO: prec: Default it to None would make us compute it.
        self._Phis = Phis    #should we create a copy of Phis, in case Phis changes?
        self._coefficient_ring = Phis.base_ring()
        self._base_ring = PowerSeriesRing(self._coefficient_ring, var)    #What precision?
    
    def base_ring(self):
        return self._base_ring
    
    def coefficient_ring(self):
        return self._coefficient_ring
    
    def variables(self):
        #returns (T, w)
        return (self._base_ring.gens()[0], self._coefficient_ring.gens()[0])
    
    def _basic_integral(self, a, j):
        #if Phis is fixed for this p-adic L-function, we should make this method cached
        pass
    
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
        pass
