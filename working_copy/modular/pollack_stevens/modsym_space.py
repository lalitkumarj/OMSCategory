from sage.modules.module import Module
from sage.structure.parent import Parent
from sage.categories.modsym_space_category import ModularSymbolSpaces

from sage.modular.pollack_stevens.fund_domain import ManinRelations
from sage.modular.pollack_stevens.modsym_element import ModularSymbolElement_generic

class ModularSymbolSpace_generic(Module):
    def __init__(self, group, coefficients, sign=0, element_class = ModularSymbolElement_generic):
        R = coefficients.base_ring()
        Parent.__init__(self, base = R, category = ModularSymbolSpaces(R))
        if sign not in [0,-1,1]:
            # sign must be be 0, -1 or 1
            raise ValueError("sign must be 0, -1, or 1")
        self._group = group
        self._coefficients = coefficients
        self.Element = element_class
        self._sign = sign
        self._source = ManinRelations(group.level())
        #other stuff
    
    # Category framework
    def coefficient_module(self):
        r"""
        A space of modular symbols is a certain set of homormophisms. This
        returns the target space, which should be in the category
        :class:`MSCoefficientModules` over the same base ring as self.
        """
        return self._coefficients
    
    def character(self):
        #CURRENT GENERIC DOES NOT CONTAIN THIS
        """
        
        """
        return self.coefficient_module()._character
    
    def source(self):
        """
        
        """
        return self._source
    
    def hecke(self, ell, algorithm = "prep"):
        r"""
        Returns the Hecke operator `T_\ell` as endomorphism of self.
        """
        #Can this be written at this level?
        pass
    
    #def prime(self):
    #    """
    #    
    #    """
    #    return self.coefficient_module()._p
    # End category framework
    
    #def _repr_(self):
    #    pass
    