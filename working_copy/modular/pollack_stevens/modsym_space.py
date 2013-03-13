from sage.modules.module import Module
from sage.structure.parent import Parent
from sage.categories.modsym_space_category import ModularSymbolSpaces

from sage.modular.pollack_stevens.fund_domain import ManinRelations
from sage.modular.pollack_stevens.modsym_element import ModularSymbolElement_generic, ModSymAction
from sage.modular.pollack_stevens.sigma0 import Sigma0

class ModularSymbolSpace_generic(Module):
    def __init__(self, group, coefficients, sign=0, element_class = ModularSymbolElement_generic):
        R = coefficients.base_ring()
        Parent.__init__(self, base = R, category = ModularSymbolSpaces(R))
        if sign not in (0,-1,1):
            # sign must be be 0, -1 or 1
            raise ValueError("sign must be 0, -1, or 1")
        self._group = group
        self._coefficients = coefficients
        self.Element = element_class
        self._sign = sign
        self._source = ManinRelations(group.level())
        try:
            action = ModSymAction(Sigma0(self.prime()), self)
        except AttributeError:
            action = ModSymAction(Sigma0(1), self)
        self._populate_coercion_lists_(action_list=[action])
    
    # Category framework
    def coefficient_module(self):
        r"""
        A space of modular symbols is a certain set of homormophisms. This
        returns the target space, which should be in the category
        :class:`MSCoefficientModules` over the same base ring as self.
        """
        return self._coefficients
    
    def weight(self):
        return self.coefficient_module().weight()
    
    def character(self):
        #CURRENT GENERIC DOES NOT CONTAIN THIS
        """
        
        """
        return self.coefficient_module().character()
    
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
    
    def _element_constructor(self, data):
        if isinstance(data, ModularSymbolSpace_generic):
            data = data._map
        elif isinstance(data, ManinMap):
            pass
        else:
            # a dict, or a single distribution specifying a constant symbol, etc
            data = ManinMap(self._coefficients, self._source, data)

        if data._codomain != self._coefficients:
            data = data.extend_codomain(self._coefficients)
 
        return self.element_class(data, self, construct=True)
    
    def _repr_(self):
        return "Generic space of modular symbols for %s with sign %s and values in %s"%(self.group(), self.sign(), self.coefficient_module())
    
    def group(self):
        return self._group
    
    def sign(self):
        #ADD THIS TO CATEGORY SETUP?
        return self._sign
    
    def ngens(self):
        return len(self.source().indices())
    
    def number_of_coset_reps(self):
        return len(self.source().reps())
    
    def level(self):
        return self.source().level()
    
    #def prime(self):
    #    """
    #    
    #    """
    #    return self.coefficient_module()._p
    # End category framework
    
    #def _repr_(self):
    #    pass
    