from sage.modules.module import Module
from sage.modular.pollack_stevens.coeffmod_element import *
from sage.structure.parent import Parent
from sage.categories.modsym_coefficient_module_category import MSCoefficientModules

class CoefficientModule_generic(Module):
    def __init__(self, k, base=None, character=None, \
                 adjuster=None, act_on_left=False, dettwist=None, \
                 action_class = WeightKAction_generic, \
                 element_class = CoefficientModuleElement_generic, padic=False):
        self.Element = element_class
        Parent.__init__(self, base, category=MSCoefficientModules(base))
        self._k = k
        self._character = character
        self._adjuster=adjuster
        self._dettwist=dettwist
        #print action_class
        self._act = action_class(self, character, adjuster, act_on_left, dettwist, padic=padic)
        self._populate_coercion_lists_(action_list=[self._act])
    
    # Category framework
    def weight(self):
        """ """
        return self._k
    
    def character(self):
        """ """
        return self._character
    
    def action(self):
        """Returns the action that acts on this coefficient module."""
        return self._act
    # End category framework
    
    def change_ring(self, new_base_ring):
        pass #TODO
    
    def base_extend(self, new_base_ring):
        pass #TODO