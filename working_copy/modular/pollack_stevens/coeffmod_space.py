from sage.modules.module import Module
from sage.modular.pollack_stevens.coeffmod_element import *

class CoefficientModule_generic(Module):
    def __init__(self, k, base=None, character=None, \
                 adjuster=None, act_on_left=False, dettwist=None, \
                 action_class = WeightKAction_generic, padic=False):
        self.Element = CoefficientModuleElement_generic
        Parent.__init__(self, base, category=MSCoefficientModules(base))
        self._k = k
        self._character = character
        self._adjuster=adjuster
        self._dettwist=dettwist
        self._act = action_class(self, character, adjuster, act_on_left, dettwist, padic=padic)
    
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