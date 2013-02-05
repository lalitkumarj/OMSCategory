class CoefficientModule_generic(Module):
    def __init__(self, k, base=None, character=None, \
                 adjuster=None, act_on_left=False, dettwist=None):
        Parent.__init__(self, base, category=MSCoefficientModules(base))
    
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