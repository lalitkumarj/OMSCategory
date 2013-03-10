import operator

from sage.categories.action import Action
from sage.structure.element import ModuleElement

from sage.modular.pollack_stevens.manin_map import ManinMap

class ModSymAction(Action):
    def __init__(self, actor, MSspace):
        Action.__init__(self, actor, MSspace, False, operator.mul)

    def _call_(self, symb, g):
        return symb.__class__(symb._map * g, symb.parent(), construct=True)

class ModularSymbolElement_generic(ModuleElement):
    def __init__(self, map_data, parent, construct = False):
        ModuleElement.__init__(self, parent)
        if construct:
            self._map = map_data
        else:
            self._map = ManinMap(parent._coefficients, parent._source, map_data)
    
    #def _repr_(self):
    #    pass
    
    #def __call__(self):
    #    pass
    
    