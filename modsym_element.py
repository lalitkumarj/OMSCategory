class PSModularSymbolsElement(ModuleElement):
    def __init__(self, map_data, parent, construct = False):
        ModuleElement.__init__(self, parent)
        if construct:
            self._map = map_data
        else:
            self._map = ManinMap(parent._coefficients, parent._source, map_data)
    
    def _repr_(self):
        pass
    
    