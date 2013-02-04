class CoefficientModule_generic(Module):
    def __init__(self, k, base=None, character=None, \
                 adjuster=None, act_on_left=False, dettwist=None):
        Parent.__init__(self, base, category=Modules(base))
    
    def