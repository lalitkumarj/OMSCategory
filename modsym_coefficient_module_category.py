from sage.categories.category_types import Category_module
#from category_types import Category_module (this file will be in the Categories
# folder when we're done)

class MSCoefficientModules(Category_module):
    def super_categories(self):
        from sage.categories.modules import Modules
        return [Modules(self.base_ring())]
    
    #def __repr__(self):
    #    return "Category of Modular Symbol Coefficient Modules"
    
    class ParentMethods:
        def involution(self):
            """
            This returns the involution on this space of modular symbols. This
            depends on how the action is defined.
            """
            return self.action().involution()
        
        @abstract_method
        def weight(self):
            """ """
        
        @abtract_method
        def character(self):
            """ """
        
        #Not sure if we need to require a definition for action on each side.
        @abstract_method
        def action(self):
            """Returns the action that acts on this coefficient module."""
        
        @abstract_method(optional = True)
        def prime(self):
            """
            An optional method to return p. Currently, the convention is to set
            the prime to 0 if over a ring such as QQ or CC, but this could
            change.
            """
        
    class ElementMethods:
        pass
