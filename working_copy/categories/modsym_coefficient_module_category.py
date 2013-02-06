#from sage.categories.category_types import Category_module
#(this file will be in the Categories folder when we're done)
from sage.categories.category_types import Category_module
from sage.misc.abstract_method import abstract_method, abstract_methods_of_class

class MSCoefficientModules(Category_module):
    r"""
    The category of modular symbol coefficient modules (over a given base ring).
    
    Objects in this category are modules `D` that appear as codomains of spaces
    of modular symbols (a la Mazur), i.e. spaces of the form
    `\Hom_\Gamma(\text{Div}^0\mathbf{P}^1(\QQ),D)`. As such, they have an action
    by some monoid of matrices. Related to this action, there is generally a
    weight and a character.
    
    EXAMPLES::
    
        sage: from sage.categories.modsym_coefficient_module_category import MSCoefficientModules
        sage: MSCoefficientModules(ZZ)
        Category of modular symbol coefficient modules over Integer Ring
    """
    def __init__(self, R):
        """
        TESTS::
        
            sage: from sage.categories.modsym_coefficient_module_category import MSCoefficientModules
            sage: D = MSCoefficientModules(QQ); TestSuite(D).run()
        """
        Category_module.__init__(self, R, 'modular symbol coefficient modules')
    
    def super_categories(self):
        """
        EXAMPLES::
        
            sage: from sage.categories.modsym_coefficient_module_category import MSCoefficientModules
            sage: MSCoefficientModules(QQ).super_categories()
            [Category of vector spaces over Rational Field]
            sage: MSCoefficientModules(ZZ).super_categories()
            [Category of modules over Integer Ring]
        """
        from sage.categories.modules import Modules
        return [Modules(self.base_ring())]
    
    #def __repr__(self):
    #    return "Category of Modular Symbol Coefficient Modules"
    
    class ParentMethods:
        #def involution(self):
        #    """
        #    This returns the involution on this space of modular symbols. This
        #    depends on how the action is defined.
        #    """
        #    return self.action().involution()
        
        @abstract_method
        def action(self):
            """
            Returns the action (generally as an :class:`Action` object) that acts
            on this coefficient module.
            """
            
        @abstract_method
        def weight(self):
            r"""
            Returns the weight associated to the action. This could be an integer
            describing the weight `k` action of `SL(2,\ZZ)` on `\text{Sym}^k(\ZZ^2)`
            or some more general notion.
            """
        
        @abstract_method
        def character(self):
            """
            Returns the character assocaited to the action, such as in cases
            where there is a Nebentypus.
            """
        
        @abstract_method(optional = True)
        def prime(self):
            r"""
            An optional method to return the prime of the residue field of the
            base ring, if that makes sense. Currently, the convention is to set
            the prime to 0 if over a ring such as `\QQ` or `\CC`, but this could
            change.
            """
        
    class ElementMethods:
        pass
