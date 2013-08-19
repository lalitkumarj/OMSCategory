#from sage.categories.category_types import Category_module
# (this file will be in the Categories folder when we're done)
from sage.categories.category_types import Category_module
from sage.misc.abstract_method import abstract_method, abstract_methods_of_class

class ModularSymbolSpaces(Category_module):
    r"""
    The category of modular symbols (over a given base ring).
    
    Objects in this category are spaces of modular symbols (a la Mazur), i.e.
    spaces of the form `\Hom_\Gamma(\text{Div}^0\mathbf{P}^1(\QQ),D)`. The codomain
    `D` is called a coefficient module and is a object in the category described in
    :class:`MSCoefficientModules`. It comes with an action by some monoid of matrices
    which endows the modular symbol space with a Hecke action. Related to this action,
    there is generally a weight and a character.
    
    EXAMPLES::
    
        sage: from sage.categories.modsym_space_category import ModularSymbolSpaces
        sage: ModularSymbolSpaces(ZZ)
        Category of modular symbol spaces over Integer Ring
    """
    
    def super_categories(self):
        """
        EXAMPLES::
        
            sage: from sage.categories.modsym_space_category import ModularSymbolSpaces
            sage: ModularSymbolSpaces(QQ).super_categories()
            [Category of vector spaces over Rational Field]
            sage: ModularSymbolSpaces(ZZ).super_categories()
            [Category of modules over Integer Ring]
        """
        from sage.categories.modules import Modules
        return [Modules(self.base_ring())]
    
    #def required_methods(self):
    #    return { "parent"  : abstract_methods_of_class(self.parent_class),
    #             "element" : abstract_methods_of_class(self.element_class) }
    
    # This might be done automagically!
    #def __repr__(self):
    #    return "Category of Modular Symbol Spaces"
    
    class ParentMethods:
        @abstract_method
        def source(self):
            """
            A space of modular symbols is a certain set of homormophisms. This
            returns the domain.
            """
        
        @abstract_method
        def coefficient_module(self):
            r"""
            A space of modular symbols is a certain set of homormophisms. This
            returns the codomain, which should be in the category
            :class:`MSCoefficientModules` over the same base ring as self.
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
            
        @abstract_method
        def hecke(self, ell):
            """
            Returns a Hecke endomorphism of self.
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

        @abstract_method
        def __call__(self, a):
            """
            A modular symbol is a homomorphism, so one should be able to evaluate
            it. This is accomplished by this function.
            """
        
        #def plus_part(self):
        #    r"""
        #    Returns the plus part of self -- i.e. self + self | [1,0,0,-1].
        #
        #    Note that we haven't divided by 2.  Is this a problem?
        #
        #    OUTPUT:
        #
        #    - self + self | [1,0,0,-1]
        #
        #    EXAMPLES::
        #
        #        sage: E = EllipticCurve('11a')
        #        sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
        #        sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
        #        [-1/5, 3/2, -1/2]
        #        sage: (phi.plus_part()+phi.minus_part()) == 2 * phi
        #        True
        #    """
        #    return self + self.action().involution()
        #
        #def minus_part(self):
        #    r"""
        #    Returns the minus part of self -- i.e. self - self | [1,0,0,-1]
        #
        #    Note that we haven't divided by 2.  Is this a problem?
        #
        #    OUTPUT:
        #
        #    - self - self | [1,0,0,-1]
        #
        #    EXAMPLES::
        #
        #        sage: E = EllipticCurve('11a')
        #        sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
        #        sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
        #        [-1/5, 3/2, -1/2]
        #        sage: (phi.plus_part()+phi.minus_part()) == phi * 2
        #        True
        #    """
        #    return self - self.action().involution()
