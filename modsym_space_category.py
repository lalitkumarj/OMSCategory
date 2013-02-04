from sage.categories.category_types import Category_module

class ModularSymbolSpaces(Category_module):
    def super_categories(self):
        from sage.categories.modules import Modules
        return [Modules(self.base_ring())]

    def __repr__(self):
        return "Category of Modular Symbol Spaces"


    class ParentMethods:
        @abstract_method
        def coefficient_module(self):
            """ """
        @abstract_method(optional = True)
        def prime(self):
           """ """
        @abstract_method
        def character(self):
            """REALLY??? CURRENT GENERIC DOES NOT CONTAIN THIS """
        @abstract_method
        def source(self):
            """ """
        @abstract_method
        def hecke(self):
            """ NOT IMPLEMENTED IN PARENT CURRENTLY. Implemented in PSmodularsymbolelement. Should this be an element method?"""
        
    class ElementMethods:

        @abstract_method
        def call(self):
            """ There is a _call in the action ... is this the same thing?"""

        
        def plus_part(self):
            r"""
            Returns the plus part of self -- i.e. self + self | [1,0,0,-1].

            Note that we haven't divided by 2.  Is this a problem?

            OUTPUT:

            - self + self | [1,0,0,-1]

            EXAMPLES::

                sage: E = EllipticCurve('11a')
                sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
                sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
                [-1/5, 3/2, -1/2]
                sage: (phi.plus_part()+phi.minus_part()) == 2 * phi
                True
            """
            #?????????????
            return self.parent.involution() + self

        def minus_part(self):
            r"""
            Returns the minus part of self -- i.e. self - self | [1,0,0,-1]

            Note that we haven't divided by 2.  Is this a problem?

            OUTPUT:

            - self - self | [1,0,0,-1]

            EXAMPLES::

                sage: E = EllipticCurve('11a')
                sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
                sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
                [-1/5, 3/2, -1/2]
                sage: (phi.plus_part()+phi.minus_part()) == phi * 2
                True
            """
            return self - self.parent.involution()

