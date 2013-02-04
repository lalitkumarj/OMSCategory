from sage.categories.all import Category

class MSCoefficientModule(Category):
    def super_categories(self):
        return []

    def __repr__(self):
        return "Category of Modular Symbol Coefficient Modules"

    class ParentMethods:
        @abstract_method
        def involution(self):
            """ """

        @abstract_method(optional = True)
        def prime(self):
            """An optional method to return p."""

        @abstract_method
        def weight(self):
            """ """

        @abtract_method
        def character(self):
            """ """

        #Not sure if we need to require a definition for action on each side.
        @abstract_method
        def action_right(self, g_right):
            """Act by g_right on the right."""

        @abstract_method
        def action_left(self, g_left):
            """Act by g_left on the left."""

        @abtract_method
        def random_element():
            """ """

        @abstract_method
        def zero_element(self, M=None):
            """ """


    class ElementMethods:

        @abstract_method
        def is_zero():
            """ """

