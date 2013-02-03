class AA(Category):
    def super_categories(self):
        return []

    def __repr__(self):
        return "AAAAAAAAAAAAAAA"
    
    class ParentMethods:
        def helloz(self):
            return 21

    class ElementMethods:
        def blah(self):
            return 22

class EL:
    def __init__(self):
	print("hi")


class PQ(Parent):
    def __init__(self):
        Parent.__init__(self,facade=None,category=AA())
    
    def hello(self):
        return 34
    Element = EL

