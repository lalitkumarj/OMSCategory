from sage.modules.module import Module
from sage.structure.parent import Parent
from sage.rings.padics.factory import ZpCA
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method
from sage.categories.action import PrecomposedAction
from sage.structure.coerce_actions import LeftModuleAction, RightModuleAction
from sage.matrix.all import MatrixSpace
from sage.rings.fast_arith import prime_range
from sage.modular.pollack_stevens.dist import get_dist_classes, Dist_long, iScale
from sage.structure.factory import UniqueFactory
from sage.structure.unique_representation import UniqueRepresentation
import operator
import sage.rings.ring as ring


class Distributions_factory(UniqueFactory):
 


Distributions_Families = Distributions_Families_factory('Distributions')

class Families_class(Distributions_abstract):

