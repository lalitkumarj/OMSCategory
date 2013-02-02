import MSCoefficientModule
#from sage.modules.module import Module
from sage.structure.parent import Parent
from sage.rings.padics.factory import ZpCA
#from sage.rings.padics.padic_generic import pAdicGeneric
#from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
#from sage.misc.cachefunc import cached_method
#from sage.categories.action import PrecomposedAction
#from sage.structure.coerce_actions import LeftModuleAction, RightModuleAction
#from sage.matrix.all import MatrixSpace
#from sage.rings.fast_arith import prime_range
from sage.modular.pollack_stevens.dist import get_dist_classes, Dist_long, iScale
#import sage.rings.ring as ring


class DistributionsSpace(Parent):
    r"""
    INPUT:
        
    - `k` -- nonnegative integer
    - `p` -- prime number or None
    - ``prec_cap`` -- positive integer or None
    - ``base`` -- ring or None
    - ``symk`` -- bool or None
    - ``character`` -- a dirichlet character or None
    - ``tuplegen`` -- None or callable that turns 2x2 matrices into a 4-tuple
    - ``act_on_left`` -- bool (default: False)
        
    EXAMPLES:: #must be updated
        
    sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
    sage: Distributions(20, 3, 10)              # indirect doctest
    Space of 3-adic distributions with k=20 action and precision cap 10
    """
    
    def __init__(self, k, p=None, prec_cap=None, base=None, symk=None, character=None, act_on_left=False):
        
        """
        See ``DistributionsSpace`` for full documentation.
        """
        
        Parent.__init__(self,category=MSCoefficientModule)
        Element = DistributionElementPy #do we want elements to be DistributionElementPy or DistributionElementBase

        k = ZZ(k)
        if p is None:
            try:
                p = base.prime()
            except AttributeError:
                raise ValueError("You must specify a prime")
        else:
            p = ZZ(p)
        if base is None:
            if prec_cap is None:
                base = ZpCA(p)
            else:
                base = ZpCA(p, prec_cap)
        if prec_cap is None:
            try:
                prec_cap = base.precision_cap()
            except AttributeError:
                raise ValueError("You must specify a base or precision cap")
        return (k, p, prec_cap, base, character, tuplegen, act_on_left, symk)

    def _repr_(self):
        """
        EXAMPLES::

        sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        sage: Distributions(0, 5, 10)._repr_()
        'Space of 5-adic distributions with k=0 action and precision cap 10'
        sage: Distributions(0, 5, 10)
        Space of 5-adic distributions with k=0 action and precision cap 10
        sage: Symk(0)
        Sym^0 Q^2
        """
        # TODO: maybe account for character, etc.
        return "Space of %s-adic distributions with k=%s action and precision cap %s"%(self._p, self._k, self._prec_cap)

    def change_ring(self, new_base_ring):
        """
        Return space of distributions like this one, but with the base ring changed.

        EXAMPLES::

        sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        sage: D = Distributions(0, 7, 4); D
        Space of 7-adic distributions with k=0 action and precision cap 4
        sage: D.base_ring()
        7-adic Ring with capped absolute precision 4
        sage: D2 = D.change_ring(QpCR(7)); D2
        Space of 7-adic distributions with k=0 action and precision cap 4
        sage: D2.base_ring()
        7-adic Field with capped relative precision 20
        """
        return DistributionsSpace(k=self._k, p=self._p, prec_cap=self._prec_cap, base=new_base_ring, character=self._character, tuplegen=self._act._tuplegen, act_on_left=self._act.is_left())

    def specialize(self, new_base_ring=None):
        """
        Return distribution space got by specializing to Sym^k, over
        the new_base_ring.  If new_base_ring is not given, use current
        base_ring.

        EXAMPLES::

        sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        sage: D = Distributions(0, 7, 4); D
        Space of 7-adic distributions with k=0 action and precision cap 4
        sage: D.is_symk()
        False
        sage: D2 = D.specialize(); D2
        Sym^0 Z_7^2
        sage: D2.is_symk()
        True
        sage: D2 = D.specialize(QQ); D2
        Sym^0 Q^2
        """
        if self._character is not None:
            raise NotImplementedError
        if new_base_ring is None:
            new_base_ring = self.base_ring()
        return SymkSpace(k=self._k, base=new_base_ring, tuplegen=self._act._tuplegen, act_on_left=self._act.is_left())

