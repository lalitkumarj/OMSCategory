import MSCoefficientModule
#from sage.modules.module import Module
from sage.structure.parent import Parent
#from sage.rings.padics.factory import ZpCA
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
#from sage.misc.cachefunc import cached_method
#from sage.categories.action import PrecomposedAction
#from sage.structure.coerce_actions import LeftModuleAction, RightModuleAction
#from sage.matrix.all import MatrixSpace
#from sage.rings.fast_arith import prime_range
#from sage.modular.pollack_stevens.dist import get_dist_classes, Dist_long, iScale
#import sage.rings.ring as ring

class SymkSpace(Parent):
    r"""
    INPUT:
        
    - `k` -- nonnegative integer
    - ``base`` -- ring or None
    - ``
        
    """
    
    def __init__(self, k, base, character, tuplegen, act_on_left):
        Parent.__init__(self,category=MSCoefficientModule, )
        Element = DistributionElementPy
        if hasattr(base, 'prime'):
            p = base.prime()
        else:
            p = None
        DistributionsSpace.__init__(self, k, p, k+1, base, character, tuplegen, act_on_left)

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Symk
            sage: Distributions(0, 5, 10)._repr_()
            'Space of 5-adic distributions with k=0 action and precision cap 10'
            sage: Distributions(0, 5, 10)
            Space of 5-adic distributions with k=0 action and precision cap 10
            sage: Symk(0)
            Sym^0 Q^2
        """
        # TODO: maybe account for character, etc.
        if self.base_ring() is QQ:
            V = 'Q^2'
        elif self.base_ring() is ZZ:
            V = 'Z^2'
        elif isinstance(self.base_ring(), pAdicGeneric) and self.base_ring().degree() == 1:
            if self.base_ring().is_field():
                V = 'Q_%s^2'%(self._p)
            else:
                V = 'Z_%s^2'%(self._p)
        else:
            V = '(%s)^2'%(self.base_ring())
        return "Sym^%s %s"%(self._k, V)

    def change_ring(self, new_base_ring):
        """
        Return a Symk with the same k but a different base ring.

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
        return Symk(k=self._k, base=new_base_ring, character=self._character, tuplegen=self._act._tuplegen, act_on_left=self._act.is_left())

    def lift(self, p=None, M=None, new_base_ring=None):
        """
        Return distribution space that contains lifts with given p,
        precision cap M, and base ring new_base_ring.

        INPUT:

        - `p` -- prime or None
        - `M` -- nonnegative integer or None
        - ``new_base_ring`` -- ring or None

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Symk(0, Qp(7)); D
            Sym^0 Q_7^2
            sage: D.lift(M=20)
            Space of 7-adic distributions with k=0 action and precision cap 20
            sage: D.lift(p=7, M=10)
            Space of 7-adic distributions with k=0 action and precision cap 10
            sage: D.lift(p=7, M=10, new_base_ring=QpCR(7,15)).base_ring()
            7-adic Field with capped relative precision 15
        """
        if self._character is not None:
            # need to change coefficient ring for character
            raise NotImplementedError
        if M is None:
            M = self._prec_cap + 1
        if new_base_ring is None:
            new_base_ring = self.base_ring()
        if p is None:
            try:
                p = new_base_ring.prime()
            except AttributeError:
                raise ValueError("You must specify a prime")
        return DistributionsSpace(k=self._k, p=p, prec_cap=M, base=new_base_ring, character=self._character, tuplegen=self._act._tuplegen, act_on_left=self._act.is_left())
