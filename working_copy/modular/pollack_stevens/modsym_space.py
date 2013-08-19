from sage.misc.misc import verbose
from sage.misc.cachefunc import cached_method
from sage.modules.module import Module
from sage.structure.parent import Parent
from sage.categories.modsym_space_category import ModularSymbolSpaces

from sage.modular.pollack_stevens.fund_domain import ManinRelations
from sage.modular.pollack_stevens.modsym_element import ModularSymbolElement_generic, ModSymAction
from sage.modular.pollack_stevens.sigma0 import Sigma0

class ModularSymbolSpace_generic(Module):
    def __init__(self, group, coefficients, sign=0, element_class = ModularSymbolElement_generic):
        R = coefficients.base_ring()
        Parent.__init__(self, base = R, category = ModularSymbolSpaces(R))
        if sign not in (0,-1,1):
            # sign must be be 0, -1 or 1
            raise ValueError("sign must be 0, -1, or 1")
        self._group = group
        self._coefficients = coefficients
        self.Element = element_class
        self._sign = sign
        self._source = ManinRelations(group.level())
        try:
            action = ModSymAction(Sigma0(self.prime()), self)
        except AttributeError:
            action = ModSymAction(Sigma0(1), self)
        self._populate_coercion_lists_(action_list=[action])
    
    # Category framework
    def coefficient_module(self):
        r"""
        A space of modular symbols is a certain set of homormophisms. This
        returns the target space, which should be in the category
        :class:`MSCoefficientModules` over the same base ring as self.
        """
        return self._coefficients
    
    def weight(self):
        return self.coefficient_module().weight()
    
    def character(self):
        #CURRENT GENERIC DOES NOT CONTAIN THIS
        """
        
        """
        return self.coefficient_module().character()
    
    def source(self):
        """
        
        """
        return self._source
    
    def hecke(self, ell, algorithm = "prep"):
        r"""
        Returns the Hecke operator `T_\ell` as endomorphism of self.
        """
        #Can this be written at this level?
        pass
    
    def _element_constructor(self, data):
        if isinstance(data, ModularSymbolSpace_generic):
            data = data._map
        elif isinstance(data, ManinMap):
            pass
        else:
            # a dict, or a single distribution specifying a constant symbol, etc
            data = ManinMap(self._coefficients, self._source, data)

        if data._codomain != self._coefficients:
            data = data.extend_codomain(self._coefficients)
 
        return self.element_class(data, self, construct=True)
    
    def _repr_(self):
        return "Generic space of modular symbols for %s with sign %s and values in %s"%(self.group(), self.sign(), self.coefficient_module())
    
    def group(self):
        return self._group
    
    def sign(self):
        #ADD THIS TO CATEGORY SETUP?
        return self._sign
    
    def ngens(self):
        return len(self.source().indices())
    
    def number_of_coset_reps(self):
        return len(self.source().reps())
    
    def level(self):
        return self.source().level()
    
    #@cached_method
    def dimension_of_ordinary_subspace(self, p=None, cusp=False):
        """
        If ``cusp`` is ``True``, return dimension of cuspidal ordinary
        subspace. This does a weight 2 computation with sage's ModularSymbols.
        
        EXAMPLES::
        
            sage: M = OverconvergentModularSymbols(11, 0, sign=-1, p=3, prec_cap=4, base=ZpCA(3, 8))
            sage: M.dimension_of_ordinary_subspace()
            2
            sage: M.dimension_of_ordinary_subspace(cusp=True)
            2
            sage: M = OverconvergentModularSymbols(11, 0, sign=1, p=3, prec_cap=4, base=ZpCA(3, 8))
            sage: M.dimension_of_ordinary_subspace(cusp=True)
            2
            sage: M.dimension_of_ordinary_subspace()
            4
            sage: M = OverconvergentModularSymbols(11, 0, sign=0, p=3, prec_cap=4, base=ZpCA(3, 8))
            sage: M.dimension_of_ordinary_subspace()
            6
            sage: M.dimension_of_ordinary_subspace(cusp=True)
            4
            sage: M = OverconvergentModularSymbols(11, 0, sign=1, p=11, prec_cap=4, base=ZpCA(11, 8))
            sage: M.dimension_of_ordinary_subspace(cusp=True)
            1
            sage: M.dimension_of_ordinary_subspace()
            2
            sage: M = OverconvergentModularSymbols(11, 2, sign=1, p=11, prec_cap=4, base=ZpCA(11, 8))
            sage: M.dimension_of_ordinary_subspace(cusp=True)
            0
            sage: M.dimension_of_ordinary_subspace()
            1
            sage: M = OverconvergentModularSymbols(11, 10, sign=1, p=11, prec_cap=4, base=ZpCA(11, 8))
            sage: M.dimension_of_ordinary_subspace(cusp=True)
            1
            sage: M.dimension_of_ordinary_subspace()
            2
        
        An example with odd weight and hence non-trivial character::
        
            sage: K = Qp(11, 6)
            sage: DG = DirichletGroup(11, K)
            sage: chi = DG([K(378703)])
            sage: MM = FamiliesOfOMS(chi, 1, p=11, prec_cap=[4, 4], base_coeffs=ZpCA(11, 4), sign=-1)
            sage: MM.dimension_of_ordinary_subspace()
            1
        """
        try:
            p = self.prime()
        except AttributeError:
            if p is None:
                raise ValueError("If self doesn't have a prime, must specify p.")
        try:
            return self._ord_dim_dict[(p, cusp)]
        except AttributeError:
            self._ord_dim_dict = {}
        except KeyError:
            pass
        from sage.modular.dirichlet import DirichletGroup
        from sage.rings.finite_rings.constructor import GF
        try:
            chi = self.character()
        except AttributeError:
            chi = DirichletGroup(self.level(), GF(p))[0]
        if chi is None:
            chi = DirichletGroup(self.level(), GF(p))[0]
        
        from sage.modular.modsym.modsym import ModularSymbols
        r = self.weight() % (p-1)
        if chi.is_trivial():
            N = chi.modulus()
            if N % p != 0:
                N *= p
            else:
                e = N.valuation(p)
                N.divide_knowing_divisible_by(p ** (e-1))
            chi = DirichletGroup(N, GF(p))[0]
        elif chi.modulus() % p != 0:
            chi = DirichletGroup(chi.modulus() * p, GF(p))(chi)
        DG = DirichletGroup(chi.modulus(), GF(p))
        if r == 0:
            from sage.modular.arithgroup.congroup_gamma0 import Gamma0_constructor as Gamma0
            verbose("in dim: %s, %s, %s"%(self.sign(), chi, p))
            M = ModularSymbols(DG(chi), 2, self.sign(), GF(p))
        else:
            psi = [GF(p)(u) ** r for u in DG.unit_gens()]    #mod p Teichmuller^r
            psi = DG(psi)
            M = ModularSymbols(DG(chi) * psi, 2, self.sign(), GF(p))
        if cusp:
            M = M.cuspidal_subspace()
        hecke_poly = M.hecke_polynomial(p)
        verbose("in dim: %s"%(hecke_poly))
        x = hecke_poly.parent().gen()
        d = hecke_poly.degree() - hecke_poly.ord(x)
        self._ord_dim_dict[(p, cusp)] = d
        return d
    #def prime(self):
    #    """
    #    
    #    """
    #    return self.coefficient_module()._p
    # End category framework
    
    #def _repr_(self):
    #    pass
    
