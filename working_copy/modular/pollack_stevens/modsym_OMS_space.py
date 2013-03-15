from sage.structure.factory import UniqueFactory
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.functions.other import ceil
from sage.modular.arithgroup.all import Gamma0
from sage.modular.pollack_stevens.modsym_space import ModularSymbolSpace_generic
from sage.modular.pollack_stevens.coeffmod_OMS_space import OverconvergentDistributions
from sage.modular.pollack_stevens.modsym_OMS_element import ModSym_OMS_element

class ModSym_OMS_factory(UniqueFactory):
    def create_key(self, group, weight=None, sign=0, p=None, prec_cap=None, base=None, coefficients=None):
        if sign not in (-1,0,1):
            raise ValueError("sign must be -1, 0, 1")
        if isinstance(group, (int, Integer)):
            group = Gamma0(group)
        if coefficients is None:
            #WHERE TO GET CHARACTER?
            coeffcients = OverconvergentDistributions(weight, p, prec_cap, base, character)
        return (group, coefficients, sign)
    
    def create_object(self, version, key):
        return ModSym_OMS_space(*key)

OverconvergentModularSymbols = ModSym_OMS_factory('ModSym_OMS_space')

class ModSym_OMS_space(ModularSymbolSpace_generic):
    def __init__(self, group, coefficients, sign=0):
        ModularSymbolSpace_generic.__init__(self, group, coefficients, sign=sign, element_class=ModSym_OMS_element)
    
    def _has_coerce_map_from_(self, other):
        if isinstance(other, ModularSymbolSpace_generic):
            if other.group() == self.group() \
                and self.coefficient_module().has_coerce_map_from(other.coefficient_module()):
                return True
        else:
            return False
    
    def _repr_(self):
        return "Space of overconvergent modular symbols for %s with sign %s and values in %s"%(self.group(), self.sign(), self.coefficient_module())
    
    def precision_cap(self):
        return self.coefficient_module().precision_cap()
    
    def prime(self):
        return self.coefficient_module().prime()
    
    def change_ring(self, new_base_ring):
        return OverconvergentModularSymbols(self.group(), coefficients=self.coefficient_module().change_ring(new_base_ring), sign=self.sign())
    
    def _an_element_(self):
        #CHANGE THIS TO AN ACTUAL ELEMENT?
        return self(self.coefficient_module().an_element())
    
    def random_element(self, M=None):
        if M is None:
            M = self.precision_cap()
        #if M == 1?
        p = self.prime()
        k = self.weight()
        M_in = _prec_for_solve_diff_eqn(M, p)
        CM = self.coefficient_module().change_precision(M_in)
        R = CM.base_ring()
        manin = self.source()
        
        ## this loop runs thru all of the generators (except (0)-(infty)) and randomly chooses a distribution 
        ## to assign to this generator (in the 2,3-torsion cases care is taken to satisfy the relevant relation)
        D = {}
        t = CM(0)
        for g in manin.gens()[1:]:
            #print "CM._prec_cap", CM.precision_cap()
            D[g] = CM.random_element(M_in)
#            print "pre:",D[g]
            if g in manin.reps_with_two_torsion():
                if g in manin.reps_with_three_torsion():
                    raise ValueError("Level 1 not implemented")
                gamg = manin.two_torsion_matrix(g)
                D[g] = D[g] - D[g] * gamg
                t -= D[g]
            else:
                if g in manin.reps_with_three_torsion():
                    gamg = manin.three_torsion_matrix(g)
                    D[g] = 2*D[g] - D[g] * gamg - D[g] * gamg**2
                    t -= D[g]
                else:
                    t += D[g] * manin.gammas[g] - D[g]
        
        ## If k = 0, then t has total measure zero.  However, this is not true when k != 0  
        ## (unlike Prop 5.1 of [PS1] this is not a lift of classical symbol).  
        ## So instead we simply add (const)*mu_1 to some (non-torsion) v[j] to fix this
        ## here since (mu_1 |_k ([a,b,c,d]-1))(trival char) = chi(a) k a^{k-1} c , 
        ## we take the constant to be minus the total measure of t divided by (chi(a) k a^{k-1} c)

        if k != 0:
            #TODO: simplify this by having the ManinRelations object compute once and for all a non-torsion generator if it exists (or does Lalit's fix make this all uneccessary?)
            j = 1
            g = manin.gens()[j]
            while (g in manin.reps_with_two_torsion()) or (g in manin.reps_with_three_torsion()) and (j < len(manin.gens())):
                j += 1
                g = manin.gens()[j]
            if j == len(manin.gens()):
                raise ValueError("everything is 2 or 3 torsion!  NOT YET IMPLEMENTED IN THIS CASE")

            gam = manin.gammas[g]
            a = gam.matrix()[0,0]
            c = gam.matrix()[1,0]

            if CM._character != None:
                chara = CM._character(a)
            else:
                chara = 1
            err = -t.moment(0)/(chara*k*a**(k-1)*c)
            v = [R(0)] * M_in
            v[1] = R(1)
            mu_1 = err * CM(v)
            D[g] += mu_1
#            print "Modifying: ",D[g]
            t = t + mu_1 * gam - mu_1
        
        Id = manin.gens()[0]
        mu = t.solve_diff_eqn()
        D[Id] = -mu
        return self(D)

#@cached_method
def _prec_for_solve_diff_eqn(M, p):
    r"""
        A helper function for determining the (relative) precision of the input
        to solve_diff_eqn required in order obtain an answer with (relative)
        precision ``M``. The parameter ``p`` is the prime and ``k`` is the weight.
        
        Given input precision `M_\text{in}`, the output has precision
        
        .. MATH::
            
            M = M_\text{in} - \lceil\log_p(M_\text{in}) - 3.
        
        ::EXAMPLES:
        
            sage: [_prec_for_solve_diff_eqn(M, p) for p in [2,3,5] for M in [1,3,10,20]]
            [7, 10, 18, 28, 6, 8, 16, 26, 5, 8, 15, 25]

    """
    # Do we need the weight?
    # A good guess to begin:
    if M < 1:
        raise ValueError("Desired precision M(=%s) must be at least 1."%(M))
    Min = ZZ(3 + M + ceil(ZZ(M).log(p)))
    # It looks like usually there are no iterations
    # For low M, there can be 1 or 2
    while M > Min - ceil(Min.log(p)) - 3:
        Min += 1
        #print("An iteration in _prec_solve_diff_eqn")
    return Min