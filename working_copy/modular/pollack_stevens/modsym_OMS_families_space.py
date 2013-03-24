from sage.structure.factory import UniqueFactory
from sage.misc.misc import verbose
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.functions.other import ceil
from sage.functions.log import log
from sage.modular.arithgroup.all import Gamma0
from sage.modular.pollack_stevens.modsym_space import ModularSymbolSpace_generic
from sage.modular.pollack_stevens.coeffmod_OMS_families_space import FamiliesOfOverconvergentDistributions
from sage.modular.pollack_stevens.modsym_OMS_families_element import ModSym_OMS_Families_element

class ModSym_OMS_Families_factory(UniqueFactory):
    def create_key(self, group, weight=None, sign=0, p=None, prec_cap=None, base=None, base_coeffs=None, coefficients=None):
        if sign not in (-1,0,1):
            raise ValueError("sign must be -1, 0, 1")
        if isinstance(group, (int, Integer)):
            character = None
        if coefficients is None:
            #WHERE TO GET CHARACTER?
            coeffcients = FamiliesOfOverconvergentDistributions(weight, p, prec_cap, base, base_coeffs, character)
        if isinstance(group, (int, Integer)):
            p = coefficients.prime()
            if group % p != 0:
                group *= p
            group = Gamma0(group)
        return (group, coefficients, sign)
    
    def create_object(self, version, key):
        return ModSym_OMS_Families_space(*key)

FamiliesOfOMS = ModSym_OMS_Families_factory('ModSym_OMS_Families_space')

class ModSym_OMS_Families_space(ModularSymbolSpace_generic):
    def __init__(self, group, coefficients, sign=0):
        ModularSymbolSpace_generic.__init__(self, group, coefficients, sign=sign, element_class=ModSym_OMS_Families_element)
    
    def _has_coerce_map_from_(self, other):
        if isinstance(other, ModularSymbolSpace_generic):
            if other.group() == self.group() \
                and self.coefficient_module().has_coerce_map_from(other.coefficient_module()):
                return True
        else:
            return False
    
    def _repr_(self):
        return "Space of families of overconvergent modular symbols for %s with sign %s and values in %s"%(self.group(), self.sign(), self.coefficient_module())
    
    def precision_cap(self):
        return self.coefficient_module().precision_cap()
    
    def prime(self):
        return self.coefficient_module().prime()
    
    def change_ring(self, new_base_ring):
        return FamiliesOfOMS(self.group(), coefficients=self.coefficient_module().change_ring(new_base_ring), sign=self.sign())
    
    def _an_element_(self):
        #CHANGE THIS TO AN ACTUAL ELEMENT?
        return self(self.coefficient_module().an_element())
    
    def random_element(self, M=None):
        if M is None:
            M = self.precision_cap()
        #if M == 1?
        p = self.prime()
        k = self.weight()
        M_in = _prec_for_solve_diff_eqn_families(M[0], p)
        #print "M_in", M_in, "var_prec", M[1]
        CM = self.coefficient_module().change_precision([M_in, M[1]])
        R = CM.base_ring()
        manin = self.source()
        gens = manin.gens()
        gammas = manin.gammas
        
        ## this loop runs thru all of the generators (except (0)-(infty)) and randomly chooses a distribution 
        ## to assign to this generator (in the 2,3-torsion cases care is taken to satisfy the relevant relation)
        D = {}
        t = CM(0)
        for g in gens[1:]:
            verbose("Looping over generators. At generator %s"%(g))
            #print "CM._prec_cap", CM.precision_cap()
            D[g] = CM.random_element([M_in, M[1]])
#            print "pre:",D[g]
            if g in manin.reps_with_two_torsion():
                if g in manin.reps_with_three_torsion():
                    raise ValueError("Level 1 not implemented")
                verbose("Generator is two-torsion")
                gamg = manin.two_torsion_matrix(g)
                D[g] = D[g] - D[g] * gamg
                t -= D[g]
            else:
                if g in manin.reps_with_three_torsion():
                    verbose("Generator is three-torsion")
                    gamg = manin.three_torsion_matrix(g)
                    D[g] = 2*D[g] - D[g] * gamg - D[g] * gamg**2
                    t -= D[g]
                else:
                    verbose("Generator is non-torsion")
                    t += D[g] * gammas[g] - D[g]
        
        ## Fill in this comment?

        Id = gens[0]
        if len(gammas) > 1:
            #There is a non-torsion generator
            gam_keys = gammas.keys()
            gam_keys.remove(Id)
            g = gam_keys[0]
            gam = gammas[g]
        else:
            verbose("All generators are torsion.")
            g = gens[1]
            if g in manin.reps_with_two_torsion():
                gam = manin.two_torsion_matrix(g)
            else:
                gam = manin.three_torsion_matrix(g)
        
        a = gam.matrix()[0,0]
        c = gam.matrix()[1,0]
        
        if CM._character != None:
            raise(NotImplementedError)
            chara = CM._character(a)
        else:
            chara = 1
        from sage.modular.pollack_stevens.families_util import automorphy_factor_vector
        R = self.base_ring()
        verbose("Compute automorphy factor.")
        K = automorphy_factor_vector(p, a, c, k, CM._character, M_in, M[1], R)  #Maybe modify aut... to only return 2 first coeffs?
        if k != 0:
            err = -t.moment(0) / (K[0] - 1)
            v = [err] + [R(0)] * (M_in - 1)
            err = CM(v)
        else:
            err = -t.moment(0) / (K[1])
            v = [R(0), err] + [R(0)] * (M_in - 2)
            err = CM(v)
        
        if g in manin.reps_with_two_torsion() or g in manin.reps_with_three_torsion():
            err = err * gam - err
            D[g] += err
            t += err
        else:
            D[g] += err
            t += err * gam - err
        
        verbose("Solve difference equation.")
        mu = t.solve_diff_eqn()
        mu_pr = mu.precision_relative()
        if mu_pr[0] < M[0] or mu_pr[1] < M[1]:
            raise ValueError("Insufficient precision after solving the difference equation.")
        D[Id] = -mu
        ret = self(D)
        if self.sign() == 1:
            return ret.plus_part()
        if self.sign() == -1:
            return ret.minus_part()
        return ret

#@cached_method
def _prec_for_solve_diff_eqn_families(M, p):
    #UPDATE THIS with valuation of K[0]-1 and K[1]
    r"""
        A helper function for determining the (relative) precision of the input
        to solve_diff_eqn required in order obtain an answer with (relative)
        precision ``M``. The parameter ``p`` is the prime and ``k`` is the weight.
        
        Given input precision `M_\text{in}`, the output has precision
        
        .. MATH::
            
            M = M_\text{in} - \lceil\log_p(M_\text{in}) - 3.

    """
    # Do we need the weight?
    # A good guess to begin:
    if M < 1:
        raise ValueError("Desired precision M(=%s) must be at least 1."%(M))
    cp = (p - 2) / (p - 1)
    Min = ZZ(3 + M + ceil(ZZ(M).log(p)))
    # It looks like usually there are no iterations
    # For low M, there can be 1 or 2
    while M > Min*cp - ceil(log((Min * cp),p)) - 3: #THINK ABOUT THIS MORE
        Min += 1
        #print("An iteration in _prec_solve_diff_eqn")
    return Min