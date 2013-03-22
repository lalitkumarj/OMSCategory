from sage.structure.factory import UniqueFactory
from sage.misc.misc import verbose
from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.finite_rings.constructor import GF
from sage.functions.other import ceil
from sage.matrix.constructor import Matrix
from sage.modular.arithgroup.all import Gamma0
from sage.modular.pollack_stevens.modsym_space import ModularSymbolSpace_generic
from sage.modular.pollack_stevens.coeffmod_OMS_space import OverconvergentDistributions
from sage.modular.pollack_stevens.modsym_OMS_element import ModSym_OMS_element
from sage.interfaces.gp import gp

class ModSym_OMS_factory(UniqueFactory):
    def create_key(self, group, weight=None, sign=0, p=None, prec_cap=None, base=None, coefficients=None):
        if sign not in (-1,0,1):
            raise ValueError("sign must be -1, 0, 1")
        if isinstance(group, (int, Integer)):
            character = None
        if coefficients is None:
            #WHERE TO GET CHARACTER?
            #OverconvergentDistributions can handle prec_cap and base being None
            coefficients = OverconvergentDistributions(weight, p, prec_cap, base, character)
        if isinstance(group, (int, Integer)):
            p = coefficients.prime()
            if group % p != 0:
                group *= p
            group = Gamma0(group)
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
        r"""
            EXAMPLES::
            
                sage: D = OverconvergentDistributions(0, 7, base=Qp(7,5))
                sage: MS = OverconvergentModularSymbols(14, coefficients=D)
                sage: Phi = MS.random_element()
                sage: Phi._consistency_check()
                This modular symbol satisfies the manin relations
                sage: D = OverconvergentDistributions(2, base=ZpCA(5,10))
                sage: MS = OverconvergentModularSymbols(3, coefficients=D)
                sage: Phi = MS.random_element()
                sage: Phi._consistency_check()
                This modular symbol satisfies the manin relations
        """
        if M is None:
            M = self.precision_cap()
        #if M == 1?
        p = self.prime()
        k = self.weight()
        manin = self.source()
        gens = manin.gens()
        Id = gens[0]
        gammas = manin.gammas
        gam_shift = 0
        if k != 0:
            ## We will now a pick a generator to later alter in order to solve the difference equation
            if len(gammas) > 1:
                verbose("There is a non-torsion generator")
                gam_keys = gammas.keys()
                gam_keys.remove(Id)
                g0 = gam_keys[0]
                gam0 = gammas[g0]
            else:
                verbose("All generators are torsion")
                g0 = gens[1]
                if g0 in manin.reps_with_two_torsion():
                    verbose("Using a two torsion generator")
                    gam0 = manin.two_torsion_matrix(g0)
                else:
                    verbose("Using a three torsion generator")
                    gam0 = manin.three_torsion_matrix(g0)
            a = gam0.matrix()[0,0]
            c = gam0.matrix()[1,0]
            if g0 in manin.reps_with_three_torsion():
                aa = (gam0**2).matrix()[0,0]
                cc = (gam0**2).matrix()[1,0]
                gam_shift = max(c.valuation(p),cc.valuation(p))
            else:
                gam_shift = c.valuation(p)
        
        M_in = _prec_for_solve_diff_eqn(M, p, k) + gam_shift
        verbose("Working with precision %s (M, p, k, gam_shift) = (%s, %s, %s, %s)"%(M_in, M, p, k, gam_shift))
        CM = self.coefficient_module().change_precision(M_in)
        R = CM.base_ring()
        #verbose("M_in, new base ring R = %s, %s"%(M_in, R))
        
        ## this loop runs thru all of the generators (except (0)-(infty)) and randomly chooses a distribution 
        ## to assign to this generator (in the 2,3-torsion cases care is taken to satisfy the relevant relation)
        D = {}
        t = CM(0)
        for g in gens[1:]:
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
                    t += D[g] * gammas[g] - D[g]
        
        ## If k = 0, then t has total measure zero.  However, this is not true when k != 0  
        ## (unlike Prop 5.1 of [PS1] this is not a lift of classical symbol).  
        ## So instead we simply add (const)*mu_1 to some (non-torsion) v[j] to fix this
        ## here since (mu_1 |_k ([a,b,c,d]-1))(trival char) = chi(a) k a^{k-1} c , 
        ## we take the constant to be minus the total measure of t divided by (chi(a) k a^{k-1} c)
        ## Something a little different is done in the two and three torsion case. 

        shift = 0
        if k != 0:
            if CM._character != None:
                chara = CM._character(a)
            else:
                chara = 1
            if not g0 in manin.reps_with_three_torsion():
                err = -t.moment(0)/(chara * k * a**(k-1) * c)
            else:
                err = -t.moment(0)/(chara * k * (a**(k-1) * c + aa**(k-1) * cc))
            err_val = err.valuation()
            if err_val < 0:
                shift -= err_val
                verbose("err_val: %s, shift: %s, err: %s"%(err_val, shift, err))
                err = err << shift
                verbose("New err: %s"%(err))
                t.ordp += shift
            v = [R(0)] * M_in
            v[1] = R(err)
            err = CM(v)

            ## In the two and three torsion case, we now adjust err to make it satisfy the torsion Manin relations
            if g0 in manin.reps_with_two_torsion():
                err = err - err * gam0
                t -= err
            elif g0 in manin.reps_with_three_torsion():
                err = 2 * err - err * gam0 - err * gam0**2
                t -= err
            else:
                t += err * gam0 - err

        verbose("The parent of this dist is %s"%(t.parent()))
        verbose("Before adjusting: %s, %s"%(t.ordp, t._moments))
        #try:
        #    mu = t.solve_diff_eqn()
        #except PrecisionError:
        #    verbose("t"%(t.ordp, t._moments, t.precision_absolute()
        verbose("Shift before mu = %s"%(shift))
        if shift > 0:
            t = t.reduce_precision(t.precision_relative() - k.valuation(p) - gam_shift)
        verbose("About to solve diff_eqn with %s, %s"%(t.ordp, t._moments))
        t.normalize()
        verbose("After normalize: About to solve diff_eqn with %s, %s"%(t.ordp, t._moments))
        mu = t.solve_diff_eqn()
        verbose("Check difference equation (right after): %s"%(mu * gammas[Id] - mu - t))
        mu_val = mu.valuation()
        verbose("mu_val, mu_ordp, mu_moments and mu: %s, %s, %s, %s"%(mu_val, mu.ordp, mu._moments, mu))
        if mu_val < 0:
            shift -= mu_val
            mu.ordp -= mu_val
            if k != 0:
                err.ordp -= mu_val
        verbose("Desired M, mu's M: %s, %s"%(M, mu.precision_relative()))
        verbose("mu.ordp, mu._moments and mu: %s, %s, %s"%(mu.ordp, mu._moments, mu))
        mu = mu.reduce_precision(M)
        mu.normalize()
        verbose("Desired M, mu's M: %s, %s"%(M, mu.precision_relative()))
        verbose("mu.ordp, mu._moments: %s, %s"%(mu.ordp, mu._moments))
        if mu.precision_relative() < M:
            raise ValueError("Insufficient precision after solving the difference equation.")
        D[Id] = -mu
        if shift > 0:
            for h in gens[1:]:
                D[h].ordp += shift
        if k != 0:
            D[g0] += err
        ret = self(D)
        verbose("ret_mu.ordp, ret_mu._moments, ret_mu._prec_rel: %s, %s, %s"%(ret._map[Id].ordp, ret._map[Id]._moments, ret._map[Id].precision_relative()))
        t.ordp -= mu_val    #only for verbose
        verbose("Check difference equation (at end): %s"%(mu * gammas[Id] - mu - t.reduce_precision(M).normalize()))
        return ret

    def is_start_of_basis(self,list):
        r"""
        Determines if the inputed list of OMS's can be extended to a basis of this space

        INPUT:

        - ``list`` -- a list of OMS's

        OUTPUT:

        - True/False
        """
        for j in range(len(list)):
            assert list[j].valuation() >= 0, "Symbols must be integral"

        R = self.base()
        ## THIS IS NOT GREAT
        M = list[0].values()[0].moment(0).precision_relative()
        p = self.prime()
        d = len(list)
            
        list = [list[a].list_of_total_measures() for a in range(len(list))]
        A = Matrix(GF(p),len(list),len(list[0]),list)

        return A.rank() == len(list)


    def linear_relation(self,list):
        r"""
        Finds a linear relation between the given list of OMSs.  If they are LI, returns a list of all 0's.
    
        INPUT:

        - ``list`` -- a list of OMSs

        OUTPUT:

        - A list of p-adic numbers describing the linear relation of the list of OMSs
        """
        for j in range(len(list)):
            assert list[j].valuation() >= 0, "Symbols must be integral"

        R = self.base()
        ## NEED A COMMAND FOR RELATIVE PRECISION OF OMS
        M = list[0].values()[0].precision_relative()
        p = self.prime()
        d = len(list)
        V = R**d
            
        if d == 1:
            if list[0].is_zero():
                return [R(1)]
            else:
                return [R(0)]
        s = '['
        for c in range(len(list)-1):
            v = list[c].list_of_total_measures()
            for r in range(len(v)):
                s += str(ZZ(v[r]))
                if r < len(v) - 1:
                    s += ','
            if c < len(list) - 2:
                s += ';'
        s = s + ']'

        A = gp(s)
        if len(list) == 2:
            A = A.Mat()

        s = '['
        v = list[len(list)-1].list_of_total_measures()
        for r in range(len(v)):
            s += str(ZZ(v[r]))
            if r < len(v) - 1:
                s += ','
        s += ']~'

        B = gp(s)

        v = A.mattranspose().matsolvemod(p**M,B)

        if v == 0:
            return [R(0) for a in range(len(v))]
        else:
            ## Move back to SAGE from Pari
            v = [R(v[a]) for a in range(1,len(v)+1)]
            return v + [R(-1)]

    def basis_of_ordinary_subspace(self,d,sign=0):
        r"""
        Finds a basis of the ordinary subspace of this space.
    
        INPUT:

        - ``d`` -- integer equal to the dimension of the ordinary subspace

        - ``sign`` -- optional variable which if 1 or -1 restricts to the plus or minus subspace

        OUTPUT:

        - A list of OMS's which form the desired basis
        """
        basis = []
        done = (d <= len(basis))
        M = self.precision_cap()
        p = self.prime()

        while not done:
            verbose("Forming a random symbol")
            Phi = self.random_element()
            if sign == 1:
                Phi = Phi.plus_part()
            elif sign == -1:
                Phi = Phi.minus_part()
                verbose("Projecting to ordinary subspace")
            for a in range(M+2):
                Phi = Phi.hecke(p)
            ## Should really check here that we are ordinary
            verbose("Forming U_p-span of this symbol")
            Phi_span = [Phi]
            LI = self.is_start_of_basis(Phi_span)
            if self.is_start_of_basis(basis + [Phi_span[len(Phi_span)-1]]) and LI:
                basis += [Phi_span[len(Phi_span)-1]]
                verbose("basis now has size %s"%(len(basis)))
            done = (d <= len(basis))
            while LI and (not done):
                Phi_span += [Phi_span[len(Phi_span)-1].hecke(p)]
                LI = self.is_start_of_basis(Phi_span)
                if self.is_start_of_basis(basis + [Phi_span[len(Phi_span)-1]]) and LI:
                    basis += [Phi_span[len(Phi_span)-1]]
                done = (d <= len(basis))
            
        return basis
                 

    def hecke_matrix(self,q,basis):
        r"""
        Finds the matrix of T_q wrt to the given basis
    
        INPUT:

        - ``q`` -- a prime

        - ``basis`` -- basis of a T_q-stable subspace

        OUTPUT:

        - d x d matrix where d is the length of basis
        """
        d = len(basis)
        T = []
        for r in range(d):
            h = basis[r].hecke(q)
            row = self.linear_relation(basis + [h])
            ## Probably should put some check here that it really worked.
            row.pop()
            T += [row]
        return Matrix(T).transpose()
            
    
#@cached_method
def _prec_for_solve_diff_eqn(M, p, k):
    r"""
        A helper function for determining the (relative) precision of the input
        to solve_diff_eqn required in order obtain an answer with (relative)
        precision ``M``. The parameter ``p`` is the prime and ``k`` is the weight.
        
        Given input precision `M_\text{in}`, the output has precision
        
        .. MATH::
            
            M = M_\text{in} - \lceil\log_p(M_\text{in}) - 3 - v_p(k),
        
        where the latter term only appears if ``k`` is not 0.
        
        ::EXAMPLES:
        
            sage: from sage.modular.pollack_stevens.modsym_OMS_space import _prec_for_solve_diff_eqn
            sage: [_prec_for_solve_diff_eqn(M, p, k) for p in [2,3,11] for M in [1,3,10,20] for k in [0, 2, 6]]
            [7, 8, 8, 10, 11, 11, 18, 19, 19, 28, 29, 29, 6, 6, 7, 8, 8, 9, 16, 16, 17, 26, 26, 27, 5, 5, 5, 7, 7, 7, 15, 15, 15, 25, 25, 25]
    """
    # A good guess to begin:
    val_k = ZZ(k).valuation(p) if k != 0 else 0
    if M < 1:
        raise ValueError("Desired precision M(=%s) must be at least 1."%(M))
    Min = ZZ(3 + M + ceil(ZZ(M).log(p)) + val_k)
    # It looks like usually there are no iterations
    # For low M, there can be 1 or 2
    while M > Min - ceil(Min.log(p)) - 3 - val_k:
        Min += 1
        #print("An iteration in _prec_solve_diff_eqn")
    return Min
