from sage.structure.factory import UniqueFactory
from sage.misc.misc import verbose
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.padics.factory import Qp
from sage.functions.other import ceil
from sage.functions.log import log
from sage.modular.arithgroup.all import Gamma0
from sage.matrix.constructor import Matrix
from sage.modular.pollack_stevens.modsym_space import ModularSymbolSpace_generic
from sage.modular.pollack_stevens.coeffmod_OMS_families_space import FamiliesOfOverconvergentDistributions
from sage.modular.pollack_stevens.modsym_OMS_families_element import ModSym_OMS_Families_element
from sage.modular.pollack_stevens.modsym_OMS_space import _prec_for_solve_diff_eqn
from sage.interfaces.gp import gp

class ModSym_OMS_Families_factory(UniqueFactory):
    """
    TESTS::
    
        sage: D = FamiliesOfOverconvergentDistributions(4, prec_cap=[5,3], base_coeffs=ZpCA(11))
        sage: MS = FamiliesOfOMS(11, coefficients=D)
    """
    def create_key(self, group, weight=None, sign=0, p=None, prec_cap=None, base=None, base_coeffs=None, coefficients=None):
        if sign not in (-1,0,1):
            raise ValueError("sign must be -1, 0, 1")
        if isinstance(group, (int, Integer)):
            character = None
        if coefficients is None:
            #WHERE TO GET CHARACTER?
            coefficients = FamiliesOfOverconvergentDistributions(weight, p, prec_cap, base, base_coeffs, character)
        if isinstance(group, (int, Integer)):
            p = coefficients.prime()
            if group % p != 0:
                group *= p
            group = Gamma0(group)
        return (group, coefficients, sign)
    
    def create_object(self, version, key):
        return ModSym_OMS_Families_space(*key)

FamiliesOfOMS = ModSym_OMS_Families_factory('FamiliesOfOMS')

class ModSym_OMS_Families_space(ModularSymbolSpace_generic):
    def __init__(self, group, coefficients, sign=0):
        """
        TEST::
        
            sage: from sage.modular.pollack_stevens.modsym_OMS_families_space import FamiliesOfOMS
            sage: MS = FamiliesOfOMS(14, 0, -1, 3, [5, 4], base_coeffs=ZpCA(3))
            sage: TestSuite(MS).run()
        """
        ModularSymbolSpace_generic.__init__(self, group, coefficients, sign=sign, element_class=ModSym_OMS_Families_element)
    
    def _has_coerce_map_from_(self, other):
        if isinstance(other, ModularSymbolSpace_generic):
            if other.group() == self.group() \
                and self.coefficient_module().has_coerce_map_from(other.coefficient_module()):
                return True
        else:
            return False
    
    #def __reduce__(self):
    #    return FamiliesOfOMS.reduce_data(self)
    
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
        r"""
        EXAMPLES::
        
            sage: DD = FamiliesOfOverconvergentDistributions(4, base_coeffs=Qp(11, 6), prec_cap=[6,4])
            sage: MM = FamiliesOfOMS(11, 4, coefficients=DD)
            sage: Phis = MM.random_element()
            sage: Phis._consistency_check()
            This modular symbol satisfies the manin relations
            sage: DD = FamiliesOfOverconvergentDistributions(0, base_coeffs=ZpCA(11, 4), prec_cap=[4,4])
            sage: MM = FamiliesOfOMS(11, 0, coefficients=DD, sign=-1)
            sage: Phis = MM.random_element()
            sage: Phis._consistency_check()
            This modular symbol satisfies the manin relations
        """
        if M is None:
            M = self.precision_cap()
        #if M == 1?
        p = self.prime()
        k = self.weight()
        manin = self.source()
        gens = manin.gens()
        gammas = manin.gammas
        Id = gens[0]
        g0, gam0, gam_shift = manin._nice_gamma(p, k)
        
        # RP: _prec_for_solve... isn't working right
        #M_in = _prec_for_solve_diff_eqn_families(M[0], p)
        #print "M_in", M_in
        #M_in += gam_shift
        #print "M_in", M_in
        ADD = 0    # Trying to fix precision issue!
        #M_in = ZZ(1 + M[0] + ceil(ZZ(M[0]).log(p))) + gam_shift + ADD  #fix this
        M_in = _prec_for_solve_diff_eqn(M[0], p) + gam_shift + ADD
        #print "M[0]", M[0], "M_in", M_in, "var_prec", M[1]
        #print "We'l get", M_in - 1 - ceil()
        CM = self.coefficient_module().change_precision([M_in, M[1]+1])

        R = CM.base_ring()
        
        ## this loop runs thru all of the generators (except (0)-(infty)) and randomly chooses a distribution 
        ## to assign to this generator (in the 2,3-torsion cases care is taken to satisfy the relevant relation)
        D = {}
        t = CM(0)
        for g in gens[1:]:
            verbose("Looping over generators. At generator %s"%(g))
            #print "CM._prec_cap", CM.precision_cap()
            D[g] = CM.random_element([M_in, M[1]])
            #print "pre:",D[g]
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
        
        a = gam0.matrix()[0,0]
        c = gam0.matrix()[1,0]
        
        shift = 0
        if CM._character != None:
            raise NotImplementedError
            chara = CM._character(a)
        else:
            chara = 1
        from sage.modular.pollack_stevens.families_util import automorphy_factor_vector
        R = CM.base_ring()
        verbose("Compute automorphy factor.")
        K = automorphy_factor_vector(p, a, c, k, CM._character, CM.length_of_moments(M_in), M_in, R)  #Maybe modify aut... to only return 2 first coeffs?
        #K = automorphy_factor_vector(p, a, c, k, CM._character, M_in, M[1], R) #Should this be it instead

        if g0 in manin.reps_with_three_torsion():
            aa = (gam0**2).matrix()[0,0]
            cc = (gam0**2).matrix()[1,0]
            KK = automorphy_factor_vector(p, aa, cc, k, CM._character, CM.length_of_moments(M_in), M_in, R) 
            
        if k != 0:
            if g0 in manin.reps_with_three_torsion():
                err = -t.moment(0)/ (2 - K[0] - KK[0])
            else:
                err = -t.moment(0) / (K[0] - 1)
            v = [err] + [R(0)] * (CM.length_of_moments(M_in) - 1)
            err = CM(v)
        else:
            from sage.modular.pollack_stevens.coeffmod_OMS_families_element import _padic_val_of_pow_series, _shift_coeffs
            if g0 in manin.reps_with_three_torsion():
                err = - t.moment(0) / (K[1] + KK[1])
            else:
                err = -t.moment(0) / (K[1])
            err_val = _padic_val_of_pow_series(err) ###
            if err_val < 0:
                shift -= err_val
                err = _shift_coeffs(err, shift, right=False)
                t.ordp += shift
            v = [R(0), err] + [R(0)] * (CM.length_of_moments(M_in) - 2)
            err = CM(v)
        
        if g0 in manin.reps_with_two_torsion():
            err = err - err * gam0
            t -= err
        elif g0 in manin.reps_with_three_torsion():
            err = 2 * err - err * gam0 - err * gam0**2
            t -= err
        else:
            t += err * gam0 - err
        
        verbose("Solve difference equation.")

        #print "t",t
        #Are the following two lines even necessary?
        err_pa = err.precision_absolute()
        err.reduce_precision_absolute([err_pa[0] - gam_shift - ADD, err_pa[1]])
        #if shift > 0:
        #    t_pr = t.precision_relative()
        #    t = t.reduce_precision([t_pr[0] - gam_shift - ADD, t_pr[1] - ADD])
        #else:
        #    t_pr = t.precision_relative()
        #    t = t.reduce_precision([t_pr[0] - ADD, t_pr[1] - ADD])
        t_pa = t.precision_absolute()
        t = t.reduce_precision_absolute([t_pa[0] - gam_shift - ADD, M[1]])

        t.normalize()
        #print "M_in =", M_in
        #print "t[0] =", t.moment(0)
        mu = t.solve_diff_eqn()
        mu_val = mu.valuation()
        #print "Shift:", shift
        #print "mu_val:", mu_val
        if mu_val < 0:
            shift -= mu_val
            mu.ordp -= mu_val
            err.ordp -= mu_val
        mu.reduce_precision_absolute(M)
        mu.normalize()
        # RP: commented out these lines as precision isn't set up to work properly yet
        #        if mu_pr[0] < M[0] or mu_pr[1] < M[1]:
        #            raise ValueError("Insufficient precision after solving the difference equation.")
        #print "mu.ordp:", mu.ordp
        D[Id] = -mu
        #print "D[Id]:",
        #print "Shift:", shift
        #for h in gens[1:]:
        #    print h, D[h].ordp
        if shift > 0:
            for h in gens[1:]:
                D[h].ordp += shift
        #for h in gens[1:]:
        #    print h, D[h].ordp
        D[g0] += err
        #ret = self(D).reduce_precision(M)
        #ret = ModSym_OMS_Families_element(D, self, construct=True)
        ret = self(D)
        #for h in gens[1:]:
        #    print h, ret._map[h].ordp
        #t.ordp -= mu_val    #only for verbose
        #print "Check difference equation (at end): %s"%self.coefficient_module()(mu * gammas[Id] - mu - t)
        if self.sign() == 1:
            return ret.plus_part()
        if self.sign() == -1:
            return ret.minus_part()
        return ret

    def is_start_of_basis(self, list):
        r"""
        Determines if the inputed list of OMS families can be extended to a basis of this space

        INPUT:

        - ``list`` -- a list of OMS's

        OUTPUT:

        - True/False
        """
        for Phi in list:
            assert Phi.valuation() >= 0, "Symbols must be integral"
        R = self.base().base()
        list = [Phi.list_of_total_measures_at_fixed_weight() for Phi in list]
        d = len(list)
        A = Matrix(R.residue_field(), d, len(list[0]), list)
        return A.rank() == d

    #@cached_method
    def basis_of_ordinary_subspace(self, d=None):
        r"""
        Finds a basis of the ordinary subspace of this space.
    
        INPUT:

        - ``d`` -- (optional) integer equal to the dimension of the ordinary subspace; otherwise this number is just computed via Hida theory

        - ``sign`` -- optional variable which if 1 or -1 restricts to the plus or minus subspace

        OUTPUT:

        - A list of OMS families which form the desired basis
        """
        if d is None:
            d = self.dimension_of_ordinary_subspace()
        basis = []
        done = (d <= len(basis))
        M = self.precision_cap()[0]
        p = self.prime()

        while not done:
            print "basis has size %s out of predicted %s"%(len(basis),d)
            verbose("Forming a random symbol")
            print "-----------------------"
            print "Forming a random symbol"
            Phi = self.random_element()

            verbose("Projecting to ordinary subspace")
            print "projecting"
            for a in range(M + 2):
                print "%s out of %s"%(a,M+1)
                Phi = Phi.hecke(p)
            ## Should really check here that we are ordinary

            verbose("Forming U_p-span of this symbol")
            print "Forming U_p-span"
            Phi_span = [Phi]
            LI = self.is_start_of_basis(Phi_span)
            if LI and self.is_start_of_basis(basis + [Phi_span[-1]]):
                basis += [Phi_span[-1]]
                verbose("basis now has size %s"%(len(basis)))
            done = (d <= len(basis))
            while LI and (not done):
                Phi_span += [Phi_span[-1].hecke(p)]
                LI = self.is_start_of_basis(Phi_span)
                if LI and self.is_start_of_basis(basis + [Phi_span[-1]]):
                    basis += [Phi_span[-1]]
                done = (d <= len(basis))
        return basis

    def linear_relation(self, List):
        r"""
        Finds a linear relation between the given list of OMSs.  If they are LI, returns a list of all 0's.
    
        INPUT:

        - ``List`` -- a list of OMSs

        OUTPUT:

        - A list of power series describing the linear relation of the list of OMS families
        """
        for Phi in List:
            assert Phi.valuation() >= 0, "Symbols must be integral"

        vs = [Phi.list_of_total_measures() for Phi in List]
        R = vs[0][0].parent()
        w = R.gen()
        p = self.prime()
        DD = self.coefficient_module()
        M = List[0].precision_absolute()[0]
        var_prec = List[0].precision_absolute()[1]

        vs_coef = [[vs[a][b].padded_list()[0] for b in range(len(vs[a]))] for a in range(len(vs))]

        c = find_linear_relation(vs_coef,p,M)
        c = [R(c[a]).add_bigoh(var_prec) for a in range(len(c))]
        
        if c == []:
            return []
        for j in range(1,var_prec):
            print "Working at coefficient %s"%(j)
            v = [0 for k in range(len(vs[0]))]
            for i in range(len(vs)):
                for a in range(j):
                    for r in range(len(vs[0])):
                        v[r] += c[i].padded_list()[a] * vs[i][r].padded_list()[j-a]
            c_coef = find_linear_relation(vs_coef + [v],p,M)
            #print "Found relation %s",%(c_coef)
            temp = [c_coef[r]/c_coef[len(c_coef)-1] for r in range(len(c_coef)-1)]
            c_coef = temp
            for r in range(len(c)):
                c[r] += c_coef[r] * w**j

        return c
    
    def linear_relation_new(self, List):
        pass
    
    def hecke_matrix(self, q, basis):
        r"""
        Finds the matrix of T_q wrt to the given basis
    
        INPUT:

        - ``q`` -- a prime

        - ``basis`` -- basis of a T_q-stable subspace

        OUTPUT:

        - d x d matrix where d is the length of basis
        """
        #d = len(basis)
        T = []
        r = 1
        for Phi in basis:
            print "At %s-th basis element"%(r)
            r = r + 1
            h = Phi.hecke(q)
            row = self.linear_relation(basis + [h])
            row = [-row[a]/row[len(row)-1] for a in range(len(row)-1)]
            ## Probably should put some check here that it really worked.
            T.append(row)
        return Matrix(T).transpose()

            






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

def find_linear_relation(vs, p, M):
    r"""
    Finds a linear relation between a given list of vectors over Z/p^MZ.  If they are LI, returns an empty list.
    
    INPUT:

    - ``vs`` -- a list of vectors over Z/p^MZ
    - ``p`` -- a prime
    - ``M`` -- positive integer
    
    OUTPUT:

    - A list of p-adic numbers describing the linear relation of the vs
    """
    d = len(vs)
    R = Qp(p,M)
    V = R**d
        
    if d == 1:
        z = True
        for r in range(len(vs[0])):
            if vs[0][r] != 0:
                z = False
        if z:
            return [R(1)]
        else:
            return []
        # Would be better to use the library as follows. Unfortunately, matsolvemod
        # is not included in the gen class!!
        #from sage.libs.pari.gen import pari
        #cols = [List[c].list_of_total_measures() for c in range(len(List) - 1)]
        #A = pari(Matrix(ZZ, len(cols[0]), d - 1, lambda i, j : cols[j][i].lift()))
        #aug_col = pari([ZZ(a) for a in List[-1].list_of_total_measures()]).Col()
        #v = A.matsolvemod(p ** M, aug_col)

        ## Stupid hack here to deal with the fact that I can't get
        ## matsolvemod to work if B=0
    z = True
    for r in range(len(vs[d-1])):
        if vs[d-1][r] != 0:
            z = False

    if z:
        return [R(0) for a in range(len(vs)-1)] + [1]
        
    s = '['
    for c in range(d-1):
        v = vs[c]
        for r in range(len(v)):
            s += str(ZZ(v[r]))
            if r < len(v) - 1:
                s += ','
        if c < d - 2:
            s += ';'
    s = s + ']'
        
    verbose("s = %s"%(s))
        
    A = gp(s)
    verbose("A = %s"%(A))
    if len(vs) == 2:
        A = A.Mat()
        
    s = '['
    v = vs[d-1]
    for r in range(len(v)):
        s += str(ZZ(v[r]))
        if r < len(v) - 1:
            s += ','
    s += ']~'
        
    verbose("s = %s"%(s))
        
    B = gp(s)
        
    verbose("B = %s"%(B))

    v = A.mattranspose().matsolvemod(p**M,B)
        
    verbose("v = %s"%(v))
        
    if v == 0:
        return []
    else:
            ## Move back to SAGE from Pari
        v = [R(v[a]) for a in range(1,len(v)+1)]
        return v + [R(-1)]

def random_check(p,N,r,M,var_prec):
    DD = FamiliesOfOverconvergentDistributions(r, base_coeffs=ZpCA(p, M), prec_cap=[M,var_prec])
    MM = FamiliesOfOMS(N, r, coefficients=DD)
    Phis = MM.random_element()
    Phis._consistency_check()
    
