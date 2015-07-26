from sage.structure.factory import UniqueFactory
from sage.misc.cachefunc import cached_method
from sage.misc.misc import verbose as Verbose
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.padics.factory import Qp
from sage.functions.other import ceil
from sage.functions.log import log
from sage.modular.arithgroup.all import Gamma0
from sage.matrix.constructor import Matrix
from sage.modular.dirichlet import DirichletCharacter
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
        elif isinstance(group, DirichletCharacter):
            character = group
        if coefficients is None:
            coefficients = FamiliesOfOverconvergentDistributions(weight, p, prec_cap, base, base_coeffs, character)
        if isinstance(group, (int, Integer)):
            p = coefficients.prime()
            if group % p != 0:
                group *= p
            group = Gamma0(group)
        elif isinstance(group, DirichletCharacter):
            p = coefficients.prime()
            group = group.modulus()
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
            sage: DD = FamiliesOfOverconvergentDistributions(0, base_coeffs=ZpCA(3, 4), prec_cap=[4,4])
            sage: MM = FamiliesOfOMS(11, 0, coefficients=DD, sign=-1)
            sage: MM.random_element()._consistency_check()
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
        CM = self.coefficient_module().change_precision([M_in, M[1]+1])  ## RP: the +1 here on M[1] is only need for k=0

        R = CM.base_ring()
        
        ## this loop runs thru all of the generators (except (0)-(infty)) and randomly chooses a distribution 
        ## to assign to this generator (in the 2,3-torsion cases care is taken to satisfy the relevant relation)
        D = {}
        t = CM(0)
        for g in gens[1:]:
            Verbose("Looping over generators. At generator %s"%(g))
            #print "CM._prec_cap", CM.precision_cap()
            D[g] = CM.random_element([M_in, M[1]])
            #print "pre:",D[g]
            if g in manin.reps_with_two_torsion():
                if g in manin.reps_with_three_torsion():
                    raise ValueError("Level 1 not implemented")
                Verbose("Generator is two-torsion")
                gamg = manin.two_torsion_matrix(g)
                D[g] = D[g] - D[g] * gamg
                t -= D[g]
            else:
                if g in manin.reps_with_three_torsion():
                    Verbose("Generator is three-torsion")
                    gamg = manin.three_torsion_matrix(g)
                    D[g] = 2*D[g] - D[g] * gamg - D[g] * gamg**2
                    t -= D[g]
                else:
                    Verbose("Generator is non-torsion")
                    t += D[g] * gammas[g] - D[g]
        
        ## Fill in this comment?
        
        a = gam0.matrix()[0,0]
        c = gam0.matrix()[1,0]
        
        shift = 0
        #if CM._character != None:
        #    raise NotImplementedError
        #    chara = CM._character(a)
        #else:
        #    chara = 1
        from sage.modular.pollack_stevens.families_util import automorphy_factor_vector
        R = CM.base_ring()
        Verbose("Compute automorphy factor.")
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
			Verbose("err_val: %s"%(err_val))
            if err_val < 0:
                shift -= err_val
                err = _shift_coeffs(err, shift, right=False)
                t.ordp += shift
            Verbose("shift after err_val: %s"%(shift))
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
        
        Verbose("Solve difference equation.")
        
        #print "t",t
        #Are the following two lines even necessary?
        #        err_pa = err.precision_absolute()
        #        err.reduce_precision_absolute([err_pa[0] - gam_shift - ADD, err_pa[1]])
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
        val, val_vector = mu._valuation(val_vector=True)
        if val < 0:
            Dmu = mu.parent()
            length = Dmu.length_of_moments(M[0])
            mu_val = val_vector[length - 1] + mu.ordp   #This is >= val
            shift -= mu_val
            err.ordp -= mu_val
            Verbose("val = %s, mu_val = %s"%(val, mu_val))
            if val == mu_val:
                mu.ordp -= mu_val
                mu = mu.reduce_precision_absolute(M)
            else:
                from sage.modular.pollack_stevens.coeffmod_OMS_families_element import CoeffMod_OMS_Families_element, _shift_coeffs
                V = Dmu.approx_module(M[0], M[1])
                mu = CoeffMod_OMS_Families_element(V([_shift_coeffs(mu._moments[i], val_vector[length - 1]) for i in range(length)]), Dmu, ordp=mu_val, check=False, var_prec=mu._var_prec)
        Verbose("shift after mu_val: %s"%(shift))
        mu.normalize()
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
        Verbose("mu.valuation: %s, err.valuation: %s"%(mu.valuation(), err.valuation()))
        #ret = self(D).reduce_precision(M)
        #ret = ModSym_OMS_Families_element(D, self, construct=True)
        ret = self(D)
        #for h in gens[1:]:
        #    print h, ret._map[h].ordp
        #t.ordp -= mu_val    #only for Verbose
        #print "Check difference equation (at end): %s"%self.coefficient_module()(mu * gammas[Id] - mu - t)
        if self.sign() == 1:
            return ret.plus_part()
        if self.sign() == -1:
            return ret.minus_part()
        return ret
    
    def random_ordinary_element_with_congruences(self, M=None, data=None):
        if M is None:
            M = self.precision_cap()
        p = self.prime()
        Phis = self.random_element(M)
        for i in range(M[0]+2):
            Phis = Phis.hecke(p)
        if data is None:
            return Phis
        for q, alpha in data:
            for i in range(M[0]+2):
                Phis = Phis.hecke(q) - alpha * Phis
        return Phis

    def is_start_of_basis(self, List):
        r"""
        Determines if the inputed list of OMS families can be extended to a basis of this space.
        More precisely, it checks that the elements of ``List`` are linearly independent modulo
        the uniformizer (by checking the total measures).

        INPUT:

        - ``list`` -- a list of OMS's

        OUTPUT:

        - True/False
        """
        for Phi in List:
            assert Phi.valuation() >= 0, "Symbols must be integral"
        R = self.base().base()
        List = [Phi.list_of_total_measures_at_fixed_weight() for Phi in List]
        d = len(List)
        A = Matrix(R.residue_field(), d, len(List[0]), List)
        Verbose("A =", A)
        return A.rank() == d

    #@cached_method
    def basis_of_ordinary_subspace(self, d=None, verbose=True):
        r"""
        Finds a basis of the ordinary subspace of this space.
    
        INPUT:

        - ``d`` -- (optional) integer equal to the dimension of the ordinary subspace; otherwise this number is just computed via Hida theory

        - ``sign`` -- optional variable which if 1 or -1 restricts to the plus or minus subspace

        OUTPUT:

        - A tuple of OMS families which form the desired basis
        
        EXAMPLES:
        
        sage: MM = FamiliesOfOMS(3, 0, sign=-1, p=3, prec_cap=[4, 4], base_coeffs=ZpCA(3, 4))
        sage: MM.basis_of_ordinary_subspace()
        ()
        
        """
        if d is None:
            d = self.dimension_of_ordinary_subspace()
        try:
            return self._ord_basis
        except AttributeError:
            pass
        basis = []
        done = (d <= len(basis))
        M = self.precision_cap()[0]
        p = self.prime()

        while not done:
            if verbose:
                print "basis has size %s out of predicted %s"%(len(basis),d)
            Verbose("Forming a random symbol")
            if verbose:
                print "-----------------------"
                print "Forming a random symbol"
            Phi = self.random_element()
            val = Phi.valuation()
            if val > 0:
                for key in Phi._map._dict.keys():
                    Phi._map._dict[key].ordp -= val

            Verbose("Projecting to ordinary subspace")
            if verbose:
                print "projecting"
            for a in range(M + 2):
                if verbose:
                    print "%s out of %s"%(a,M+1)
                Phi = Phi.hecke(p)
            ## Should really check here that we are ordinary
            ## Figure out why valuation might increase (example: N=11, p=3, prec>=9)
            val = Phi.valuation()
            if val > 0:
                Verbose("Valuation drop: val=%s"%(val))
                if verbose:
                    print "Valuation drop val=%s"%(val)
                for key in Phi._map._dict.keys():
                    Phi._map._dict[key].ordp -= val

            Verbose("Forming U_p-span of this symbol")
            if verbose:
                print "Forming U_p-span"
            Phi_span = [Phi]
            LI = self.is_start_of_basis(Phi_span)
            if LI and self.is_start_of_basis(basis + [Phi_span[-1]]):
                basis += [Phi_span[-1]]
                Verbose("basis now has size %s"%(len(basis)))
            done = (d <= len(basis))
            while LI and (not done):
                Phi_span += [Phi_span[-1].hecke(p)]
                LI = self.is_start_of_basis(Phi_span)
                if LI and self.is_start_of_basis(basis + [Phi_span[-1]]):
                    basis += [Phi_span[-1]]
                done = (d <= len(basis))
        self._ord_basis = tuple(basis)
        return self._ord_basis
    
    def basis_of_congruence_subspace(self, data=None, verbose=True):
        r"""
        INPUT::
        
            - ``data`` -- a pair whose first entry is the dimension of the
            desired space and whose second entry is a list of pairs `[q, \alpha]`
            where `q` is a prime and `alpha` is an element of the base ring. The
            output will then be a basis of the subspace that is left over after
            killing off the forms where the `q`th Hecke operator acts by
            `\alpha`.
        """
        if data is None:
            return self.basis_of_ordinary_subspace(verbose=verbose)
        q_alpha = [tuple(dat) for dat in data[1]]
        q_alpha.sort()
        q_alpha = tuple(q_alpha)
        try:
            return self._congruence_subspace[q_alpha]
        except AttributeError:
            self._congruence_subspace = {}
        except KeyError:
            pass
        
        d = data[0]
        basis = []
        done = (d <= len(basis))
        M = self.precision_cap()[0]
        p = self.prime()
        
        while not done:
            if verbose:
                print "basis has size %s out of desired %s"%(len(basis),d)
            Verbose("Forming a random symbol")
            if verbose:
                print "-----------------------"
                print "Forming a random symbol"
            Phi = self.random_ordinary_element_with_congruences(data=q_alpha)
            val = Phi.valuation()
            if val > 0:
                for key in Phi._map._dict.keys():
                    Phi._map._dict[key].ordp -= val
            
            Verbose("Forming U_p-span of this symbol")
            if verbose:
                print "Forming U_p-span"
            Phi_span = [Phi]
            LI = self.is_start_of_basis(Phi_span)
            if LI and self.is_start_of_basis(basis + [Phi_span[-1]]):
                basis += [Phi_span[-1]]
                Verbose("basis now has size %s"%(len(basis)))
            done = (d <= len(basis))
            while LI and (not done):
                Phi_span += [Phi_span[-1].hecke(p)]
                LI = self.is_start_of_basis(Phi_span)
                if LI and self.is_start_of_basis(basis + [Phi_span[-1]]):
                    basis += [Phi_span[-1]]
                done = (d <= len(basis))
        self._congruence_subspace[q_alpha] = tuple(basis)
        return self._congruence_subspace[q_alpha]

    def linear_relation_old(self, List, verbose=True):
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
        M = List[0].precision_absolute()[0]
        var_prec = List[0].precision_absolute()[1]

        vs_coef = [[vs[a][b].padded_list()[0] for b in range(len(vs[a]))] for a in range(len(vs))]

        c = find_linear_relation(vs_coef,p,M)
        c = [R(c[a]).add_bigoh(var_prec) for a in range(len(c))]
        #print "c = "
        #print c
        
        if c == []:
            return []
        for j in range(1,var_prec):
            if verbose:
                print "Working at coefficient %s"%(j)
            v = [0 for k in range(len(vs[0]))]
            for i in range(len(vs)):
                for a in range(j):
                    for r in range(len(vs[0])):
                        v[r] += c[i].padded_list()[a] * vs[i][r].padded_list()[j-a]
            #print "j=%s, v=\n%s\n"%(j,v)
            c_coef = find_linear_relation(vs_coef + [v],p,M)
            if len(c_coef) == 0:
                raise ValueError("No linear relation exists.")
            #print "Found relation %s",%(c_coef)
            #print "c_coef[:-1], c_coef[-1] = \n%s\n%s\n"%(c_coef[:-1], c_coef[-1])
            temp = [c_coef[r]/c_coef[len(c_coef)-1] for r in range(len(c_coef)-1)]
            c_coef = temp
            for r in range(len(c)):
                c[r] += c_coef[r] * w**j
            #print "c at %s:\n%s\n"%(j, c)
        return c
    
    def linear_relation(self, List, Psi, verbose=True, prec=None):
        r"""
        INPUT::
        
            - ``List`` -- a list of families of OMSs
            - ``Psi`` -- a family of OMSs
        
        OUTPUT::
        
            - `[v`, `c`]`, where `v` is a vector (over the base ring of self,
            with ``len(List)`` entries) and c is a scalar (over the base ring
            of self) such that
            
            .. MATH::
            
                c\Psi + \sum_{i}v_i\Phi_i = 0,
            
            where `\Phi_i` is the `i`th element of ``List`` and `v_i` is the
            `i`th entry of `v`. If there is a non-trivial such relation, it is
            returned, otherwise the pair [``None``, ``None``] is returned. If ``List``
            has length 0, the first entry of the output is ``None`` while the
            second is 1 or 0 depending on whether `\Psi` is 0 or not.
        """
        for Phi in List:
            assert Phi.valuation() >= 0, "Symbols must be integral"
        assert Psi.valuation() >= 0
        R = self.base()
        Rbase = R.base()
        w = R.gen()
        d = len(List)
        if d == 0:
            if Psi.is_zero():
                return [None, R(1)]
            else:
                return [None, 0]
        if prec is None:
            M, var_prec = Psi.precision_absolute()
        else:
            M, var_prec = prec
        p = self.prime()
        V = R**d
        List_TMs = [Phi.list_of_total_measures() for Phi in List]
        #from sage.structure.sage_object import dumps
        #print repr(dumps(List_TMs))
        List_coef = [[TM.padded_list()[0] if TM != 0 else Rbase(0, M) for TM in i] for i in List_TMs]
        Psi_TMs = Psi.list_of_total_measures()
        Psi_coef = [TM.padded_list()[0] if TM != 0 else Rbase(0, M) for TM in Psi_TMs]
        vectors_TMs = List_TMs + [Psi_TMs]
        vectors_coef = List_coef + [Psi_coef]
        
        v, c = _find_linear_relation(List_coef, Psi_coef, p, M)
        if c is None: #no relation found
            return [None, None]
        v = [R(vv, var_prec) for vv in v]
        ans = v + [R(c, var_prec)]
        #print "ans = "
        #print ans
        length = len(Psi_coef)
        for j in range(1, var_prec):
            if verbose:
                print "Working at coefficient %s"%(j)
            new_vector = [sum(ans[i].padded_list()[a] * vectors_TMs[i][r].padded_list()[j-a] for i in range(d + 1) for a in range(j)) for r in range(length)]
            #print "j=%s, v=\n%s\n"%(j,new_vector)
            higher_order, extra = _find_linear_relation(vectors_coef, new_vector, p, M)
            if extra is None:
                return [None, None]
            #print "higher_order & extra =\n%s\n%s\n"%(higher_order, extra)
            for i in range(d + 1):
                ans[i] += (higher_order[i] / extra) * w**j
            #print "ans at %s=\n%s\n"%(j, ans)
        #print ans
        return [V(ans[:-1]), ans[-1]]
    
    def linear_relation_laurent_series(self, List, Psi, verbose=True, prec=None):
        for Phi in List:
            assert Phi.valuation() >= 0, "Symbols must be integral"
        assert Psi.valuation() >= 0
        R = self.base()
        Rbase = R.base()
        w = R.gen()
        d = len(List)
        if d == 0:
            if Psi.is_zero():
                return [None, R(1)]
            else:
                return [None, 0]
        if prec is None:
            M, var_prec = Psi.precision_absolute()
        else:
            M, var_prec = prec
        p = self.prime()
        V = R**d
        RSR = LaurentSeriesRing(Rbase, R.variable_name())
        VSR = RSR**self.source().ngens()
        List_TMs = [VSR(Phi.list_of_total_measures()) for Phi in List]
        Psi_TMs = VSR(Psi.list_of_total_measures())
        A = Matrix(RSR, List_TMs).transpose()
        try:
            sol = V(A.solve_right(Psi_TMs))
        except ValueError:
            #try "least squares"
            if verbose:
                print "Trying least squares."
            sol = (A.transpose() * A).solve_right(A.transpose() * Psi_TMs)
            #check precision (could make this better, checking each power of w)
            p_prec = M
            diff = Psi_TMs - sum([sol[i] * List_TMs[i] for i in range(len(List_TMs))])
            for i in diff:
                for j in i.list():
                    if p_prec > j.valuation():
                        p_prec = j.valuation()
            if verbose:
                print "p-adic precision is now", p_prec
            #Is this right?
            sol = V([R([j.add_bigoh(p_prec) for j in i.list()]) for i in sol])
        return [sol, R(-1)]
    
    #@cached_method
    def hecke_matrix(self, q, basis=None, verbose=True):
        r"""
        Finds the matrix of T_q wrt to the given basis
    
        INPUT:

        - ``q`` -- a prime

        - ``basis`` -- basis of a T_q-stable subspace

        OUTPUT:

        - d x d matrix where d is the length of basis
        
        EXAMPLES::
        
            sage: MM = FamiliesOfOMS(3, 0, sign=-1, p=3, prec_cap=[4, 4], base_coeffs=ZpCA(3, 4))
            sage: MM.hecke_matrix(7)
            []
        """
        try:
            return self._hecke_matrices[(q, tuple(basis) if basis is not None else None)]
        except AttributeError:
            self._hecke_matrices = {}
        except KeyError:
            pass
        basis_was_None = basis is None
        if basis_was_None:
            basis = self.basis_of_ordinary_subspace(verbose=verbose)
        if len(basis) == 0:
            self._hecke_matrices[()] = Matrix(self.base_ring(), 0)  #Store answer for empty basis independently of q
            return self._hecke_matrices[()]
        basis = list(basis)
        M, var_prec = basis[0].precision_absolute()
        #d = len(basis)
        T = []
        r = 1
        for Phi in basis:
            if verbose:
                print "At %s-th basis element"%(r)
            r = r + 1
            h = Phi.hecke(q)
            extra = None
            p_prec_dec = False
            while extra is None and M > 0:
                row, extra = self.linear_relation(basis, h, verbose=verbose, prec=[M, var_prec])
                if extra is None:
                    print "PRECISION LOSS!!!!"
                    if p_prec_dec:
                        M -= 1
                        if verbose:
                            print "Reducing p-adic precision to %s."%(M)
                        p_prec_dec = False
                    else:
                        var_prec -= 1
                        if verbose:
                            print "Reducing w-adic precision to %s."%(M)
                        p_prec_dec = True
            if extra is None:
                raise ValueError("basis does not span a T_q-stable subspace")
            row = [-a / extra for a in row]
            ## Probably should put some check here that it really worked.
            T.append(row)
        self._hecke_matrices[(q, tuple(basis))] = Matrix(self.base_ring(), T).transpose()
        if basis_was_None:
            self._hecke_matrices[(q, None)] = self._hecke_matrices[(q, tuple(basis))]
        return self._hecke_matrices[(q, tuple(basis))]
    
    #@cached_method
    def hecke_polynomial(self, q, var='x', basis=None, verbose=True):
        r"""
        Note: I (RH) think that the p-adic precision on the base ring of self should be > prec_cap[0] + prec_cap[1]
        because w = pT. (Or using relative precision p-adics could work too I guess).
        EXAMPLES::
        
            sage: MM = FamiliesOfOMS(3, 0, sign=-1, p=3, prec_cap=[4, 4], base_coeffs=ZpCA(3, 8))
            sage: MM.hecke_polynomial(29)
            1 + O(3^4) + O(w^4)
            sage: MM = FamiliesOfOMS(11, 0, sign=-1, p=3, prec_cap=[4, 4], base_coeffs=ZpCA(3, 8))
            sage: MM.hecke_polynomial(3, verbose=False)
            (1 + O(3^8))*x^2 + (2 + 2*3 + 3^2 + O(3^4) + (2*3 + 2*3^2 + O(3^4))*w + O(3^4)*w^2 + (3^3 + O(3^4))*w^3 + O(w^4))*x + 1 + 2*3 + 3^2 + O(3^4) + O(3^4)*w + (2*3^2 + 3^3 + O(3^4))*w^2 + (3^3 + O(3^4))*w^3 + O(w^4)
        """
        #TODO: create a cache for this
        HM = self.hecke_matrix(q, basis, verbose=verbose)
        if HM.nrows() == 0:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = self.base_ring()
            prec_cap = self.precision_cap()
            return PolynomialRing(R ,var)(R(R.base_ring()(1, prec_cap[0]), prec_cap[1]))
        return self.hecke_matrix(q, basis).charpoly(var)
    
    def hecke_polynomial_in_T_variable(self, q, var='x', basis=None, verbose=True):
        r"""
        The function hecke_polynomial returns a polynomial whose coefficients
        are power series in the variable `w`, which represents an element in the
        disc of radius `1/p`. This function instead uses the more standard
        variable `T`, which represents an element in the disc of radius `1`.
        
        EXAMPLES::
        
            sage: MM = FamiliesOfOMS(11, 0, sign=-1, p=3, prec_cap=[4, 4], base_coeffs=ZpCA(3, 8))
            sage: HP = MM.hecke_polynomial_in_T_variable(3, verbose=False); HP
            (1 + O(3^8))*x^2 + (2 + 2*3 + 3^2 + O(3^4) + (2 + 2*3 + O(3^3))*T + O(3^2)*T^2 + (1 + O(3))*T^3 + O(T^4))*x + 1 + 2*3 + 3^2 + O(3^4) + O(3^3)*T + (2 + 3 + O(3^2))*T^2 + (1 + O(3))*T^3 + O(T^4)
        """
        HPw = self.hecke_polynomial(q, var, basis, verbose)
        from sage.rings.power_series_ring import PowerSeriesRing
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        v_prec = self.precision_cap()[1]
        RT = PowerSeriesRing(self.base_ring().base_ring(), 'T', default_prec=v_prec)
        R = PolynomialRing(RT, var)
        poly_coeffs = []
        for c in HPw.padded_list():
            prec = c.prec()
            cL = c.padded_list()
            length = len(cL)
            cL = [cL[i] >> i for i in range(length)]
            j = 0
            while j < length:
                if cL[j].precision_absolute() <= 0:
                    break
                j += 1
            poly_coeffs.append(RT(cL, j))
        poly_coeffs[-1] = RT.one()
        return R(poly_coeffs)
        
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

def find_linear_relation(vs, p, M): #This is the old version. The new one is _find_linear_relation
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
    R = Qp(p, M)
    V = R**d
        
    if d == 1:
        z = True
        for r in range(len(vs[0])):
            if vs[0][r] != 0:
                z = False
        if z:
            return [R(1, M)]
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
        return [R(0, M) for a in range(len(vs)-1)] + [1]
        
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
        
    Verbose("s = %s"%(s))
        
    A = gp(s)
    Verbose("A = %s"%(A))
    if len(vs) == 2:
        A = A.Mat()
        
    s = '['
    v = vs[d-1]
    for r in range(len(v)):
        s += str(ZZ(v[r]))
        if r < len(v) - 1:
            s += ','
    s += ']~'
        
    Verbose("s = %s"%(s))
        
    B = gp(s)
        
    Verbose("B = %s"%(B))

    v = A.mattranspose().matsolvemod(p**M,B)
        
    Verbose("v = %s"%(v))
        
    if v == 0:
        return []
    else:
            ## Move back to SAGE from Pari
        v = [R(v[a], M) for a in range(1,len(v)+1)]
        return v + [R(-1, M)]

def _find_linear_relation(Alist, B, p, M):    #new!
    r"""
    Finds a linear relation between a given list of vectors over Z/p^MZ.  If they are LI, returns an empty list.
    
    INPUT:

    - ``vs`` -- a list of vectors over Z/p^MZ
    - ``p`` -- a prime
    - ``M`` -- positive integer
    
    OUTPUT:

    - A list of p-adic numbers describing the linear relation of the vs
    
    TESTS::
    
        sage: from sage.modular.pollack_stevens.modsym_OMS_families_space import _find_linear_relation
        sage: R = ZpCA(3, 4)
        sage: List = [Sequence([O(3^4), 2*3 + 3^2 + 3^3 + O(3^4), 2 + 3 + 3^2 + 2*3^3 + O(3^4), O(3^4), 1 + 3 + 3^2 + O(3^4), 2 + O(3^4), 1 + 3^2 + 3^3 + O(3^4), 2 + O(3^4), 1 + 3 + 3^2 + O(3^4)], universe=R), Sequence([O(3^4), 1 + 3^2 + O(3^4), 2 + 3 + 3^3 + O(3^4), O(3^4), 1 + 3 + 2*3^2 + 3^3 + O(3^4), 1 + 2*3^2 + O(3^4), 1 + 2*3 + 2*3^2 + 3^3 + O(3^4), 1 + 2*3^2 + O(3^4), 1 + 3 + 2*3^2 + 3^3 + O(3^4)], universe=R), Sequence([O(3^4), 2 + O(3^4), 3 + 3^2 + 3^3 + O(3^4), O(3^4), 2*3 + 3^2 + 3^3 + O(3^4), 1 + 2*3 + 2*3^2 + 3^3 + O(3^4), 3^3 + O(3^4), 1 + 2*3 + 2*3^2 + 3^3 + O(3^4), 2*3 + 3^2 + 3^3 + O(3^4)], universe=R)]
        sage: _find_linear_relation(List[:-1], List[-1], 3, 4)
        ([1 + 3^2 + O(3^4), 2 + 3 + 2*3^2 + O(3^4)], 2 + 2*3 + 2*3^2 + 2*3^3 + O(3^4))
        sage: _find_linear_relation(List[:-1], List[-1], 3, 3)
        ([1 + 3^2 + O(3^3), 2 + 3 + 2*3^2 + O(3^3)], 2 + 2*3 + 2*3^2 + O(3^3))
        sage: List[-1][1] -= 1
        sage: _find_linear_relation(List[:-1], List[-1], 3, 4)
        ([], None)
        sage: List[-1][1] += 1 + 3^3
        sage: _find_linear_relation(List[:-1], List[-1], 3, 4)
        ([], None)
        sage: _find_linear_relation(List[:-1], List[-1], 3, 3)
        ([1 + 3^2 + O(3^3), 2 + 3 + 2*3^2 + O(3^3)], 2 + 2*3 + 2*3^2 + O(3^3))
    """
    vs = list(Alist) + [B]
    d = len(vs)
    R = Qp(p,M)
    V = R**d
        
    if d == 1:
        z = True
        for r in range(len(vs[0])):
            if vs[0][r] != 0:
                z = False
        if z:
            return [R(1, M)]
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
        return [R(0, M) for a in range(len(vs)-1)], 1
        
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
        
    Verbose("s = %s"%(s))
        
    A = gp(s)
    Verbose("A = %s"%(A))
    if len(vs) == 2:
        A = A.Mat()
        
    s = '['
    v = vs[d-1]
    for r in range(len(v)):
        s += str(ZZ(v[r]))
        if r < len(v) - 1:
            s += ','
    s += ']~'
        
    Verbose("s = %s"%(s))
        
    B = gp(s)
        
    Verbose("B = %s"%(B))

    v = A.mattranspose().matsolvemod(p**M,B)
        
    Verbose("v = %s"%(v))
        
    if v == 0:
        return [], None
    else:
            ## Move back to SAGE from Pari
        v = [R(v[a], M) for a in range(1,len(v)+1)]
        return v, R(-1, M)

def random_check(p,N,r,M,var_prec):
    DD = FamiliesOfOverconvergentDistributions(r, base_coeffs=ZpCA(p, M), prec_cap=[M,var_prec])
    MM = FamiliesOfOMS(N, r, coefficients=DD)
    Phis = MM.random_element()
    Phis._consistency_check()

def Iwasawa_invariants(F):
    r"""
    Computes the Iwasawa invariants of the power series `F`. It returns a pair
    `(\mu, \lambda)` where `\mu` is the `\mu`-invariant and `\lambda` is the
    `\lambda`-invariant.
    
    The `\mu`-invariant is is maximum power of `p` dividing the coefficients of
    `F` and the `\lambda`-invariant is the least power of `T` exactly divisible
    by `p^\mu`.
    
    Since the power series is represented by a finite amount of data, what is
    returned is an upper-bound on `\mu` and a lower-bound on `\lambda`.
    
    EXAMPLES::
    
        sage: from sage.modular.pollack_stevens.modsym_OMS_families_space import Iwasawa_invariants
        sage: R = PowerSeriesRing(ZpCA(3, 4), 'T')
        sage: F = R([1,0,1])
        sage: Iwasawa_invariants(F)
        (0, 0)
        sage: F = R([0,0,0])
        sage: Iwasawa_invariants(F)
        (+Infinity, +Infinity)
        sage: F = R([3, 9, 3])
        sage: Iwasawa_invariants(F)
        (1, 0)
        sage: F = R([9, 9, 3])
        sage: Iwasawa_invariants(F)
        (1, 2)
        sage: F = R(0, 3)
        sage: Iwasawa_invariants(F)
        (+Infinity, 3)
        sage: F = R([9,9,9], 4)
        sage: Iwasawa_invariants(F)
        (2, 0)
    """
    from sage.rings.infinity import Infinity
    if F == 0:
        return (Infinity, F.prec())
    Flist = F.padded_list()
    mu = Infinity
    Lambda = -1
    for i in range(len(Flist)):
        val = Flist[i].valuation() if Flist[i] !=0 else Infinity
        if val < mu:
            mu = val
            Lambda = ZZ(i)
    return (mu, Lambda)
