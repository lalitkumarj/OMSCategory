# This should be cythoned once it's done.

#from sage.rings.power_series_ring import PowerSeriesRing

from copy import copy, deepcopy   #Not necessary in cython version
from sage.misc.misc import verbose
from sage.rings.infinity import Infinity
from sage.rings.padics.precision_error import PrecisionError
from sage.modular.pollack_stevens.coeffmod_element import CoefficientModuleElement_generic
from sage.modular.pollack_stevens.coeffmod_OMS_element import CoeffMod_OMS_element
from sage.modular.pollack_stevens.modsym_OMS_space import _prec_for_solve_diff_eqn
from sage.rings.integer_ring import ZZ
from sage.rings.padics.factory import Qp
from sage.functions.other import ceil

class CoeffMod_OMS_Families_element(CoefficientModuleElement_generic):
    r"""
        Fill this in.
        
        EXAMPLES::
        
            sage: D8 = FamiliesOfOverconvergentDistributions(0, prec_cap = [8 ,4], base_coeffs=ZpCA(3, 8))
            sage: mu8 = D8([1,2,3,4,5,6,7,8,9,10]); mu8 #the before last moment should be O(3) but power series don't work well
            (1 + O(3^5), 2 + O(3^5), 3 + O(3^4), 1 + 3 + O(3^4), 2 + 3 + O(3^3), 2*3 + O(3^3), 1 + 2*3 + O(3^2), 2 + 2*3 + O(3^2), 0, 1 + O(3)) + O(w^4)
            sage: D4 = FamiliesOfOverconvergentDistributions(0, prec_cap = [4 ,4], base_coeffs=ZpCA(3, 4))
            sage: mu4 = D4([1,2,3,4,5,6,7,8,9,10]); mu4
            (1 + O(3^4), 2 + O(3^4), 3 + O(3^3), 1 + 3 + O(3^3), 2 + 3 + O(3^2), 2*3 + O(3^2), 1 + O(3), 2 + O(3)) + O(w^4)
            sage: D4(mu8)
            (1 + O(3^4), 2 + O(3^4), 3 + O(3^3), 1 + 3 + O(3^3), 2 + 3 + O(3^2), 2*3 + O(3^2), 1 + O(3), 2 + O(3)) + O(w^4)
            sage: mu4 == D4(mu8)
            True
            sage: D42 = FamiliesOfOverconvergentDistributions(0, prec_cap = [4 ,2], base_coeffs=ZpCA(3, 4))
            sage: mu42 = D42([1,2,3,4,5,6,7,8,9,10]); mu42 
            (1 + O(3^4), 2 + O(3^4), 3 + O(3^3), 1 + 3 + O(3^3), 2 + 3 + O(3^2), 2*3 + O(3^2), 1 + O(3), 2 + O(3)) + O(w^2)
            sage: D42(mu8)
            (1 + O(3^4), 2 + O(3^4), 3 + O(3^3), 1 + 3 + O(3^3), 2 + 3 + O(3^2), 2*3 + O(3^2), 1 + O(3), 2 + O(3)) + O(w^2)
            sage: mu42 == D42(mu8)
            True
            sage: D42(15)   #the last moment should be O(3) but power series don't work well
            3 * (2 + O(3), 0) + O(w^2)
            sage: D = FamiliesOfOverconvergentDistributions(2, prec_cap = [4 ,4], base_coeffs=ZpCA(11, 4))
            sage: R = D.base_ring(); K = R.base_extend(R.base_ring().fraction_field())
            sage: v = [K([1,2,11]) / 11, K([1]), K([11,1,1])]
            sage: D(v)
            11^-1 * (1 + O(11^2) + (2 + O(11^2))*w + (11 + O(11^2))*w^2, 11 + O(11^2), 0) + O(w^4)
            sage: D(0)
            11^4 * () + O(w^4)
    
    TEST::
    
        sage: TestSuite(mu4).run()
    """
    def __init__(self, moments, parent, ordp=0, check=True, var_prec=None):
        CoefficientModuleElement_generic.__init__(self, parent)
        #TODO: for sage 6.2+, will need to make parent the argument after self
        #TODO: finish this
        #if var_prec == None:
        #    var_prec = parent._prec_cap[1]
        #self._var_prec = var_prec
        if check:
            # Need to check/construct: moments, ordp, and var_prec
            par_p_prec, par_v_prec = parent.precision_cap()
            if var_prec is None:
                var_prec = par_v_prec
            else:
                var_prec = min(var_prec, par_v_prec)
            if isinstance(moments, CoeffMod_OMS_Families_element):
                ordp = moments.ordp
                var_prec = min(var_prec, moments._var_prec)
                # Note: imposing var_prec precision ins't done in this method.
                # You should call normalize if you want that.
                moments = moments._moments[:parent.length_of_moments(par_p_prec)]
                moments = moments.change_ring(parent.base_ring())
            elif isinstance(moments, CoeffMod_OMS_element):
                # Coerce in using constant family
                # Need to truncate
                ordp = moments.ordp
                r_prec = moments.precision_relative()
                p_prec = min(parent.length_reverse_lookup(r_prec), par_p_prec)
                moments = parent.approx_module(p_prec, var_prec)(moments._moments[:parent.length_of_moments(p_prec)])
            elif hasattr(moments, '__len__'):
                #Need to modify if uniformiser is not p
                #Need to deal with var_prec
                if len(moments) == 0 or moments == 0:   #should do something with how "zero" moments is
                    V = parent.approx_module(0, var_prec) #var_prec?
                    moments = V([])
                    ordp = par_p_prec
                else:
                    R = parent.base_ring()
                    K = R.base_extend(R.base_ring().fraction_field())
                    p_prec = min(parent.length_reverse_lookup(len(moments)), par_p_prec)
                    #figure out var_prec from moments
                    V = parent.approx_module(p_prec, var_prec)
                    VK = V.base_extend(K)
                    moments = VK(moments[:parent.length_of_moments(p_prec)])
                    ordp = min([_padic_val_of_pow_series(a) for a in moments])
                    moments = V([_shift_coeffs(a, ordp) for a in moments])
            elif moments == 0:  #should do something with how "zero" moments is
                V = parent.approx_module(0, var_prec)
                moments = V([])
                ordp = par_p_prec
            else:
                R = parent.base_ring()
                Rbase = R.base_ring()
                K = R.base_extend(R.base_ring().fraction_field())
                V = parent.approx_module(1, var_prec)
                moments = K(moments)
                ordp = _padic_val_of_pow_series(moments)
                if ordp != 0:
                    moments = V([_shift_coeffs(moments, ordp), Rbase(0, 1)])
                else:
                    moments = V([moments, Rbase(0, 1)])
        self._moments = moments
        #if var_prec is None:
        #    self._var_prec = parent._prec_cap[1]
        #else:
        #    self._var_prec = var_prec
        self.ordp = ordp
        self._var_prec = var_prec
    
    #def __reduce__(self):
#        """
#        TESTS::
#        
#            sage: D = FamiliesOfOverconvergentDistributions(2, prec_cap=[8,5], base_coeffs=ZpCA(3))
#            sage: mu = D.random_element()
#            sage: loads(dumps(mu)) == mu
#            True
#        """
#        from sage.modular.pollack_stevens.coeffmod_OMS_families_element import create__CoeffMod_OMS_Families_element
#        return (create__CoeffMod_OMS_Families_element, (self._moments, self.parent(), self.ordp, self._var_prec))
    
    #def _relprec(self):
    #    return len(self._moments)
    def _repr_(self):
        #RH: adapted from coeffmod_OMS_element.py
        self.normalize()
        valstr = "("
        if self.ordp == 1:
            valstr = "%s * "%(self.parent().prime()) + valstr
        elif self.ordp != 0:
            valstr = "%s^%s * "%(self.parent().prime(), self.ordp) + valstr
        n = len(self._moments)
        if n == 0:
            return valstr + ") + O(w^%s)"%(self._var_prec)
        for i in range(n):
            mom_str = repr(self._moments[i])
            big_oh = mom_str.rfind(" + O(%s"%(self._moments[i].variable()))
            if big_oh > -1:
                valstr += mom_str[:big_oh] + ", "
            else:
                valstr += "0, "
        return valstr[:-2] + ") + O(w^%s)"%(self._var_prec)
    
    def character(self):
        """
        Returns the tame character of self.
        """
        return self._character
    
    def _add_(self, right):
        #RH: copied from coeffmod_OMS_element.py
        return self._addsub(right, False)
    
    def _sub_(self, right):
        #RH: copied from coeffmod_OMS_element.py
        return self._addsub(right, True)
    
    def _addsub(self, right, negate):
        #RH: adapted from coeffmod_OMS_element.py
        ans = self.parent()(0)
        saprec = self.precision_absolute()
        raprec = right.precision_absolute()
        aprec = min(saprec[0], raprec[0])
        ans.ordp = min(self.ordp, right.ordp)
        rprec = aprec - ans.ordp
        var_prec = min(saprec[1], raprec[1])
        V = ans.parent().approx_module(rprec, var_prec)
        R = V.base_ring()
        smoments = copy(self._moments)
        rmoments = copy(right._moments)
        length = self.parent().length_of_moments(rprec)
        if smoments.parent() is not V:
            smoments = V(smoments.list(copy=False)[:length] + ([R(0)] * (length - len(smoments)) if length > len(smoments) else []))
        if rmoments.parent() is not V:
            rmoments = V(rmoments.list(copy=False)[:length] + ([R(0)] * (length - len(rmoments)) if length > len(rmoments) else []))
        # We multiply by the relative power of p
        if self.ordp > right.ordp:
            for i in range(length):
                smoments[i] = _shift_coeffs(smoments[i], self.ordp - right.ordp, right=False)
        elif self.ordp < right.ordp:
            for i in range(length):
                rmoments[i] = _shift_coeffs(rmoments[i], right.ordp - self.ordp, right=False)
        if negate:
            rmoments = -rmoments
        ans._moments = smoments + rmoments
        return ans
    
    def _lmul_(self, right):
        #RH: adapted from coeffmod_OMS_element.py
        """
        Scalar multiplication self*right.
        """
        ans = CoeffMod_OMS_Families_element(None, self.parent(), None, False, var_prec=self._var_prec)
        if right.is_zero():
            ans._moments = self.parent().approx_module(0, self._var_prec)([])
            ans.ordp = min(self.parent().precision_cap()[0], right.valuation()+self.ordp)
        else:
            v, u = _padic_val_unit_of_pow_series(right)
            ans._moments = self._moments * u
            ans.ordp = self.ordp + v
        return ans
    
    def _rmul_(self, left):
        #RH: copied from coeffmod_OMS_element.py
        """
        Scalar multiplication left*self.
        """
        return self._lmul_(left)
    
    def __cmp__(left, right):
        r"""
        TESTS::
        
            sage: DD = FamiliesOfOverconvergentDistributions(4, prec_cap=[4,6], base_coeffs=ZpCA(7))
            sage: mus = DD(0)
            sage: mus == DD.zero()
            True
            sage: mus = DD(2)
            sage: nus = DD(1)
            sage: nus == mus
            False
            sage: mus == nus
            False
            sage: K = Qp(11, 6)
            sage: DG = DirichletGroup(11, K)
            sage: chi = DG([K(378703)]) 
            sage: DD = FamiliesOfOverconvergentDistributions(1, p=11, prec_cap=[6, 4], base_coeffs=ZpCA(11, 6), character=chi)
            sage: mus = DD(0)
            sage: mus == DD.zero()
            True
            sage: mus = DD(2)
            sage: nus = DD(1)
            sage: nus == mus
            False
            sage: mus == nus
            False
        """
        #RH: adapted from coeffmod_OMS_element.py
        #left = deepcopy(left)
        #right = deepcopy(right)
        left = CoeffMod_OMS_Families_element(left._moments, left.parent(), ordp=left.ordp, check=False, var_prec=left._var_prec)
        right = CoeffMod_OMS_Families_element(right._moments, right.parent(), ordp=right.ordp, check=False, var_prec=right._var_prec)
        left.normalize()
        right.normalize()
        lrprec = left.precision_relative()[0]
        rrprec = right.precision_relative()[0]
        if lrprec == 0:
            if rrprec == 0:
                return 0
            else:
                return -1
        elif rrprec == 0:
            return 1
        c = cmp(left.ordp, right.ordp)
        if c: return c
        length = left.parent().length_of_moments(min(lrprec, rrprec))
        #Will this work? Not sure power series compare well (see e.g. trac 9457)
        #So, ues a hack
        #return cmp(left._moments[:length], right._moments[:length])
        for i in range(length):
            c = cmp(left._moments[i].padded_list(), right._moments[i].padded_list())
            if c: return c
        return 0
#
#        p = left.parent().prime()
#        if left.ordp > right.ordp:
#            shift = p ** (left.ordp - right.ordp)
#            for i in range(rprec):
#                L = (shift * left._unscaled_moment(i)).padded_list()
#                R = (right._unscaled_moment(i)).padded_list()
#                c = cmp(L, R)
#                if c: return c
#        elif left.ordp < right.ordp:
#            shift = p ** (right.ordp - left.ordp)
#            for i in range(rprec):
#                L = (left._unscaled_moment(i)).padded_list()
#                R = (shift * right._unscaled_moment(i)).padded_list()
#                c = cmp(L, R)
#                if c: return c
#        else:
#            for i in range(rprec):
#                L = (left._unscaled_moment(i)).padded_list()
#                R = (right._unscaled_moment(i)).padded_list()
#                c = cmp(L, R)
#                if c: return c
#        return 0

#        parent_prec = self.parent()._prec_cap
#        var_prec = min(self._var_prec, parent_prec[1], other._var_prec)
#        p_prec = min(len(self._moments), parent_prec[0], len(other._moments))
#        for i in range(p_prec):
#            selflist = _add_big_ohs_list(self._moments[i], [((p_prec - i) * self._cp).ceil(), var_prec])
#            otherlist = _add_big_ohs_list(other._moments[i], [((p_prec - i) * self._cp).ceil(), var_prec])
#            c = cmp(selflist, otherlist)
#            if c:
#                return c
#        return 0
    
#    def __nonzero__(self):
#        """
#        Checks that self is non-zero up to precision ...
#        """
#        parent_prec = self.parent()._prec_cap
#        var_prec = min(self._var_prec, parent_prec[1])
#        p_prec = min(len(self._moments), parent_prec[0])
#        for i in range(p_prec):
#            selflist = _add_big_ohs_list(self._moments[i], [((p_prec - i) * self._cp).ceil(), var_prec])
#            for c in selflist:
#                if c != 0:
#                    return True
#        return False
    
    def is_zero(self, prec=None):
        self.normalize()
        n, v_prec = self.precision_relative()
        if n == 0:
            return True
        else:
            return False
        aprec, v_aprec = self.precision_absolute()
        if prec is None:
            prec = [n, v_prec]
        elif not hasattr(prec, '__len__'):
            prec = [ZZ(prec), v_prec]
        elif prec[0] > aprec or prec[1] > v_aprec:
            return False    #Should this raise a PrecisionError instead
        p_precs = self.parent().filtration_precisions(prec[0])
        for a in xrange(len(p_precs)):
            if not self._unscaled_moment(a)._is_zero_padic_power_series([p_precs[a], prec[1]]):
                return False
        return True
        
#        if prec is None:
#            prec = s_prec_rel
#        elif len(prec) > 2:
#            raise TypeError("prec must have length at most 2.")
#        if prec[0] > s_prec_rel[0] or prec[1] > s_prec_rel[1]:
#            return False
#        num_moments = min(prec[0], s_prec_rel[0]) if prec[0] is not None else s_prec_rel[0]
#        var_prec = min(prec[1], s_prec_rel[1]) if prec[1] is not None else s_prec_rel[1]
#        p_precs = self.parent().filtration_precisions(num_moments)
#        for i in range(num_moments):
#            selflist = _add_big_ohs_list(self._moments[i], [p_precs[i], var_prec])
#            for c in selflist:
#                if c != 0:
#                    return False
#        return True
    
    def find_scalar(self, other, M = None, check=True):
        r"""
            M[0] is the desired p-adic precision.
        """
        #self.normalize()
        #other.normalize()
        n, s_var_prec = self.precision_relative()
        other_pr, other_var_prec = other.precision_relative()
        if n == 0:
            raise ValueError("self is zero")
        if M is None:
            M = [n, s_var_prec]
        elif not hasattr(moments, '__len__'):
            M = [ZZ(M), ZZ(min(s_var_prec, other_var_prec))]
        elif len(M) > 2:
            raise TypeError("prec must have length at most 2.")
        #elif not isinstance(M, (list, tuple)):
        #    M = [M, None]
        #if min(s_var_prec, other_var_prec) < M[1]:
        #    raise ValueError("Insufficient precision in variable (%s requested, min(%s, %s) obtained)"%(M[1], s_var_prec, other_var_prec))
        i = 0
        verbose("n = %s"%n)
        verbose("moment 0")
        a = self._unscaled_moment(i)
        verbose("a, val = %s, %s"%(a, self.ordp))
        p = self.parent().prime()
        v = _padic_val_of_pow_series(a, p)
        R = self.parent().base_ring()
        p_precs = self.parent().filtration_precisions(n)
        while v >= p_precs[i]:
            i += 1
            verbose("p moment %s"%i)
            try:
                a = self._unscaled_moment(i)
            except IndexError:
                raise ValueError("self is zero")
            v = _padic_val_of_pow_series(a, p)
            verbose("a, val = %s, %s"%(a, self.ordp))
        relprec = p_precs[i] - v
        var_prec = min(s_var_prec, other_var_prec)  #should this depend on w-adic valuation?
        Rbase = R.base_ring()
        RK = R.change_ring(Rbase.fraction_field())
        other_length = self.parent().length_of_moments(other_pr)
        if i < other_length:
            #Hack
            verbose("val, other._unscaled_moment(%s) = %s, %s"%(i, other.ordp, other._moments[i]))
            #alpha = _sanitize_alpha(RK(other._unscaled_moment(i)) / RK(a))
            alpha = RK(other._unscaled_moment(i)) / RK(a)
            alpha = _truncate_low_prec(alpha)
            #alpha, var_prec_loss = _custom_ps_div(other._unscaled_moment(i), a)
            #alpha = R(_add_big_ohs_list(alpha, [ceil((n - i) * self._cp), min(s_var_prec, other_var_prec)]))
            #alpha = R(_add_big_ohs_list(other._unscaled_moment(i) / a, [ceil((n - i) * self._cp), min(s_var_prec, other_var_prec)]))
        else:
            #Fix var_prec??
            alpha = RK([Rbase(0, p_precs[i])], var_prec)
        verbose("alpha = %s"%(alpha))
        
        while i < other_length-1:
            i += 1
            verbose("comparing p moment %s"%i)
            a = self._unscaled_moment(i)
            if check:
#                   verbose("self.moment=%s, other.moment=%s"%(a, other._unscaled_moment(i)))
                #if other._unscaled_moment(i) != _add_big_ohs_list(alpha * a, [ceil((n - i) * self._cp), var_prec]):
                verbose("val, self._unscaled_moment(%s) = %s, %s"%(i, self.ordp, self._moments[i]))
                verbose("val, other._unscaled_moment(%s) = %s, %s"%(i, other.ordp, other._moments[i]))
                if other._unscaled_moment(i) != alpha * a:
                    var_prec = (other._unscaled_moment(i) - alpha * a).valuation()
                    if var_prec < 1:
                        print i,"\n"
                        print other._unscaled_moment(i),"\n"
                        print alpha * a,"\n"
                        print other._unscaled_moment(i) - alpha * a,"\n"
                        print var_prec,"\n"
                        raise ValueError("not a scalar multiple")
                    alpha = alpha.add_bigoh(var_prec)
                    verbose("Resetting alpha (lowering variable prec)")
            v = _padic_val_of_pow_series(a, p)
            if p_precs[i] - v > relprec:
                verbose("Resetting alpha: relprec=%s, i=%s, p_precs[i]=%s, v=%s"%(relprec, i,p_precs[i], v))
                relprec = p_precs[i] - v
                #Should we alter var_prec, too?
                #Hack
                alpha = RK(other._unscaled_moment(i)) / RK(a)
                alpha = _truncate_low_prec(alpha)   #Is this necessary?
                #alpha = R(_add_big_ohs_list(other._unscaled_moment(i) / a, [ceil((n - i) * self._cp), min(s_var_prec, other_var_prec)]))
                verbose("alpha=%s"%(alpha))
        if relprec < M[0]:  #May need to change this
            raise ValueError("result not determined to high enough precision")
        #if var_prec < M[1]:
        
        alpha = alpha * self.parent().prime()**(other.ordp - self.ordp)
        verbose("alpha=%s"%(alpha))
        #Hack to workaround bug with power series
        #RK = PowerSeriesRing(R.base_ring().fraction_field(), R.variable_name(), default_prec=R.default_prec()) #Hack because of bugs in power series
        #alpha_list = _add_big_ohs_list(alpha, [relprec, M[1]])
        #try:
        #    return R(alpha)    #if alpha = O(p^7) + (... + O(p^7)*w + ..., then R(alpha) = O(p^8) + w * ()... if say 8 is prec of Rbase. BUG IN SAGE!
        #    #return RK(alpha_list) #Do we need R(*) here?
        #except ValueError:
        #    return alpha
        return alpha
    
    def precision_relative(self):
        return [self.parent().length_reverse_lookup(len(self._moments)), self._var_prec]
    
    def precision_absolute(self):
        return [ZZ(self.parent().length_reverse_lookup(len(self._moments)) + self.ordp), self._var_prec]
    
    def valuation(self):
        r"""
        Returns the `\varpi`-adic valuation of this family of distributions.
        
        .. WARNING::

            This function modifies the distribution in place since it calls
            :meth:`~sage.modular.pollack_stevens.coeffmod_OMS_families_element.CoeffMod_OMS_Families_element.normalize`.
        """
        return self.normalize().ordp
    
    def _valuation(self, val_vector=False):
        n, v_prec = self.precision_relative()
        if n == 0:
            if val_vector:
                return [self.ordp, []]
            return self.ordp
        length = self.parent().length_of_moments(n)
        p_precs = self.parent().filtration_precisions(n)
        #print "ordp, moms:", self.ordp, "\n", self._moments
        if val_vector:
            cur = _padic_val_of_pow_series(self._unscaled_moment(0), var_prec=self._var_prec)
            min_val = cur# if not cur.is_zero() else Infinity
            vv = [min_val]
            for a in range(1, length):
                cur_mom = self._unscaled_moment(a)
                cur = _padic_val_of_pow_series(cur_mom, var_prec=self._var_prec) if not cur_mom.is_zero() else n - p_precs[a] + _padic_val_of_pow_series(cur_mom, var_prec=self._var_prec)
                if cur < min_val:
                    min_val = cur
                vv.append(min_val)
            ret = self.ordp + min_val#min(n, min_val)
            if ret == Infinity: #This can happen because power series over ZpCA don't treat the zero power series properly
                ret = self.ordp + n
            #verbose("ret %s"%(ret), level=2)
            #verbose("\n******** end valuation ********", level=2)
            return [ret, vv]
        ret = self.ordp + min([_padic_val_of_pow_series(self._unscaled_moment(a), var_prec=self._var_prec) if not self._unscaled_moment(a).is_zero() else n - p_precs[a] + _padic_val_of_pow_series(self._unscaled_moment(a), var_prec=self._var_prec) for a in range(length)])
        if ret == Infinity: #This can happen because power series over ZpCA don't treat the zero power series properly
            ret = self.ordp + n
        return ret
    
    def normalize(self):
        r"""
        TESTS::
        
            sage: from sage.modular.pollack_stevens.coeffmod_OMS_families_element import CoeffMod_OMS_Families_element
            sage: DD = FamiliesOfOverconvergentDistributions(0, base_coeffs=ZpCA(5, 4), prec_cap=[4,2])
            sage: V = DD.approx_module()
            sage: R = DD.base_ring()
            sage: w = R.gen()
            sage: mus = V((4*5 + 5^3 + O(5^4) + (1 + 3*5 + 5^2 + 4*5^3 + O(5^4))*w + O(w^2), 1 + 5 + 4*5^2 + 5^3 + O(5^4) + (5 + 3*5^2 + 5^3 + O(5^4))*w + O(w^2), 4 + 2*5 + 2*5^2 + 4*5^3 + O(5^4) + (1 + 2*5^2 + 5^3 + O(5^4))*w + O(w^2), 1 + 5 + O(5^4) + (3 + 5 + 4*5^2 + O(5^4))*w + O(w^2), 4 + 3*5 + 5^2 + 2*5^3 + O(5^4) + (1 + 3*5 + 4*5^2 + O(5^4))*w + O(w^2), 2*5 + 4*5^2 + 2*5^3 + O(5^4) + (2*5 + 3*5^2 + 5^3 + O(5^4))*w + O(w^2)))
            sage: mus = CoeffMod_OMS_Families_element(mus, DD, check=False, var_prec=2)
            sage: mus   #Indirect doctest
            (4*5 + 5^3 + O(5^4) + (1 + 3*5 + 5^2 + 4*5^3 + O(5^4))*w, 1 + 5 + 4*5^2 + 5^3 + O(5^4) + (5 + 3*5^2 + 5^3 + O(5^4))*w, 4 + 2*5 + 2*5^2 + O(5^3) + (1 + 2*5^2 + O(5^3))*w, 1 + 5 + O(5^2) + (3 + 5 + O(5^2))*w, 4 + O(5) + (1 + O(5))*w, 0) + O(w^2)
            sage: mus = V((4*5 + 5^3 + O(5^4) + (1 + 3*5 + 5^2 + 4*5^3 + O(5^4))*w + O(w^2), 1 + 5 + O(5^2) + (5 + 3*5^2 + 5^3 + O(5^4))*w + O(w^2), 4 + 2*5 + 2*5^2 + 4*5^3 + O(5^4) + (1 + 2*5^2 + 5^3 + O(5^4))*w + O(w^2), 1 + 5 + O(5^4) + (3 + 5 + 4*5^2 + O(5^4))*w + O(w^2), 4 + 3*5 + 5^2 + 2*5^3 + O(5^4) + (1 + 3*5 + 4*5^2 + O(5^4))*w + O(w^2), 2*5 + 4*5^2 + 2*5^3 + O(5^4) + (2*5 + 3*5^2 + 5^3 + O(5^4))*w + O(w^2)))
            sage: mus = CoeffMod_OMS_Families_element(mus, DD, check=False, var_prec=2); mus   #Indirect doctest
            (4*5 + O(5^2) + (1 + 3*5 + O(5^2))*w, 1 + 5 + O(5^2) + (5 + O(5^2))*w, 4 + O(5) + (1 + O(5))*w) + O(w^2)
            sage: mus = V((4*5^2 + 5^3 + O(5^4) + (5^2 + 4*5^3 + O(5^4))*w + O(w^2), 4*5^2 + 5^3 + O(5^42) + (3*5^2 + 5^3 + O(5^4))*w + O(w^2), 2*5^2 + 4*5^3 + O(5^4) + (2*5^2 + 5^3 + O(5^4))*w + O(w^2), O(5^4) + (4*5^2 + O(5^4))*w + O(w^2), 5^2 + 2*5^3 + O(5^4) + (4*5^2 + O(5^4))*w + O(w^2), 4*5^2 + 2*5^3 + O(5^4) + (3*5^2 + 5^3 + O(5^4))*w + O(w^2)))
            sage: mus = CoeffMod_OMS_Families_element(mus, DD, check=False, var_prec=2); mus    #Indirect doctest
            5^2 * (4 + 5 + O(5^2) + (1 + 4*5 + O(5^2))*w, 4 + 5 + O(5^2) + (3 + 5 + O(5^2))*w, 2 + O(5) + (2 + O(5))*w) + O(w^2)
            sage: mus = V((4*5^2 + 5^3 + O(5^4) + (5^2 + 4*5^3 + O(5^4))*w + O(w^2), 4*5^2 + 5^3 + O(5^42) + (3*5^2 + 5^3 + O(5^4))*w + O(w^2), 2*5^2 + 4*5^3 + O(5^4) + (2*5^2 + 5^3 + O(5^4))*w + O(w^2), O(5^4) + (4*5^2 + O(5^4))*w + O(w^2), 5^2 + 2*5^3 + O(5^4) + (4*5^2 + O(5^4))*w + O(w^2), 4*5^2 + 2*5^3 + O(5^4) + (2*5 + 3*5^2 + 5^3 + O(5^4))*w + O(w^2)))
            sage: mus = CoeffMod_OMS_Families_element(mus, DD, check=False, var_prec=2); mus   #Indirect doctest
            5^2 * (4 + 5 + O(5^2) + (1 + 4*5 + O(5^2))*w, 4 + 5 + O(5^2) + (3 + 5 + O(5^2))*w, 2 + O(5) + (2 + O(5))*w) + O(w^2)
        """
        n, v_prec = self.precision_relative()
        if n == 0:
            return self
        adjust_moms = 0
        p_precs = self.parent().filtration_precisions(n)
        length = len(p_precs)
        for i in range(length):
            adjust_moms = max(adjust_moms, p_precs[i] - _padic_abs_prec_of_pow_series(self._moments[i], v_prec))
        #print "adjust_moms:", adjust_moms
        if adjust_moms >= n:
            assert False    #Deal with this later...
        R = self.parent().base_ring()
        if adjust_moms == 0:
            for i in range(length):
                self._moments[i] = R(_add_big_ohs_list(self._moments[i], [p_precs[i], self._var_prec]), self._var_prec)
        else:
            n -= adjust_moms
            p_precs = self.parent().filtration_precisions(n)
            length = len(p_precs)
            V = self.parent().approx_module(n, self._var_prec)
            self._moments = V([R(_add_big_ohs_list(self._moments[i], [p_precs[i], self._var_prec]), self._var_prec) for i in range(length)])
        val_diff = self._valuation() - self.ordp
        #print "val_diff:", val_diff
        if val_diff > 0:
            n -= val_diff
            p_precs = self.parent().filtration_precisions(n)
            length = len(p_precs)
            self.ordp += val_diff
            V = self.parent().approx_module(n, self._var_prec)
            self._moments = V([_shift_coeffs(self._moments[i], val_diff) for i in range(length)])
        return self
    
    def normalize_old(self):
        #RH: adapted from coeffmod_OMS_element.py
        #Not tested
        V = self._moments.parent()
        R = V.base_ring()
        n, v_prec = self.precision_relative()
        if n == 0:
            return self
        self_aprec = self.precision_absolute()[0] 
        self_ordp = self.ordp
        self_val, val_vector = self._valuation(val_vector=True)
        while True:
            ## RP -- just added these lines to stop crashing, but doesn't quite work -- takes more than once to fully normalize to 0.
            if self_val == Infinity:
                V = self.parent().approx_module(0, self._var_prec)
                self._moments = V([])
                self.ordp = self_aprec
                return self
            shift = self_val - self.ordp
            self.ordp = self_val
            #Factor out powers of uniformizer and check precision
            adjust_moms = 0
            p_precs = self.parent().filtration_precisions(n)
            length = len(p_precs)
            verbose("n: %s; shift: %s; _mom: %s"%(n, shift, self._moments), level=2)
            if shift > 0:
                for i in range(length):
                    self._moments[i] = _shift_coeffs(self._moments[i], shift)
                    adjust_moms = max(adjust_moms, p_precs[i] - _padic_abs_prec_of_pow_series(self._moments[i], v_prec))
            elif shift == 0:
                for i in range(length):
                    adjust_moms = max(adjust_moms, p_precs[i] - _padic_abs_prec_of_pow_series(self._moments[i], v_prec))
            else:
                raise NotImplementedError("Currently only deals with the case where the base ring is a ring of integers.")
            #Cut down moments because of precision loss
            verbose("adjust_moms: %s\nn %s\n_moms: %s"%(adjust_moms, n, self._moments), level=2)
            if adjust_moms >= n:
                V = self.parent().approx_module(0, self._var_prec)
                self._moments = V([])
                #self.ordp = adjust_moms    #should we take min with parent().precision_cap()?
                verbose("adjust_moms %s, n %s, self.ordp %s"%(adjust_moms, n, self.ordp)) 
                verbose("output: ordp %s, _moments %s"%(self.ordp, self._moments), level=2)
                verbose("\n******** end normalize ********", level=2)
                return self
            if adjust_moms == 0:
                for i in range(length):
                    self._moments[i] = R(_add_big_ohs_list(self._moments[i], [p_precs[i], v_prec]), v_prec)
                verbose("adjust_moms %s, n %s, self.ordp %s"%(adjust_moms, n, self.ordp)) 
                verbose("output: ordp %s, _moments %s"%(self.ordp, self._moments), level=2)
                verbose("\n******** end normalize ********", level=2)
                return self
            n -= adjust_moms
            p_precs = self.parent().filtration_precisions(n)
            length = len(p_precs)
            val_diff = val_vector[length-1] - shift
            if val_diff == 0:
                #p_precs = self.parent().filtration_precisions(n)
                V = self.parent().approx_module(n, self._var_prec)
                self._moments = V([R(_add_big_ohs_list(self._moments[i], [p_precs[i], v_prec]), v_prec) for i in range(length)])
                verbose("output: ordp %s, _moments %s"%(self.ordp, self._moments), level=2)
                verbose("\n******** end normalize ********", level=2)
                return self
            self_val = val_vector[length-1] + self_ordp
            val_vector = [val_vector[i] - shift for i in range(length)]
            self_ordp = self.ordp
            verbose("\nn: %s, self_val: %s, val_vector: %s, self.ordp: %s, self_ordp: %s"%(n, self_val, val_vector, self.ordp, self_ordp))
        
#        p = self.parent().prime()
#        self_val = self._valuation()
#        shift = self_val - self.ordp
#        self.ordp = self_val
#        #Factor out powers of uniformizer and check precision
#        p_precs = self.parent().filtration_precisions(n)
#        adjust_moms = 0
#        verbose("n: %s; shift: %s; _mom: %s\np_precs: %s"%(n, shift, self._moments, p_precs), level=2)
#        if shift > 0:
#            for i in range(n):
#                self._moments[i] = _shift_coeffs(self._moments[i], shift)
#                adjust_moms = max(adjust_moms, p_precs[i] - _padic_abs_prec_of_pow_series(self._moments[i], v_prec))
#        elif shift == 0:
#            for i in range(n):
#                adjust_moms = max(adjust_moms, p_precs[i] - _padic_abs_prec_of_pow_series(self._moments[i], v_prec))
#        else:
#            raise NotImplementedError("Currently only deals with the case where the base coefficients are a ring of integers.")
#        #Cut down moments because of precision loss
#        verbose("adjust_mom: %s\nn %s \n_moms: %s"%(adjust_moms, n, self._moments), level=2)
#        if adjust_moms >=n:
#            V = self.parent().approx_module(0)
#            self._moments = V([])
#            #self.ordp = adjust_moms
#            verbose("adjust_moms %s, \nn %s, \nself.ordp %s"%(adjust_moms, n, self.ordp))
#        elif adjust_moms > 0:
#            n -= adjust_moms    #Is this going to give the correct precision?
#            p_precs = self.parent().filtration_precisions(n)
#            verbose("new p_precs: %s"%(p_precs), level=2)
#            R = self.parent().base_ring()
#            for i in range(n):
#                self._moments[i] = R(_add_big_ohs_list(self._moments[i], [p_precs[i], v_prec]), v_prec)
#        else:
#            R = self.parent().base_ring()
#            for i in range(n):
#                self._moments[i] = R(_add_big_ohs_list(self._moments[i], [p_precs[i], v_prec]), v_prec)
#        return self

    def moment(self, n):
        #RH: "copied" from dist.pyx
        r"""
        Returns the `n`-th moment.
        """
        if self.ordp == 0:
            return self._unscaled_moment(n)
        else:
            return _shift_coeffs(self._unscaled_moment(n), self.ordp, right=False)
    
    def _unscaled_moment(self, n):
        #RH: "copied" from coeffmod_OMS_element.py
        r"""
        Returns the `n`-th moment, unscaled by the overall power of p stored in self.ordp.
        """
        return self._moments[n]
    
    def reduce_precision(self, new_prec):
        rprec = self.precision_relative()
        try:
            p_prec, var_prec = new_prec
        except TypeError:   #if new_prec isn't iterable (wrong size is ValueError)
            p_prec = new_prec
            var_prec = rprec[1]
        if p_prec > rprec[0] or var_prec > rprec[1]:
            raise ValueError("Precisions specified must be less than current precisions.")
        moments = self._moments[:self.parent().length_of_moments(p_prec)]
        ordp = self.ordp
        return CoeffMod_OMS_Families_element(moments, self.parent(), ordp=ordp, check=False, var_prec=var_prec)
    
    def reduce_precision_absolute(self, new_prec):
        aprec = self.precision_absolute()
        try:
            p_prec, var_prec = new_prec
        except TypeError:   #if new_prec isn't iterable (wrong size is ValueError)
            p_prec = new_prec
            var_prec = aprec[1]
        if p_prec > aprec[0] or var_prec > aprec[1]:
            raise ValueError("Precisions specified must be less than current precisions.")
        if p_prec == aprec[0] and var_prec == aprec[1]:
            return self
        ordp = self.ordp
        if p_prec - ordp <= 0:
            moments = self.parent().approx_module(0)([])
            ordp = new_prec
        else:
            moments = self._moments[:self.parent().length_of_moments(p_prec - ordp)]
        return CoeffMod_OMS_Families_element(moments, self.parent(), ordp, check=False, var_prec=var_prec)
    
    def solve_diff_eqn(self):
        r"""
        EXAMPLES::
        
            sage: from sage.modular.pollack_stevens.coeffmod_OMS_families_element import CoeffMod_OMS_Families_element
            sage: DD = FamiliesOfOverconvergentDistributions(0, base_coeffs=ZpCA(11, 8), prec_cap=[8, 8])
            sage: V = DD.approx_module()
            sage: R = DD.base_ring()
            sage: w = R.gen()
            sage: mu = V((R(0,8), R(5 + 2*11 + 8*11^2 + 2*11^3 + 10*11^4 + 5*11^5 + 10*11^6 + O(11^8) + (9 + 7*11 + 3*11^2 + 3*11^3 + 5*11^4 + 11^5 + 6*11^6 + O(11^8))*w + (2 + 5*11^2 + 5*11^4 + 10*11^5 + 6*11^6 + 5*11^7 + O(11^8))*w^2 + (2 + 6*11 + 10*11^2 + 10*11^3 + 7*11^4 + 10*11^5 + 9*11^6 + 8*11^7 + O(11^8))*w^3 + (6*11^2 + 9*11^3 + 4*11^5 + 8*11^6 + 3*11^7 + O(11^8))*w^4 + (6 + 3*11 + 8*11^2 + 3*11^3 + 11^4 + 10*11^5 + 2*11^6 + 2*11^7 + O(11^8))*w^5 + (3*11 + 7*11^2 + 4*11^3 + 7*11^4 + 10*11^5 + 5*11^7 + O(11^8))*w^6 + (7 + 4*11 + 6*11^2 + 9*11^4 + 10*11^5 + 5*11^6 + 4*11^7 + O(11^8))*w^7 + O(w^8)), R(1 + 10*11 + 10*11^2 + 6*11^3 + 8*11^4 + 6*11^5 + 2*11^6 + O(11^7) + (2 + 7*11 + 10*11^2 + 5*11^3 + 4*11^5 + 3*11^6 + O(11^7))*w + (9 + 5*11 + 6*11^2 + 7*11^3 + 8*11^5 + 7*11^6 + O(11^7))*w^2 + (5 + 2*11 + 7*11^2 + 10*11^3 + 4*11^4 + 5*11^5 + 3*11^6 + O(11^7))*w^3 + (10 + 5*11 + 10*11^2 + 9*11^3 + 10*11^4 + 5*11^5 + 3*11^6 + O(11^7))*w^4 + (1 + 10*11 + 3*11^2 + 3*11^3 + 8*11^4 + 7*11^5 + 8*11^6 + O(11^7))*w^5 + (9 + 10*11^2 + 10*11^3 + 3*11^4 + 5*11^5 + 5*11^6 + O(11^7))*w^6 + (3 + 6*11 + 11^2 + 10*11^3 + 10*11^4 + 3*11^6 + O(11^7))*w^7 + O(w^8)), R(5 + 10*11 + 11^2 + 11^3 + 4*11^4 + 11^5 + O(11^6) + (3 + 5*11 + 8*11^2 + 4*11^3 + 6*11^4 + 5*11^5 + O(11^6))*w + (10 + 2*11 + 8*11^2 + 3*11^4 + 9*11^5 + O(11^6))*w^2 + (10 + 7*11 + 8*11^3 + 3*11^4 + 2*11^5 + O(11^6))*w^3 + (3 + 5*11 + 8*11^2 + 7*11^3 + 7*11^4 + O(11^6))*w^4 + (6 + 3*11 + 9*11^3 + 8*11^4 + 7*11^5 + O(11^6))*w^5 + (10 + 9*11 + 5*11^2 + 2*11^3 + 6*11^4 + 11^5 + O(11^6))*w^6 + (7 + 8*11 + 9*11^2 + 3*11^3 + 7*11^4 + 8*11^5 + O(11^6))*w^7 + O(w^8)), R(9 + 9*11^2 + 2*11^3 + 4*11^4 + O(11^5) + (9*11 + 10*11^2 + 4*11^3 + 7*11^4 + O(11^5))*w + (1 + 8*11 + 2*11^2 + 6*11^3 + 6*11^4 + O(11^5))*w^2 + (2 + 3*11^2 + 3*11^3 + 11^4 + O(11^5))*w^3 + (2 + 2*11 + 3*11^2 + 8*11^3 + 4*11^4 + O(11^5))*w^4 + (8*11 + 5*11^2 + 8*11^3 + 3*11^4 + O(11^5))*w^5 + (9 + 6*11 + 6*11^2 + 4*11^3 + 10*11^4 + O(11^5))*w^6 + (3 + 5*11 + 11^3 + 4*11^4 + O(11^5))*w^7 + O(w^8)), R(1 + 3*11 + 4*11^2 + 2*11^3 + O(11^4) + (10*11 + 7*11^2 + 9*11^3 + O(11^4))*w + (3 + 11 + 7*11^2 + 4*11^3 + O(11^4))*w^2 + (9 + 9*11 + 5*11^3 + O(11^4))*w^3 + (4 + 7*11 + 6*11^2 + O(11^4))*w^4 + (7 + 9*11 + 3*11^2 + 5*11^3 + O(11^4))*w^5 + (3 + 8*11 + 11^2 + 11^3 + O(11^4))*w^6 + (9 + 3*11 + 9*11^2 + 3*11^3 + O(11^4))*w^7 + O(w^8)), R(3 + 7*11 + 10*11^2 + O(11^3) + (9 + 2*11 + 10*11^2 + O(11^3))*w + (7 + 11 + 10*11^2 + O(11^3))*w^2 + (4 + 4*11 + 2*11^2 + O(11^3))*w^3 + (7 + 9*11 + 6*11^2 + O(11^3))*w^4 + (4 + 8*11 + 3*11^2 + O(11^3))*w^5 + (3 + 11 + 2*11^2 + O(11^3))*w^6 + (9 + 4*11 + 8*11^2 + O(11^3))*w^7 + O(w^8)), R(4 + 3*11 + O(11^2) + (7 + 8*11 + O(11^2))*w + (6 + 4*11 + O(11^2))*w^2 + (4 + 3*11 + O(11^2))*w^3 + (5 + O(11^2))*w^4 + (4*11 + O(11^2))*w^5 + (3*11 + O(11^2))*w^6 + (8 + 9*11 + O(11^2))*w^7 + O(w^8)), R(10 + O(11) + (9 + O(11))*w + (3 + O(11))*w^2 + (7 + O(11))*w^3 + (1 + O(11))*w^4 + (3 + O(11))*w^5 + (7 + O(11))*w^6 + (9 + O(11))*w^7 + O(w^8))))
            sage: mu = CoeffMod_OMS_Families_element(mu, DD, check=False, var_prec=8)
            sage: nu = mu.solve_diff_eqn()
            sage: S0 = DD.action().actor()
            sage: Delta_mat = S0([1,1,0,1])
            sage: nu * Delta_mat - nu - mu
            11^7 * () + O(w^8)
            sage: DD = FamiliesOfOverconvergentDistributions(0, base_coeffs=ZpCA(5, 8), prec_cap=[8, 8])
            sage: V = DD.approx_module()
            sage: R = DD.base_ring()
            sage: w = R.gen()
            sage: mu = V((O(w^8), 4 + 5 + 5^2 + 4*5^3 + 2*5^4 + O(5^8) + (3 + 2*5 + 2*5^2 + 2*5^4 + 4*5^5 + 4*5^6 + 5^7 + O(5^8))*w + (3*5 + 4*5^2 + 4*5^4 + 2*5^5 + 5^6 + 4*5^7 + O(5^8))*w^2 + (2 + 3*5^2 + 2*5^3 + 5^4 + 2*5^6 + 4*5^7 + O(5^8))*w^3 + (1 + 4*5 + 4*5^2 + 3*5^3 + 4*5^4 + 2*5^5 + 2*5^6 + 5^7 + O(5^8))*w^4 + (3 + 2*5^3 + 5^4 + 4*5^5 + 2*5^6 + 5^7 + O(5^8))*w^5 + (1 + 2*5 + 5^2 + 3*5^3 + 4*5^4 + 5^5 + 3*5^7 + O(5^8))*w^6 + (3 + 2*5 + 3*5^3 + 3*5^4 + O(5^8))*w^7 + O(w^8), 2 + 5 + 2*5^2 + 4*5^3 + 5^5 + 5^6 + O(5^8) + (4 + 3*5^3 + 2*5^4 + 4*5^5 + 4*5^6 + O(5^8))*w + (1 + 3*5 + 3*5^2 + 4*5^3 + 2*5^4 + 5^6 + 5^7 + O(5^8))*w^2 + (4 + 3*5 + 3*5^3 + 2*5^5 + 2*5^6 + 5^7 + O(5^8))*w^3 + (4 + 5 + 2*5^2 + 3*5^3 + 4*5^4 + O(5^8))*w^4 + (1 + 3*5 + 4*5^3 + 2*5^4 + 2*5^6 + 3*5^7 + O(5^8))*w^5 + (3 + 2*5 + 5^2 + 2*5^3 + 3*5^4 + 4*5^5 + 5^6 + 4*5^7 + O(5^8))*w^6 + (3 + 4*5 + 5^2 + 3*5^3 + 5^4 + 4*5^5 + 4*5^6 + 4*5^7 + O(5^8))*w^7 + O(w^8), 3 + 3*5 + 5^2 + 4*5^3 + 4*5^4 + 4*5^7 + O(5^8) + (4 + 2*5 + 2*5^2 + 5^3 + 4*5^4 + 3*5^5 + 3*5^6 + 4*5^7 + O(5^8))*w + (5^2 + 3*5^3 + 4*5^4 + 5^5 + 5^6 + 2*5^7 + O(5^8))*w^2 + (3 + 4*5 + 2*5^3 + 2*5^4 + 5^5 + 2*5^6 + 2*5^7 + O(5^8))*w^3 + (4*5 + 5^2 + 4*5^3 + 4*5^4 + 4*5^6 + 3*5^7 + O(5^8))*w^4 + (4 + 4*5 + 3*5^2 + 5^4 + 4*5^5 + 4*5^6 + 5^7 + O(5^8))*w^5 + (2 + 4*5 + 4*5^2 + 4*5^3 + 2*5^5 + 3*5^6 + 5^7 + O(5^8))*w^6 + (2 + 2*5^2 + 3*5^3 + 4*5^5 + 3*5^6 + 3*5^7 + O(5^8))*w^7 + O(w^8), 2 + 2*5 + 5^2 + 2*5^3 + 5^4 + 2*5^5 + 2*5^7 + O(5^8) + (4 + 3*5^2 + 3*5^3 + 3*5^4 + 3*5^6 + O(5^8))*w + (4*5 + 4*5^3 + 4*5^4 + 5^5 + 4*5^7 + O(5^8))*w^2 + (2 + 3*5 + 2*5^4 + 3*5^5 + 5^6 + O(5^8))*w^3 + (4*5 + 2*5^2 + 5^3 + 2*5^4 + 5^5 + 4*5^7 + O(5^8))*w^4 + (3 + 5 + 5^3 + 4*5^4 + 5^6 + 5^7 + O(5^8))*w^5 + (3*5 + 3*5^3 + 5^5 + 2*5^6 + 3*5^7 + O(5^8))*w^6 + (4*5 + 2*5^2 + 5^4 + 5^5 + 3*5^6 + 5^7 + O(5^8))*w^7 + O(w^8), 2 + 2*5 + 3*5^2 + 3*5^3 + 4*5^5 + 5^6 + 3*5^7 + O(5^8) + (2 + 4*5 + 2*5^4 + 5^5 + 4*5^6 + 5^7 + O(5^8))*w + (3 + 5 + 2*5^2 + 3*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 3*5^7 + O(5^8))*w^2 + (2 + 4*5 + 2*5^2 + 5^3 + 3*5^4 + 5^5 + 4*5^7 + O(5^8))*w^3 + (2 + 2*5 + 3*5^2 + 2*5^3 + 5^4 + 4*5^5 + 5^6 + 5^7 + O(5^8))*w^4 + (1 + 3*5 + 5^3 + 5^4 + 2*5^5 + 5^6 + 3*5^7 + O(5^8))*w^5 + (1 + 5 + 4*5^3 + 3*5^4 + 3*5^5 + 2*5^6 + 4*5^7 + O(5^8))*w^6 + (2*5^2 + 2*5^3 + 2*5^4 + 4*5^5 + 5^6 + 2*5^7 + O(5^8))*w^7 + O(w^8), 2 + 5 + 3*5^2 + 5^3 + 3*5^5 + 4*5^6 + O(5^8) + (1 + 2*5 + 4*5^2 + 5^4 + 4*5^5 + 2*5^6 + 2*5^7 + O(5^8))*w + (3 + 3*5^2 + 4*5^3 + 3*5^4 + 4*5^5 + 2*5^7 + O(5^8))*w^2 + (1 + 2*5^2 + 5^3 + 2*5^4 + 3*5^7 + O(5^8))*w^3 + (4 + 2*5 + 2*5^2 + 4*5^3 + 2*5^4 + 5^6 + 4*5^7 + O(5^8))*w^4 + (4 + 2*5^2 + 2*5^3 + 4*5^4 + 4*5^6 + 2*5^7 + O(5^8))*w^5 + (1 + 2*5 + 5^2 + 5^3 + 3*5^5 + 3*5^6 + O(5^8))*w^6 + (3 + 3*5 + 2*5^3 + 5^4 + O(5^8))*w^7 + O(w^8), 1 + 3*5^2 + 2*5^3 + 5^4 + 2*5^5 + 3*5^6 + 4*5^7 + O(5^8) + (2 + 2*5 + 4*5^2 + 4*5^3 + 2*5^4 + 2*5^5 + 5^7 + O(5^8))*w + (3 + 2*5 + 3*5^4 + 3*5^5 + 4*5^7 + O(5^8))*w^2 + (1 + 3*5 + 2*5^2 + 2*5^4 + 2*5^5 + 4*5^6 + 2*5^7 + O(5^8))*w^3 + (1 + 5 + 4*5^2 + 5^3 + 3*5^4 + 5^6 + O(5^8))*w^4 + (2 + 2*5 + 5^3 + 5^4 + 3*5^5 + 3*5^6 + 3*5^7 + O(5^8))*w^5 + (4 + 2*5 + 5^2 + 5^3 + 2*5^4 + 4*5^5 + 4*5^6 + 2*5^7 + O(5^8))*w^6 + (2 + 2*5 + 2*5^2 + 4*5^4 + 3*5^5 + 2*5^6 + 5^7 + O(5^8))*w^7 + O(w^8), 2*5 + 5^2 + 3*5^3 + 4*5^5 + 5^6 + 2*5^7 + O(5^8) + (2*5 + 4*5^2 + 2*5^3 + 5^4 + 2*5^5 + 3*5^6 + O(5^8))*w + (2 + 4*5 + 4*5^3 + 2*5^5 + 5^6 + 4*5^7 + O(5^8))*w^2 + (3 + 4*5 + 2*5^2 + 3*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + 2*5^7 + O(5^8))*w^3 + (4 + 4*5 + 3*5^2 + 5^4 + 2*5^5 + 5^6 + 2*5^7 + O(5^8))*w^4 + (1 + 3*5 + 5^2 + 5^3 + 4*5^4 + 5^5 + 3*5^6 + O(5^8))*w^5 + (3 + 3*5 + 5^2 + 4*5^3 + 4*5^4 + 5^5 + 2*5^6 + O(5^8))*w^6 + (1 + 5 + 2*5^3 + 2*5^4 + 4*5^5 + 5^6 + 3*5^7 + O(5^8))*w^7 + O(w^8), 2 + 2*5 + 5^2 + 4*5^3 + 2*5^5 + 3*5^6 + 5^7 + O(5^8) + (3 + 4*5 + 5^2 + 3*5^3 + 4*5^4 + 2*5^5 + 3*5^7 + O(5^8))*w + (3*5 + 2*5^2 + 2*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + O(5^8))*w^2 + (1 + 3*5 + 4*5^3 + 4*5^4 + 4*5^5 + 5^6 + 5^7 + O(5^8))*w^3 + (4 + 2*5 + 3*5^2 + 4*5^3 + 3*5^6 + 4*5^7 + O(5^8))*w^4 + (2*5 + 3*5^2 + 5^3 + 3*5^4 + 2*5^5 + 2*5^6 + 4*5^7 + O(5^8))*w^5 + (2 + 3*5 + 4*5^3 + 3*5^4 + 3*5^6 + 4*5^7 + O(5^8))*w^6 + (1 + 3*5 + 5^2 + 3*5^3 + 2*5^6 + 2*5^7 + O(5^8))*w^7 + O(w^8), 3 + 5^2 + 2*5^4 + 2*5^7 + O(5^8) + (2*5 + 5^2 + 5^4 + 5^6 + 3*5^7 + O(5^8))*w + (3 + 3*5^2 + 5^3 + 4*5^4 + 3*5^5 + 2*5^6 + 5^7 + O(5^8))*w^2 + (2 + 2*5 + 4*5^3 + 5^4 + 2*5^5 + 4*5^6 + 3*5^7 + O(5^8))*w^3 + (3 + 3*5 + 4*5^2 + 5^3 + 3*5^4 + 4*5^5 + 3*5^7 + O(5^8))*w^4 + (3 + 3*5 + 2*5^2 + 4*5^3 + 3*5^5+ 2*5^6 + O(5^8))*w^5 + (1 + 2*5 + 4*5^2 + 5^3 + 2*5^5 + 2*5^6 + 2*5^7 + O(5^8))*w^6 + (3 + 5 + 4*5^2 + 3*5^3 + 4*5^4 + 3*5^5+ 4*5^6 + O(5^8))*w^7 + O(w^8)))
            sage: mu = CoeffMod_OMS_Families_element(mu, DD, check=False, var_prec=8)
            sage: nu = mu.solve_diff_eqn()
            sage: S0 = DD.action().actor()
            sage: Delta_mat = S0([1,1,0,1])
            sage: nu * Delta_mat - nu - mu
            5^6 * () + O(w^8)
        """
        #Do something about ordp. Is this comment still relevant?
        p = self.parent().prime()
        abs_prec, var_prec = self.precision_absolute()
        if self.is_zero():
            V = self.parent().approx_module(0, var_prec)
            return CoeffMod_OMS_Families_element(V([]), self.parent(), ordp=(abs_prec - ZZ(abs_prec).exact_log(p) - 1), check=False, var_prec=var_prec)
        if not self._unscaled_moment(0).is_zero():
            raise ValueError("Family of distribution must have total measure 0 to be in image of difference operator; total measure is %s"%self.moment(0))
        M = ZZ(len(self._moments))
        if M == 2:
            if p == 2:
                raise ValueError("Not enough accuracy to return anything")
            else:
                #mu = self.parent()()
                #mu.ordp = self.ordp
                V = self.parent().approx_module(1, var_prec)
                return CoeffMod_OMS_Families_element(V([self._unscaled_moment(1), V.base_ring().base_ring()(0, 1)]), self.parent(), ordp=self.ordp, check=False, var_prec=var_prec)
        from sage.modular.pollack_stevens.coeffmod_OMS_space import OverconvergentDistributions
        DD = self.parent()
        ## Hack fix here -- the problem is that because we are using fixed weight distributions
        ## with the regular filtration, one needs more precision than the base ring has.
        M_in = _prec_for_solve_diff_eqn(M, p)
        R = Qp(DD.prime(), M_in)
        #R = self.base_ring().base_ring()
        D = OverconvergentDistributions(0, base=R, prec_cap=M_in, character=DD._character, adjuster=DD._adjuster, act_on_left=DD.action().is_left(), dettwist=DD._dettwist)
        V = D.approx_module(M_in)
        Elem = D.Element
        v = V([R.zero(), R.one()] + [R.zero()]*(M_in-2))
        mu = Elem(v, D, ordp=0, check=False)
        mus = self.moment(1) * mu.solve_diff_eqn().lift(DD)
        for j in range(2, M):
            mu._moments[j] = R.one()
            mu._moments[j-1] = R.zero()
            mus += self.moment(j) * mu.solve_diff_eqn().lift(DD)
        #prec = DD.length_reverse_lookup(M)
        #new_prec = prec - prec.exact_log(p) - 1 - mus.ordp
        #v = mus._moments[:DD.length_of_moments(new_prec)]
        #S = DD.base_ring()
        #v[len(v) - 1] = S(_add_big_ohs_list(v[len(v) - 1], [1, self._var_prec]))
        #mus = CoeffMod_OMS_Families_element(v, self.parent(), ordp=mus.ordp, check=False, var_prec=self._var_prec)
        #Should we remove precision like at end of non-family code, or is this taken care of?
        verbose("Abs_prec of self: %s"%(abs_prec))
        verbose("mus_abs_prec before reduction: %s"%(mus.precision_absolute()))
        from sage.structure.sage_object import dumps
        verbose("self dump: %s"%(repr(dumps(self))))
        mus = mus.reduce_precision_absolute([abs_prec - abs_prec.exact_log(p) - 1, var_prec])
        verbose("mus_abs_prec after reduction: %s"%(mus.precision_absolute()))
        return mus.normalize()  #Is it necessary to normalize?

def _padic_val_of_pow_series(f, p=None, var_prec=None):
    r"""
        Given a power series ``f`` return its ``p``-adic valuation, i.e. the
        minimum ``p``-adic valuation of its coefficients
    """
    if f == 0:
        return Infinity
    if var_prec is None:
        return min([coeff.valuation() if not coeff.is_zero() else coeff.precision_absolute() for coeff in f])
    flist = f.padded_list()[:var_prec]
    return min([coeff.valuation() if not coeff.is_zero() else coeff.precision_absolute() for coeff in flist])

def _padic_val_unit_of_pow_series(f, p=None):
    r"""
        Given a power series ``f``, return a tuple `(v, u)`, where `v` is the
        ``p``-adic valuation of ``f`` and `u` is the ``p``-adic unit part of
        ``f``.
    """
    if f == 0:
        return (Infinity, self.parent()(0,0))
    if p is None:
        p = f.parent().base_ring().prime()
    v = _padic_val_of_pow_series(f, p)
    u = f.parent()([(coeff / (p ** v)) for coeff in f])
    return (v, u)

def _padic_abs_prec_of_pow_series(f, var_prec=None):
    r"""
    """
    if var_prec is None:
        var_prec = f._var_prec
    flist = f.padded_list()
    if len(flist) == 0:
        return Infinity
    return min([flist[i].precision_absolute() for i in range(min(len(flist), var_prec))])

def _is_zero_padic_power_series(f, prec):
    r"""
    Determines whether f is zero with `w`-adic precision prec[1], and `p`-adic
    precision prec[0].
    
    EXAMPLES::
    
        sage: from sage.modular.pollack_stevens.coeffmod_OMS_families_element import _is_zero_padic_power_series
        sage: R = PowerSeriesRing(Zp(5), 'w')
        sage: f = R([0,1,3,-2], 6)
        sage: _is_zero_padic_power_series(f, [5,4])
        False
        sage: _is_zero_padic_power_series(f, [5,7])
        False
        sage: _is_zero_padic_power_series(f, [25,7])
        False
        sage: _is_zero_padic_power_series(f, [25,1])
        True
        sage: f2 = R([O(5^2),1,3,-2],6)
        sage: _is_zero_padic_power_series(f2, [5,1])
        False
        sage: _is_zero_padic_power_series(f2, [2,1])
        True
        sage: f3 = R(0, 3)
        sage: _is_zero_padic_power_series(f3, [100, 3])
        True
        sage: _is_zero_padic_power_series(f3, [100, 4])
        False

    """
    if f.prec() < prec[1]:
        return False
    flist = _add_big_ohs_list(f, prec)
    for c in flist:
        try:
            if not c.is_zero(prec[0]):
                return False
        except PrecisionError:
            return False
    return True

def _add_big_ohs_list(f, prec_cap):
    #There may be a problem if you pass the list [3, O(p^2)] to the power series ring. It will truncate all big-ohs after last "non-zero" entry
    r"""
    Returns a (padded) list of length (at most) ``prec_cap``[1] of the coefficients
    up to `p`-adic precision ``prec_cap``[0]. The input is checked. 
    
    INPUT:
    
        - ``f`` -- a power series over a `p`-adic ring
        - ``prec_cap`` -- a pair [``p_prec``, ``var_prec``], where ``p_prec`` is
        the desired `p`-adic precision of the coefficients of f and ``var_prec``
        is the desired variable-adic precision
    
    OUTPUT:
    
        - a list of coefficients of ``f`` of precision up to (at most)
        `O(p^{\text{p_prec}})`. The length of this list is the minimum of the
        precision of ``f`` and of ``var_prec``. The padded part should be exact zeroes.
    
    EXAMPLES::
    
        sage: from sage.modular.pollack_stevens.coeffmod_OMS_families_element import _add_big_ohs_list
        sage: R = PowerSeriesRing(Zp(3), 'w')
        sage: f = R([0,-1,2,0,-3])
        sage: _add_big_ohs_list(f, [10,7])
        [O(3^10), 2 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + O(3^10), 2 + O(3^10), O(3^10), 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + O(3^10)]
        sage: _add_big_ohs_list(f, [2,4])
        [O(3^2), 2 + 2*3 + O(3^2), 2 + O(3^2), O(3^2)]
    """
    p_prec, var_prec = prec_cap
    flist = f.padded_list()
    return [flist[i].add_bigoh(p_prec) for i in range(min(len(flist), var_prec))]

def _shift_coeffs(f, shift, right=True):
    r"""
        Given a power series ``f``, apply '>> shift' to each of its coefficients,
        i.e. divide each by shift powers of the uniformizer.
    """
    if shift == 0:
        return f
    # RP: Function w/o this fix seems to be failing when f = 0 and returning O(w^0) 
    if f == 0:
        return f
    if not right:
        shift = -shift
    flist = [a >> shift for a in f.padded_list()]
    # Hack to circumvent a bug in sage's handling of power series/polynomials
    # over p-adic rings
    R = f.parent()
    if R.base_ring().is_field():
        return R(flist)
    RP = R._poly_ring()
    from sage.rings.infinity import infinity
    absprec = min([infinity] + [a.precision_absolute() for a in flist])
    from sage.rings.polynomial.polynomial_element import Polynomial_generic_dense
    fpoly = Polynomial_generic_dense(RP, x=flist, check=False, absprec=absprec)
    return R.element_class(R, f=fpoly, prec=len(flist), check=False)

def _sanitize_alpha(alpha):
    RK = alpha.parent()
    alpha_list = alpha.list()
    len_alpha = len(alpha_list)
    #Remove O(p^negative)
    i = 0
    while i < len_alpha:
        try:
            alpha_list[i].is_zero(1)
        except PrecisionError:
            break
        i += 1
    if i == len_alpha:
        return alpha
    return RK(alpha_list[:i], i)

def _truncate_low_prec(alpha):
    alist = alpha.padded_list()
    for i in range(len(alist)):
        if alist[i].precision_absolute() <= 0:
            return alpha.add_bigoh(i)
    return alpha
#Hack to overcome lack of functionality in dividing power series
#def _custom_ps_div(num, denom):
#    #print "\n", num
#    #print denom
#    num_val = num.valuation()
#    denom_val = denom.valuation()
#    if denom_val == 0:
#        prec = min(num.prec(), denom.prec())
#        #print "\Val (0,0)"
#        #print num[0]
#        #print num
#        R = num.parent()
#        num_list = num.padded_list(prec)
#        denom_list = denom.padded_list(prec)
#        outlist = [num[0]/denom[0]]
#        for i in range(1, prec):
#            outlist.append(((num_list[i] - sum([denom_list[j] * outlist[i-j] for j in range(1, i + 1)])) / denom_list[0])) #.add_bigoh(prec - i)) #Why did I write this last part?
#        ret = R(outlist, prec)
#        #Hack around bug in power series
#        #d = ret.degree()
#        #missing = prec - d - 1
#        #if missing == 0:
#        #    return ret
#        #w = R.gen()
#        #B = R.base_ring()
#        #for i in range(missing):
#        #    shift = d + 1 + i
#        #    ret += B(0, prec - shift) * w ** (shift)
#        return [ret, prec]
#    if denom_val > num_val:
#        raise ValueError("Denominator must have smaller valuation than numerator.")
#    if denom_val == num_val:
#        alpha, new_prec = _custom_ps_div(num.valuation_zero_part(), denom.valuation_zero_part())
#        return [alpha, new_prec - denom_val]
#    val_diff = num_val - denom_val
#    alpha, new_prec = _custom_ps_div(num >> val_diff, denom >> val_diff)
#    return [alpha << val_diff, new_prec - val_diff]

#def create__CoeffMod_OMS_Families_element(moments, parent, ordp, var_prec):
#    """
#    Used for unpickling.
#    """
#    return CoeffMod_OMS_Families_element(moments, parent, ordp=ordp, check=False, var_prec=var_prec)

def test_correctness_and_precision_of_solve_diff_eqn_fam(number=20, verbosity=1):
    """
    ``number`` is how many different random distributions to check. 
    
    Currently, avoids the prime 2.
    """
    from sage.misc.prandom import randint
    from sage.rings.arith import random_prime
    from sage.rings.padics.factory import ZpCA
    from sage.modular.pollack_stevens.coeffmod_OMS_families_space import FamiliesOfOverconvergentDistributions
    from sage.structure.sage_object import dumps
    errors = []
    munus = []
    for i in range(number):
        Mspace = randint(1, 20)    #Moments of space
        M = randint(max(0, Mspace - 5), Mspace)
        var_prec = randint(1, 8)
        p = random_prime(13, lbound=3)
        k = randint(0, 6)
        Rprec = Mspace + randint(0, 5)
        R = ZpCA(p, Rprec)
        D = FamiliesOfOverconvergentDistributions(k, base_coeffs=R, prec_cap=[Mspace, var_prec])
        S0 = D.action().actor()
        Delta_mat = S0([1,1,0,1])
        mu = D.random_element(M)
        mu_save = dumps(mu)#[deepcopy(mu.ordp), deepcopy(mu._moments)]
        if verbosity > 0:
            print "\nTest #{0} data (Mspace, M, p, k, Rprec, var_prec) =".format(i+1), (Mspace, M, p, k, Rprec, var_prec)
            print "mu =", mu
        
        nu = mu * Delta_mat - mu
        nu_save = [deepcopy(nu.ordp), deepcopy(nu._moments)]
        mu2 = nu.solve_diff_eqn()
        nu_abs_prec = nu.precision_absolute()[0]
        expected = nu_abs_prec - nu_abs_prec.exact_log(p) - 1
        if M != 1:
            try:
                agree = (mu - mu2).is_zero(expected)
            except PrecisionError:
                print (Mspace, M, p, k, Rprec), mu_save._repr_(), nu_save
                assert False
        else:
            agree = mu2.is_zero(expected)
        if verbosity > 1:
            print "    Just so you know:"
            #print "     mured =", mu.reduce_precision_absolute(expected)
            print "       mu2 =", mu2
            print "        nu = ", nu
        if not agree:
            errors.append((i+1, 1))
            munus.append((mu_save, nu_save, mu2, (Mspace, M, p, k, Rprec, var_prec)))
        if verbosity > 0:
            print "    Test finding mu from mu|Delta accurate: %s"%(agree)
            print "        nu_abs_prec  soln_abs_prec_expected  actual  agree"
        mu2_abs_prec = mu2.precision_absolute()[0]
        agree = (expected == mu2_abs_prec)
        if verbosity > 0:
            print "        %s             %s                       %s      %s"%(nu_abs_prec, expected, mu2_abs_prec, agree)
        if not agree:
            errors.append((i+1, 2))
            munus.append((mu_save, nu_save, mu2, (Mspace, M, p, k, Rprec, var_prec)))
        
        if mu.precision_relative()[0] > 0:
            mu._moments[0] = R(0, mu.precision_relative()[0])
        mu_save = [deepcopy(mu.ordp), deepcopy(mu._moments)]
        if verbosity > 0:
            print "    mu modified =", mu
        nu = mu.solve_diff_eqn()
        mu_abs_prec = mu.precision_absolute()[0]
        expected = mu_abs_prec - mu_abs_prec.exact_log(p) - 1
        nud = nu * Delta_mat - nu
        nu_save = [deepcopy(nu.ordp), deepcopy(nu._moments)]
        agree = (nud - mu).is_zero(expected)
        if verbosity > 1:
            print "    Just so you know:"
            print "        mu =", mu
            #print "     mured =", mu.reduce_precision_absolute(expected)
            print "       nud =", nud
        if not agree:
            errors.append((i+1, 3))
            munus.append((mu_save, nu_save, (Mspace, M, p, k, Rprec, var_prec)))
        if verbosity > 0:
            print "    Test finding nu with nu|Delta == mu: %s"%(agree)
            print "        mu_abs_prec  soln_abs_prec_expected  actual  agree"
        nu_abs_prec = nu.precision_absolute()[0]
        agree = (expected == nu_abs_prec)
        if verbosity > 0:
            print "        %s             %s                       %s      %s"%(mu_abs_prec, expected, nu_abs_prec, agree)
        if not agree:
            errors.append((i+1, 4))
            munus.append((mu_save, nu_save, (Mspace, M, p, k, Rprec, var_prec)))
    if len(errors) == 0:
        if verbosity > 0:
            print "\nTest passed with no errors."
        return
    if verbosity > 0:
        print "\nTest failed with errors: %s\n"%(errors)
    return errors, munus
