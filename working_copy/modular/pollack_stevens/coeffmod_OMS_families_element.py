# This should be cythoned once it's done.

#from sage.rings.power_series_ring import PowerSeriesRing

from copy import copy   #Not necessary in cython version
from sage.misc.misc import verbose
from sage.rings.infinity import Infinity
from sage.modular.pollack_stevens.coeffmod_element import CoefficientModuleElement_generic
from sage.modular.pollack_stevens.coeffmod_OMS_element import CoeffMod_OMS_element
from sage.rings.integer_ring import ZZ
from sage.functions.other import ceil

class CoeffMod_OMS_Families_element(CoefficientModuleElement_generic):
    r"""
        Fill this in.
        
        EXAMPLES::
        
            sage: D8 = FamiliesOfOverconvergentDistributions(0, prec_cap = [8 ,4], base_coeffs=ZpCA(3, 8))
            sage: mu8 = D8([1,2,3,4,5,6,7,8,9,10]); mu8
            (1 + O(3^4), 2 + O(3^4), 3 + O(3^3), 1 + 3 + O(3^3), 2 + 3 + O(3^2), 2*3 + O(3^2), 1 + O(3), 2 + O(3)) + O(w^4)
            sage: D4 = FamiliesOfOverconvergentDistributions(0, prec_cap = [8 ,4], base_coeffs=ZpCA(3, 4))
            sage: mu4 = D4([1,2,3,4,5,6,7,8,9,10]); mu4
            (1 + O(3^2), 2 + O(3^2), 0, 1 + O(3)) + O(w^4)
            sage: D4(mu8)
            (1 + O(3^2), 2 + O(3^2), 0, 1 + O(3)) + O(w^4)
            sage: mu4 == D4(mu8)
            True
            sage: D42 = FamiliesOfOverconvergentDistributions(0, prec_cap = [8 ,2], base_coeffs=ZpCA(3, 4))
            sage: mu42 = D42([1,2,3,4,5,6,7,8,9,10])
            sage: mu42
            (1 + O(3^2), 2 + O(3^2), 0, 1 + O(3)) + O(w^2)
            sage: D42(mu8)
            (1 + O(3^2), 2 + O(3^2), 0, 1 + O(3)) + O(w^2)
            sage: mu42 == D42(mu8)
            True
            sage: D42(15)
            3 * (2 + O(3)) + O(w^2)
            sage: D = FamiliesOfOverconvergentDistributions(2, prec_cap = [8 ,4], base_coeffs=ZpCA(11, 4))
            sage: R = D.base_ring(); K = R.base_extend(R.base_ring().fraction_field())
            sage: v = [K([1,2,11]) / 11, K([1]), K([11,1,1])]
            sage: D(v)
            11^-1 * (1 + O(11^3) + (2 + O(11^3))*w + (11 + O(11^3))*w^2, 11 + O(11^2), 0) + O(w^4)
            sage: D(0)
            11^4 * () + O(w^4)
    """
    def __init__(self, moments, parent, ordp=0, check=True, var_prec=None):
        CoefficientModuleElement_generic.__init__(self, parent)
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
                moments = moments._moments[:par_p_prec]
                moments = moments.change_ring(parent.base_ring())
            if isinstance(moments, CoeffMod_OMS_element):
                # Coerce in using constant family
                ordp = moments.ordp
                p_prec = min(moments.precision_relative(), par_p_prec)
                moments = self.parent().approx_module(p_prec, var_prec)(moments._moments)
            elif hasattr(moments, '__len__'):
                #Need to modify if uniformiser is not p
                #Deal with var_prec
                R = parent.base_ring()
                K = R.base_extend(R.base_ring().fraction_field())
                p_prec = min(len(moments), par_p_prec)
                #figure out var_prec from moments
                V = parent.approx_module(p_prec, var_prec)
                VK = V.base_extend(K)
                moments = VK(moments[:p_prec])
                if len(moments) == 0 or moments == 0:   #should do something with how "zero" moments is
                    V = parent.approx_module(0, var_prec) #var_prec?
                    moments = V([])
                    ordp = par_p_prec
                else:
                    ordp = min([_padic_val_of_pow_series(a) for a in moments])
                    moments = V([_right_shift_coeffs(a, ordp) for a in moments])
            elif moments == 0:  #should do something with how "zero" moments is
                V = parent.approx_module(0, var_prec)
                moments = V([])
                ordp = par_p_prec
            else:
                R = parent.base_ring()
                K = R.base_extend(R.base_ring().fraction_field())
                V = parent.approx_module(1, var_prec)
                moments = K(moments)
                ordp = _padic_val_of_pow_series(moments)
                if ordp != 0:
                    moments = V([_right_shift_coeffs(moments, ordp)])
                else:
                    moments = V([moments])
        self._moments = moments
        #if var_prec is None:
        #    self._var_prec = parent._prec_cap[1]
        #else:
        #    self._var_prec = var_prec
        self.ordp = ordp
        self._var_prec = var_prec
    
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
        #print "In _addsub"
        #print saprec, raprec
        aprec = min(saprec[0], raprec[0])
        #print "aprec", aprec
        ans.ordp = min(self.ordp, right.ordp)
        #print self.ordp, right.ordp
        #print "ans.ordp", ans.ordp
        rprec = aprec - ans.ordp
        #print rprec
        #print self.parent().precision_cap()
        var_prec = min(saprec[1], raprec[1])
        V = ans.parent().approx_module(rprec, var_prec)
        R = V.base_ring()
        smoments = self._moments
        rmoments = right._moments
        if smoments.parent() is not V:
            smoments = V(smoments.list(copy=False)[:rprec] + ([R(0)] * (rprec - len(smoments)) if rprec > len(smoments) else []))
        if rmoments.parent() is not V:
            rmoments = V(rmoments.list(copy=False)[:rprec] + ([R(0)] * (rprec - len(rmoments)) if rprec > len(rmoments) else []))
        # We multiply by the relative power of p
        if self.ordp > right.ordp:
            smoments *= self.parent().prime()**(self.ordp - right.ordp)
        elif self.ordp < right.ordp:
            rmoments *= self.parent().prime()**(right.ordp - self.ordp)
        if negate:
            rmoments = -rmoments
        ans._moments = smoments + rmoments
        return ans
    
    def _lmul_(self, right):
        #RH: adapted from coeffmod_OMS_element.py
        """
        Scalar multiplication self*right.
        """
        ans = self.parent()(0)
        if right.is_zero():
            return ans
        p = self.parent().prime()
        if right.is_zero():
            ans._moments = self.parent().approx_module(0, 0)([])
            ans.ordp = self.parent().precision_cap()[0]
        else:
            v, u = _padic_val_unit_of_pow_series(right, p)
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
        #RH: adapted from coeffmod_OMS_element.py
        left = copy(left)
        right = copy(right)
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
        rprec = min(lrprec, rrprec)
        #if lrprec > rprec:
        #    left.reduce_precision(rprec)
        #elif rrprec > rprec:
        #    right.reduce_precision(rprec)
        #left.normalize()
        #right.normalize()
        p = left.parent().prime()
        if left.ordp > right.ordp:
            shift = p ** (left.ordp - right.ordp)
            for i in range(rprec):
                L = (shift * left._unscaled_moment(i)).padded_list()
                R = (right._unscaled_moment(i)).padded_list()
                c = cmp(L, R)
                if c: return c
        elif left.ordp < right.ordp:
            shift = p ** (right.ordp - left.ordp)
            for i in range(rprec):
                L = (left._unscaled_moment(i)).padded_list()
                R = (shift * right._unscaled_moment(i)).padded_list()
                c = cmp(L, R)
                if c: return c
        else:
            for i in range(rprec):
                L = (left._unscaled_moment(i)).padded_list()
                R = (right._unscaled_moment(i)).padded_list()
                c = cmp(L, R)
                if c: return c
        return 0

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
    
    def _coerce_map_from_(self, other):
        if isinstance(other, CoeffMod_OMS_Families_element) \
            and other._k  == self._k \
            and self._character == other._character \
            and self.base_ring().has_coerce_map_from(other.base_ring()):
            return True
        elif isinstance(other, CoeffMod_OMS_element):
            #kdiff = other._k - self._k
            #self._character = other._character
            #and self.base_ring().has_coerce_map_from(other.base_ring()):
            return False
        else:
            return False
    
    def is_zero(self, prec=None):
        s_prec_rel = self.precision_relative()
        if prec is None:
            prec = s_prec_rel
        elif not hasattr(moments, '__len__'):
            prec = [ZZ(prec), None]
        elif len(prec) > 2:
            raise TypeError("prec must have length at most 2.")
        if prec[0] > s_prec_rel[0] or prec[1] > s_prec_rel[1]:
            return False
        num_moments = min(prec[0], s_prec_rel[0]) if prec[0] is not None else s_prec_rel[0]
        var_prec = min(prec[1], s_prec_rel[1]) if prec[1] is not None else s_prec_rel[1]
        p_precs = self.parent().filtration_precisions(num_moments)
        for i in range(num_moments):
            selflist = _add_big_ohs_list(self._moments[i], [p_precs[i], var_prec])
            for c in selflist:
                if c != 0:
                    return False
        return True
    
    def find_scalar(self, other, M = None, check=True):
        self.normalize()
        other.normalize()
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
        verbose("a = %s"%(a))
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
        relprec = n - i - v
        var_prec = min(s_var_prec, other_var_prec)  #should this depend on w-adic valuation?
        Rbase = R.base_ring()
        RK = R.change_ring(Rbase.fraction_field())
        if i < other_pr:
            #Hack
            verbose("val, other._unscaled_moment(%s) = %s, %s"%(i, other.ordp, other._moments[i]))
            #alpha = _sanitize_alpha(RK(other._unscaled_moment(i)) / RK(a))
            alpha = RK(other._unscaled_moment(i)) / RK(a)
            #alpha, var_prec_loss = _custom_ps_div(other._unscaled_moment(i), a)
            #alpha = R(_add_big_ohs_list(alpha, [ceil((n - i) * self._cp), min(s_var_prec, other_var_prec)]))
            #alpha = R(_add_big_ohs_list(other._unscaled_moment(i) / a, [ceil((n - i) * self._cp), min(s_var_prec, other_var_prec)]))
        else:
            #Fix var_prec??
            alpha = RK([Rbase(0, p_precs[i])], var_prec)
        verbose("alpha = %s"%(alpha))
        
        while i < other_pr-1:
            i += 1
            verbose("comparing p moment %s"%i)
            a = self._unscaled_moment(i)
            if check:
#                   verbose("self.moment=%s, other.moment=%s"%(a, other._unscaled_moment(i)))
                #if other._unscaled_moment(i) != _add_big_ohs_list(alpha * a, [ceil((n - i) * self._cp), var_prec]):
                verbose("val, self._unscaled_moment(%s) = %s, %s"%(i, self.ordp, other._moments[i]))
                verbose("val, other._unscaled_moment(%s) = %s, %s"%(i, other.ordp, other._moments[i]))
                if other._unscaled_moment(i) != alpha * a:
                    raise ValueError("not a scalar multiple")
            v = _padic_val_of_pow_series(a, p)
            if n - i - v > relprec:
                verbose("Resetting alpha: relprec=%s, n-i=%s, v=%s"%(relprec, n-i, v))
                relprec = n - i - v
                #Should we alter var_prec, too?
                #Hack
                alpha = RK(other._unscaled_moment(i)) / RK(a)
                #alpha = R(_add_big_ohs_list(other._unscaled_moment(i) / a, [ceil((n - i) * self._cp), min(s_var_prec, other_var_prec)]))
                verbose("alpha=%s"%(alpha))
        if relprec < M[0]:
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
        #RH: copied from coeffmod_OMS_element.py
        return [ZZ(len(self._moments)), self._var_prec]
    
    def precision_absolute(self):
        #RH: copied from coeffmod_OMS_element.py
        return [ZZ(len(self._moments) + self.ordp), self._var_prec]
    
    def valuation(self):
        #RH: adapted from coeffmod_OMS_element.py
        p = self.parent().prime()
        n = self.precision_relative()[0]
        if n == 0:
            return self.ordp
        return self.ordp + min([n] + [_padic_val_of_pow_series(self._unscaled_moment(a), p) for a in range(n) if not self._unscaled_moment(a).is_zero()])
    
    def normalize(self):
        #RH: adapted from coeffmod_OMS_element.py
        #Not tested
        V = self._moments.parent()
        R = V.base_ring()
        n, v_prec = self.precision_relative()
        p = self.parent().prime()
        self_val = self.valuation()
        shift = self_val - self.ordp
        self.ordp = self_val
        #Factor out powers of uniformizer and check precision
        p_precs = self.parent().filtration_precisions(n)
        adjust_moms = 0
        verbose("n: %s; shift: %s; _mom: %s\np_precs: %s"%(n, shift, self._moments, p_precs), level=2)
        if shift > 0:
            for i in range(n):
                self._moments[i] = _right_shift_coeffs(self._moments[i], shift)
                adjust_moms = max(adjust_moms, p_precs[i] - _padic_abs_prec_of_pow_series(self._moments[i], v_prec))
        elif shift == 0:
            for i in range(n):
                adjust_moms = max(adjust_moms, p_precs[i] - _padic_abs_prec_of_pow_series(self._moments[i], v_prec))
        else:
            raise NotImplementedError("Currently only deals with the case where the base coefficients are a ring of integers.")
        #Cut down moments because of precision loss
        verbose("adjust_mom: %s\n_moms: %s"%(adjust_moms, self._moments), level=2)
        if adjust_moms >=n:
            V = self.parent().approx_module(0)
            self._moments = V([])
            self._ordp = adjust_moms
        elif adjust_moms > 0:
            n -= adjust_moms    #Is this going to give the correct precision?
            p_precs = self.parent().filtration_precisions(n)
            verbose("new p_precs: %s"%(p_precs), level=2)
            R = self.parent().base_ring()
            for i in range(n):
                self._moments[i] = R(_add_big_ohs_list(self._moments[i], [p_precs[i], v_prec]), v_prec)
        else:
            R = self.parent().base_ring()
            for i in range(n):
                self._moments[i] = R(_add_big_ohs_list(self._moments[i], [p_precs[i], v_prec]), v_prec)
        return self

    def moment(self, n):
        #RH: "copied" from dist.pyx
        r"""
        Returns the `n`-th moment.
        """
        if self.ordp == 0:
            return self._unscaled_moment(n)
        else:
            return self.parent().prime()**(self.ordp) * self._unscaled_moment(n)
    
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
            raise ValuError("Precisions specified must be less than current precisions.")
        moments = self._moments[:p_prec]
        ordp = self.ordp
        return CoeffMod_OMS_Families_element(moments, self.parent(), ordp=ordp, check=False, var_prec=var_prec)
    
    def solve_diff_eqn(self):
        #Do something about ordp
        if self.is_zero():
            return self.parent()(0)
        if self._moments[0] != 0:
            raise ValueError("Family of distribution must have total measure 0 to be in image of difference operator.")
        M = len(self._moments)
        if M == 1:
            return self.parent().zero()
        from sage.modular.pollack_stevens.coeffmod_OMS_space import OverconvergentDistributions
        R = self.base_ring().base_ring()
        DD = self.parent()
        D = OverconvergentDistributions(0, base=R, prec_cap=M, character=DD._character, adjuster=DD._adjuster, act_on_left=DD.action().is_left(), dettwist=DD._dettwist)
        V = D.approx_module(M)
        Elem = D.Element
        v = V([R.zero(), R.one()] + [R.zero()]*(M-2))
        mu = Elem(v, D, ordp=0, check=False)
        #Normalize mu?
        mus = self.moment(1) * mu.solve_diff_eqn().lift(DD)
        for j in range(2, M):
            mu._moments[j] = R.one()
            mu._moments[j-1] = R.zero()
            mus += self.moment(j) * mu.solve_diff_eqn().lift(DD)
        return mus

def _padic_val_of_pow_series(f, p=None):
    r"""
        Given a power series ``f`` return its ``p``-adic valuation, i.e. the
        minimum ``p``-adic valuation of its coefficients
    """
    if f == 0:
        return Infinity
    return min([coeff.valuation() if not coeff.is_zero() else coeff.precision_absolute() for coeff in f])

def _padic_val_unit_of_pow_series(f, p=None):
    r"""
        Given a power series ``f``, return a tuple `(v, u)`, where `v` is the
        ``p``-adic valuation of ``f`` and `u` is the ``p``-adic unit part of
        ``f``.
    """
    if f == 0:
        return (Infinity, self.parent()(0,0))
    if p is None:
        p = self.parent().base_ring().prime()
    v = _padic_val_of_pow_series(f, p)
    u = f.parent()([(coeff / (p ** v)) for coeff in f])
    return (v, u)

def _padic_abs_prec_of_pow_series(f, var_prec=None):
    r"""
    """
    if var_prec is None:
        var_prec = f.parent().default_prec()
    flist = f.padded_list()
    if len(flist) == 0:
        return Infinity
    return min([flist[i].precision_absolute() for i in range(min(len(flist), var_prec))])

def _add_big_ohs_list(f, prec_cap):
    #There may be a problem if you pass the list [3, O(p^2)] to the power series ring. It will truncate all big-ohs after last "non-zero" entry
    r"""
    Returns a (padded) list of length (at most) ``prec_cap``[1] of the coefficients
    up to `p`-adic precision ``prec_cap``[0]. The input is checked. 
    
    INPUT:
    
        - ``f`` -- a power series over a `p`-adic ring
        - ``prec_cap`` -- a pair [``p_prec``, ``var_prec``, where ``p_prec`` is
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

def _right_shift_coeffs(f, shift):
    r"""
        Given a power series ``f``, apply '>> shift' to each of its coefficients,
        i.e. divide each by shift powers of the uniformizer.
    """
    if shift == 0:
        return f
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