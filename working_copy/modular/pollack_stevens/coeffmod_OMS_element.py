# This should be cythoned once it's done.

from copy import copy, deepcopy   #Not necessary in cython version
from sage.misc.misc import verbose
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.padics.precision_error import PrecisionError
from sage.rings.arith import binomial, bernoulli
from sage.modular.pollack_stevens.coeffmod_element import CoefficientModuleElement_generic

#maxordp = 2 ** 14   #change this to what it is in dist.pyx

class CoeffMod_OMS_element(CoefficientModuleElement_generic):
    r"""
        Fill this in.
        
        EXAMPLES::
        
            sage: D8 = OverconvergentDistributions(0, base=ZpCA(11, 8))
            sage: mu8 = D8([1,2,3,4,5,6,7,8,9,10]); mu8
            (1 + O(11^8), 2 + O(11^7), 3 + O(11^6), 4 + O(11^5), 5 + O(11^4), 6 + O(11^3), 7 + O(11^2), 8 + O(11))
            sage: D4 = OverconvergentDistributions(0, base=ZpCA(11, 4))
            sage: mu4 = D4([1,2,3,4,5,6,7,8,9,10]); mu4
            (1 + O(11^4), 2 + O(11^3), 3 + O(11^2), 4 + O(11))
            sage: D4(mu8)
            (1 + O(11^4), 2 + O(11^3), 3 + O(11^2), 4 + O(11))
            sage: mu4 == D4(mu8)
            True
            sage: D = OverconvergentDistributions(0,base=ZpCA(3,5))
            sage: x = D([3,3,3,3]); x._moments
            (1 + O(3^5), 1 + O(3^5), 1 + O(3^5), 1 + O(3^5))
            sage: x
            3 * (1 + O(3^4), 1 + O(3^3), 1 + O(3^2), 1 + O(3))
            sage: z = 3 * x; z._moments
            (1 + O(3^4), 1 + O(3^3), 1 + O(3^2), 1 + O(3))
            sage: z
            3^2 * (1 + O(3^4), 1 + O(3^3), 1 + O(3^2), 1 + O(3))
            sage: y = D([9,9,9,9]); y._moments
            (1 + O(3^5), 1 + O(3^5), 1 + O(3^5), 1 + O(3^5))
            sage: y
            3^2 * (1 + O(3^4), 1 + O(3^3), 1 + O(3^2), 1 + O(3))
            sage: y == z
            True
            sage: z == y
            True
            sage: a = D([81,81,81,81]); a._moments
            (1 + O(3^5), 1 + O(3^5), 1 + O(3^5), 1 + O(3^5))
            sage: a
            3^4 * (1 + O(3^4), 1 + O(3^3), 1 + O(3^2), 1 + O(3))
            sage: b = 27 * x
            sage: b._moments    #because in ZpCA with prec=5 27 = 3^3 * (1 + O(3^2))
            (1 + O(3^2), 1 + O(3^2), 1 + O(3^2), 1 + O(3))
            sage: b
            3^4 * (1 + O(3^2), 1 + O(3))
            sage: a == b
            True
            sage: b == a
            True
            sage: D(15)
            3 * (2 + O(3))
            sage: D = OverconvergentDistributions(0, base=ZpCA(5, 8))
            sage: D([25+O(5^2), 5+O(5)])
            5^2 * ()
            sage: D([25, 5])
            5 * (5 + O(5^2), 1 + O(5))
    
    TEST::
    
        sage: TestSuite(mu8).run()
        sage: from sage.modular.pollack_stevens.coeffmod_OMS_element import test_correctness_and_precision_of_solve_diff_eqn
        sage: test_correctness_and_precision_of_solve_diff_eqn(10, 0)
    """
    #RH: copied from dist.pyx (fixed dealing with 0)
    def __init__(self, moments, parent, ordp=0, check=True):
        CoefficientModuleElement_generic.__init__(self, parent)
        #TODO: finish this
        if check:
            if isinstance(moments, CoeffMod_OMS_element):
                ordp = moments.ordp
                moments = moments._moments[:parent.precision_cap()]
                moments = moments.change_ring(parent.base_ring())
            elif hasattr(moments, '__len__'):
                R = parent.base_ring()
                K = R.fraction_field()
                M = min(len(moments), parent.precision_cap())
                V = parent.approx_module(M)
                VK = V.base_extend(K)
                moments = VK(moments[:M])
                if len(moments) == 0 or moments == 0:
                    V = parent.approx_module(0)
                    if len(moments) > 0:
                        #Determine ordp
                        m = len(moments)
                        diff = 0
                        for i in range(m):
                            diff = max(diff, m - i - moments[i].precision_absolute())
                        ordp = m - diff
                    else:
                        ordp = parent.precision_cap()
                    moments = V([])
                else:
                    ordp = min([a.valuation() for a in moments])
                    for i in range(M):
                        moments[i] = moments[i] >> ordp
                    moments = V(moments)
            elif moments == 0:  #should do something with how "zero" moments is
                V = parent.approx_module(0)
                moments = V([])
                ordp = parent.precision_cap()
            else:
                K = parent.base_ring().fraction_field()
                V = parent.approx_module(1)
                moments = K(moments)
                ordp = moments.valuation()
                moments = V([moments >> ordp])
        
        self._moments = moments
        #if var_prec is None:
        #    self._var_prec = parent._prec_cap[1]
        #else:
        #    self._var_prec = var_prec
        self.ordp = ZZ(ordp)
    
    #def __reduce__(self):
#        """
#        TESTS::
#        
#            sage: D = OverconvergentDistributions(4, 3, 10)
#            sage: mu = D.random_element()
#            sage: loads(dumps(mu)) == mu
#            True
#        """
#        from sage.modular.pollack_stevens.coeffmod_OMS_element import create__CoeffMod_OMS_element
#        return (create__CoeffMod_OMS_element, (self._moments, self.parent(), self.ordp))
    
    #def _relprec(self):
    #    return len(self._moments)
    def _repr_(self):
        #RH: "copied" from dist.pyx
        self.normalize()
        valstr = ""
        if self.ordp == 1:
            valstr = "%s * "%(self.parent().prime())
        elif self.ordp != 0:
            valstr = "%s^%s * "%(self.parent().prime(), self.ordp)
        return valstr + repr(self._moments)
    
    def _add_(self, right):
        #RH: "copied" from dist.pyx
        # Should check by hand that the following examples do actually return
        # the correct answer
        r"""
            EXAMPLES::
            
                sage: D = OverconvergentDistributions(2,base=ZpCA(7,5))
                sage: mu = D((4 + 5*7 + 4*7^2 + 2*7^3 + O(7^5), 3 + 4*7 + 2*7^3 + O(7^4), 3 + 7 + 5*7^2 + O(7^3), 4 + 5*7 + O(7^2), 1 + O(7)))
                sage: nu = D((5 + 2*7 + 6*7^2 + 3*7^3 + 3*7^4 + O(7^5), 6*7 + 5*7^2 + 2*7^3 + O(7^4), 4 + 4*7 + O(7^3), 3 + 5*7 + O(7^2), 4 + O(7)))
                sage: mu + nu
                (2 + 7 + 4*7^2 + 6*7^3 + 3*7^4 + O(7^5), 3 + 3*7 + 6*7^2 + 4*7^3 + O(7^4), 6*7 + 5*7^2 + O(7^3), 4*7 + O(7^2), 5 + O(7))
                sage: xi = D((3 + 4*7 + 3*7^4 + O(7^5), 3*7^-1 + 3*7^2 + 3*7^3 + O(7^4), 6 + 6*7 + 7^2 + 4*7^3 + 6*7^4 + O(7^5), 5*7 + 6*7^3 + 4*7^4 + O(7^6), 6*7 + 3*7^3 + 6*7^4 + O(7^6)))
                sage: mu + xi
                7^-1 * (3*7^2 + 5*7^3 + 2*7^4 + O(7^5), 3 + 3*7 + 4*7^2 + 3*7^3 + O(7^4), 2*7 + 7^2 + O(7^3), 4*7 + O(7^2), O(7))
                sage: eta = D((7 + 5*7^2 + 6*7^3 + 4*7^4 + O(7^5), 2*7^3 + 5*7^5 + 4*7^6 + 5*7^7 + O(7^8), 7 + 6*7^2 + 2*7^4 + 3*7^5 + O(7^6), 7^2 + 2*7^3 + O(7^7), 2*7 + 4*7^2 + 4*7^3 + 2*7^4 + 5*7^5 + O(7^6)))
                sage: eta + nu
                (5 + 3*7 + 4*7^2 + 3*7^3 + 7^4 + O(7^5), 6*7 + 5*7^2 + 4*7^3 + O(7^4), 4 + 5*7 + 6*7^2 + O(7^3), 3 + 5*7 + O(7^2), 4 + O(7))
                sage: xi + eta
                7^-1 * (3*7 + 5*7^2 + 5*7^3 + 6*7^4 + O(7^5), 3 + 3*7^3 + O(7^4), 6*7 + O(7^3), O(7^2), O(7))
        """
        return self._addsub(right, False)
    
    def _sub_(self, right):
        #RH: "copied" from dist.pyx
        return self._addsub(right, True)
    
    def _addsub(self, right, negate):
        #RH: "copied" from dist.pyx
        verbose("\n******** begin _addsub ******** with negate: %s"%(negate), level=2)
        verbose("s_ordp %s, s_moments %s"%(self.ordp, self._moments), level=2)
        verbose("r_ordp %s, r_moments %s"%(right.ordp, right._moments), level=2)
        ans = self.parent()(0)
        aprec = min(self.ordp + len(self._moments), right.ordp + len(right._moments))
        ans.ordp = min(self.ordp, right.ordp)
        rprec = aprec - ans.ordp
        V = ans.parent().approx_module(rprec)
        R = V.base_ring()
        verbose([aprec, ans.ordp, rprec], level=2)
        smoments = copy(self._moments)
        rmoments = copy(right._moments)
        if smoments.parent() is not V:
            smoments = V(smoments.list(copy=False)[:rprec] + ([R(0)] * (rprec - len(smoments)) if rprec > len(smoments) else []))
        if rmoments.parent() is not V:
            rmoments = V(rmoments.list(copy=False)[:rprec] + ([R(0)] * (rprec - len(rmoments)) if rprec > len(rmoments) else []))
        # We multiply by the relative power of p
        if self.ordp > right.ordp:
            for i in range(rprec):
                smoments[i] = smoments[i] << (self.ordp - right.ordp)
        elif self.ordp < right.ordp:
            for i in range(rprec):
                rmoments[i] = rmoments[i] << (right.ordp - self.ordp)
        if negate:
            rmoments = -rmoments
        verbose(self._moments, level=2)
        ans._moments = smoments + rmoments
        verbose("ans_ordp %s, ans_moments %s"%(ans.ordp, ans._moments), level=2)
        verbose("\n******** end _addsub ********", level=2)
        return ans
    
    def _lmul_(self, right):
        """
        Scalar multiplication self*right.
        """
        #RH: mostly "copied" from dist.pyx but then changed
        ans = CoeffMod_OMS_element(None, self.parent(), None, False)
        #p = self.parent().prime()
        if right.is_zero():
            verbose("right is zero: %s"%(right), level=2)
            ans._moments = self.parent().approx_module(0)([])
            ans.ordp = min(self.parent().precision_cap(), right.valuation()+self.ordp)
        else:
            v, u = right.val_unit()
            ans._moments = self._moments * u
            ans.ordp = self.ordp + v
        return ans
    
    def _rmul_(self, left):
        #RH: "copied" from dist.pyx
        """
        Scalar multiplication left*self.
        """
        return self._lmul_(left)
    
    #def __eq__(self, other):
    #    parent_prec = self.parent()._prec_cap
    #    #print other, type(other)
    #    prec = min(len(self._moments), parent_prec, len(other._moments))    #this is a problem for checking bla == 0
    #    for i in range(prec):
    #        if self._moments[i].add_bigoh(prec - i) != other._moments[i].add_bigoh(prec - i):
    #            return False
    #    return True
    
    #def __nonzero__(self):
    #    """
    #    Checks that self is non-zero up to precision ...
    #    """
    #    parent_prec = self.parent()._prec_cap
    #    prec = min(len(self._moments), parent_prec)
    #    for i in range(prec):
    #        if self._moments[i].add_bigoh(prec - i) != 0:
    #            return True
    #    return False
    
    def _div_(self, right):
        r"""
            Divide ``self`` by ``right``.
            
        """
        #right = self.parent().base_ring()(right)    #is this necessary? it's probably already been coerced in
        if right.is_zero():
            raise ZeroDivisionError("Cannot divide by zero.")
        v, u = right.val_unit()
        return CoeffMod_OMS_element((~u) * self._moments, self.parent(), self.ordp - v, check=False)
    
    def __lshift__(self, shift):
        r"""
            The left shift operator <<. It "multiplies" self by the uniformizer
            raised to the shift power. In reality it simply adds shift to the
            valuation of self.
            
            EXAMPLES::
            
                sage: D = OverconvergentDistributions(3,base=ZpCA(3,5))
                sage: mu = D([-2,5,3,6,7,11]); mu
                (1 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + O(3^5), 2 + 3 + O(3^4), 3 + O(3^3), 2*3 + O(3^2), 1 + O(3))
                sage: mu << 3
                3^3 * (1 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + O(3^5), 2 + 3 + O(3^4), 3 + O(3^3), 2*3 + O(3^2), 1 + O(3))
                sage: mu << -2
                3^-2 * (1 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + O(3^5), 2 + 3 + O(3^4), 3 + O(3^3), 2*3 + O(3^2), 1 + O(3))
                sage: mu << 0
                (1 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + O(3^5), 2 + 3 + O(3^4), 3 + O(3^3), 2*3 + O(3^2), 1 + O(3))
                sage: nu = D(0)
                sage: nu = D(0); nu
                3^5 * ()
                sage: nu << 3   #Is this really the desired behaviour?
                3^8 * ()
        """
        ordp = self.ordp + shift
        return CoeffMod_OMS_element(self._moments, self.parent(), ordp, check=False)

    def __rshift__(self, shift):
        r"""
            The right shift operator >>. It "divides" self by the uniformizer
            raised to the shift power. In reality it simply subtracts shift from
            the valuation of self.
            
            EXAMPLES::
            
                sage: D = OverconvergentDistributions(3,base=ZpCA(3,5))
                sage: mu = D([-2,5,3,6,7,11]); mu
                (1 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + O(3^5), 2 + 3 + O(3^4), 3 + O(3^3), 2*3 + O(3^2), 1 + O(3))
                sage: mu >> -3
                3^3 * (1 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + O(3^5), 2 + 3 + O(3^4), 3 + O(3^3), 2*3 + O(3^2), 1 + O(3))
                sage: mu >> 2
                3^-2 * (1 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + O(3^5), 2 + 3 + O(3^4), 3 + O(3^3), 2*3 + O(3^2), 1 + O(3))
                sage: mu >> 0
                (1 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + O(3^5), 2 + 3 + O(3^4), 3 + O(3^3), 2*3 + O(3^2), 1 + O(3))
                sage: mu << 5 == 0
                False
        """
        ordp = self.ordp - shift
        return CoeffMod_OMS_element(self._moments, self.parent(), ordp, check=False)
        
    def __cmp__(left, right):
        r"""
        Compare left (i.e. self) and right, which are both elements of left.parent().
        
        EXAMPLES::
        
            sage: D = OverconvergentDistributions(0,base=Zp(3,5))
            sage: mu = D([1])
            sage: nu = D(0)
            sage: mu == nu
            False
            sage: mu = D([3,3,3,3,3,3])
            sage: nu = D([3,3,3])
            sage: mu == nu
            True
            sage: mu * 3^4 == mu << 4
            True
            sage: mu * 3^3 == 0
            False
        
        Check that something with no moments is 0::
        
            sage: from sage.modular.pollack_stevens.coeffmod_OMS_element import CoeffMod_OMS_element
            sage: R = ZpCA(7, 5)
            sage: D = OverconvergentDistributions(2, base=R)
            sage: mu = CoeffMod_OMS_element((R**0)([]), D, ordp=6, check=False)
            sage: mu == 0
            True
        
        Check that deepcopy is deep enough (the last two commands fail with copy
        instead of deepcopy)::
        
            sage: mu = CoeffMod_OMS_element((R**2)([R(2 * 7, 2), R(0, 1)]), D, ordp=-1, check=False)
            sage: nu = CoeffMod_OMS_element(D.approx_module(1)([R(2, 1)]), D, ordp=0, check=False)
            sage: mu == nu
            True
            sage: mu == nu
            True
            sage: mu
            (2 + O(7))
        """
        #RH: "copied" from dist.pyx, but then changed
        left = deepcopy(left)
        right = deepcopy(right)
        left.normalize()
        right.normalize()
        lrprec = left.precision_relative()
        rrprec = right.precision_relative()
        if lrprec == 0:
            if rrprec == 0:
                return 0
            else:
                return -1   # RH: This is probably incorrect
        elif rrprec == 0:
            return 1
        rprec = min(lrprec, rrprec)
        c = cmp(left.ordp, right.ordp)
        if c: return c
        return cmp(left._moments[:rprec], right._moments[:rprec])
#        p = left.parent().prime()
#        if left.ordp > right.ordp:
#            shift = p ** (left.ordp - right.ordp)
#            for i in range(rprec):
#                c = cmp(shift * left._unscaled_moment(i), right._unscaled_moment(i))
#                if c: return c
#        elif left.ordp < right.ordp:
#            shift = p ** (right.ordp - left.ordp)
#            for i in range(rprec):
#                c = cmp(left._unscaled_moment(i), shift * right._unscaled_moment(i))
#                if c: return c
#        else:
#            for i in range(rprec):
#                c = cmp(left._unscaled_moment(i), right._unscaled_moment(i))
#                if c: return c
#        return 0
        #        return cmp([left._moments,left.ordp],[right._moments,right.ordp])
    
    def is_zero(self, prec=None):
        #RH: Mostly "copied" from dist.pyx
        self.normalize()
        n = self.precision_relative()
        if n == 0:  #Should compare absolute prec first?
            return True
        aprec = self.precision_absolute()
        if prec is None:
            prec = n
        elif prec > aprec:
            return False    #Should this raise a PrecisionError instead
        elif prec < aprec:
            n -= (aprec - prec) #Why is this done?
            prec -= self.ordp
        p = self.parent().prime()
        usearg = True
        try:
            z = self._unscaled_moment(0).is_zero(prec)
        except TypeError:
            z = self._unscaled_moment(0).is_zero()
            use_arg = False
        if not z: return False
        for a in xrange(1, prec):
            if usearg:
                try:
                    z = self._unscaled_moment(a).is_zero(prec-a)
                except PrecisionError:
                    print "self_ordp, self._moments", self.ordp, self._moments
                    print "n, prec, a, aprec", n, prec, aprec, a
                    
                    assert(False)
                #z = self._unscaled_moment(a).is_zero(prec-a)
            else:
                z = self._unscaled_moment(a).is_zero()
            if not z: return False
        return True
    
    def find_scalar(self, other, M = None, check=True):
        i = 0
        n = self.precision_relative()
        other_pr = other.precision_relative()
        if n == 0:
            raise ValueError("self is zero")
## RP: This code doesn't seem right.  For instance, if the eigenvalue has positive valuation
##     then the relative precision will go down.
##        if n != other.precision_relative():
##            raise ValueError("other should have the same number of moments")
        verbose("n = %s"%n)
        verbose("moment 0")
        a = self._unscaled_moment(i)
        verbose("a = %s"%(a))
        p = self.parent().prime()
        v = a.valuation()
        while v >= n - i:
            i += 1
            verbose("p moment %s"%i)
            try:
                a = self._unscaled_moment(i)
            except IndexError:
                raise ValueError("self is zero")
            v = a.valuation()
        relprec = n - i - v
#            verbose("p=%s, n-i=%s\nself.moment=%s, other.moment=%s"%(p, n-i, a, other._unscaled_moment(i)),level=2)
## RP: This code was crashing because other may have too few moments -- so I added this bound with other's relative precision
        if i < other_pr:
            alpha = (other._unscaled_moment(i) / a).add_bigoh(n-i)
        else:
            alpha = (0*a).add_bigoh(other_pr-i)
        verbose("alpha = %s"%(alpha))
## RP: This code was crashing because other may have too few moments -- so I added this bound with other's relative precision
        while i < other_pr-1:
            i += 1
            verbose("comparing p moment %s"%i)
            a = self._unscaled_moment(i)
            if check:
#                   verbose("self.moment=%s, other.moment=%s"%(a, other._unscaled_moment(i)))
                if other._unscaled_moment(i) != alpha * a:
                    verbose("self._unscaled_moment=%s, other._unscaled_ moment=%s"%(a, other._unscaled_moment(i)))
                    raise ValueError("not a scalar multiple")
            v = a.valuation()
            if n - i - v > relprec:
                verbose("Reseting alpha: relprec=%s, n-i=%s, v=%s"%(relprec, n-i, v))
                relprec = n - i - v
                alpha = (other._unscaled_moment(i) / a).add_bigoh(n-i)
                verbose("alpha=%s"%(alpha))
        if relprec < M:
            raise ValueError("result not determined to high enough precision")
        alpha = alpha << (other.ordp - self.ordp)
        verbose("alpha=%s"%(alpha))
        try:
            return self.parent().base_ring()(alpha)
        except ValueError:
            return alpha
        
    def precision_relative(self):
        #RH: "copied" from dist.pyx
        return ZZ(len(self._moments))
    
    def precision_absolute(self):
        #RH: "copied" from dist.pyx
        return ZZ(len(self._moments) + self.ordp)
    
    def valuation(self):
        r"""
        Returns the `\varpi`-adic valuation of this distribution.
        
        .. WARNING::

            This function modifies the distribution in place since it calls
            :meth:`~sage.modular.pollack_stevens.coeffmod_OMS_element.CoeffMod_OMS_element.normalize`.
        
        EXAMPLES::
        
            sage: from sage.modular.pollack_stevens.coeffmod_OMS_element import CoeffMod_OMS_element
            sage: R = ZpCA(5, 5)
            sage: D = OverconvergentDistributions(2, base=R)
            sage: v = (R**5)((R(O(5^5)), R(O(5^5)), R(O(5^4)), R(5^2 + O(5^3)), R(4*5 + O(5^2))))
            sage: mu = CoeffMod_OMS_element(v, D, check=False)
            sage: mu.valuation()
            5
        """
        return self.normalize().ordp
    
    def _valuation(self, val_vector=False):
        r"""
        If self is normalized, this gives the same answer as
        :meth:`~sage.modular.pollack_stevens.coeffmod_OMS_element.CoeffMod_OMS_element.valuation`.
        Otherwise, it computes what is basically self.ordp plus the minimum of
        the \varpi`-adic valuations of the moments; however, `O(\varpi^a)` in the
        `i`th moment is considered to have valuation `a+i`. This is consistent
        with the identification
        
        .. MATH::
        
            (\varpi^rD^0)/\text{Fil}^M \cong \varpi^r\left(D^0/\text{Fil}^{M-r}\right).
        
        EXAMPLES::
        
            sage: from sage.modular.pollack_stevens.coeffmod_OMS_element import CoeffMod_OMS_element
            sage: R = ZpCA(5, 5)
            sage: D = OverconvergentDistributions(2, base=R)
            sage: v = (R**5)((R(O(5^5)), R(O(5^5)), R(O(5^4)), R(5^2 + O(5^3)), R(4*5 + O(5^2))))
            sage: mu = CoeffMod_OMS_element(v, D, check=False)
            sage: mu._valuation()
            1
            sage: mu._moments
            (O(5^5), O(5^5), O(5^4), 5^2 + O(5^3), 4*5 + O(5^2))
            sage: mu.normalize()
            5^5 * ()
            sage: mu._valuation()
            5
        """
        #RH: "copied" from dist.pyx, but then adjusted
        verbose("\n******** begin valuation ********", level=2)
        verbose("ordp %s, _moments %s"%(self.ordp, self._moments), level=2)
        n = self.precision_relative()
        if n == 0:
            if val_vector:
                return [self.ordp, []]
            return self.ordp
        if val_vector:
            cur = self._unscaled_moment(0).valuation()
            min_val = cur# if not cur.is_zero() else Infinity
            vv = [min_val]
            for a in range(1, n):
                cur_mom = self._unscaled_moment(a)
                cur = cur_mom.valuation() if not cur_mom.is_zero() else a + cur_mom.valuation()
                if cur < min_val:
                    min_val = cur
                vv.append(min_val)
            ret = self.ordp + min_val#min(n, min_val)
            verbose("ret %s"%(ret), level=2)
            verbose("\n******** end valuation ********", level=2)
            return [ret, vv]
        ret = self.ordp + min([self._unscaled_moment(a).valuation() if not self._unscaled_moment(a).is_zero() else a + self._unscaled_moment(a).valuation() for a in range(n)])#min([n] + [self._unscaled_moment(a).valuation() if not self._unscaled_moment(a).is_zero() else a + self._unscaled_moment(a).valuation() for a in range(n)])
        verbose("ret %s"%(ret), level=2)
        verbose("\n******** end valuation ********", level=2)
        return ret
    
    def normalize(self):
        r"""
            EXAMPLES::
            
                sage: D = OverconvergentDistributions(0,7,base=Qp(7,5))
                sage: mu = D([1,1,1]); mu
                (1 + O(7^3), 1 + O(7^2), 1 + O(7))
                sage: mu - mu
                7^3 * ()
                sage: D = OverconvergentDistributions(0,7,base=ZpCA(7,5))
                sage: mu = D([1,1,1]); mu
                (1 + O(7^3), 1 + O(7^2), 1 + O(7))
                sage: mu - mu
                7^3 * ()
        """
        #RH: "copied" from dist.pyx, but then changed (see issue #17)
        verbose("\n******** begin normalize ********", level=2)
        verbose("input: ordp %s, _moments %s"%(self.ordp, self._moments), level=2)
        V = self._moments.parent()
        R = V.base_ring()
        n = self.precision_relative()
        if n == 0:
            return self
        self_ordp = self.ordp
        self_val, val_vector = self._valuation(val_vector=True)
        while True:
            shift = self_val - self.ordp
            self.ordp = self_val
            #Factor out powers of uniformizer and check precision
            m = n
            adjust_moms = 0
            verbose("n: %s; shift: %s; _mom: %s"%(n, shift, self._moments), level=2)
            if shift > 0:
                for i in range(n):
                    self._moments[i] = self._moments[i] >> shift
                    adjust_moms = max(adjust_moms, m - self._moments[i].precision_absolute())
                    m -= 1
            elif shift == 0:
                for i in range(n):
                    adjust_moms = max(adjust_moms, m - self._moments[i].precision_absolute())
                    m -= 1
            else:
                raise NotImplementedError("Currently only deals with the case where the base ring is a ring of integers.")
            #Cut down moments because of precision loss
            verbose("adjust_moms: %s\nn %s\n_moms: %s"%(adjust_moms, n, self._moments), level=2)
            if adjust_moms >= n:
                V = self.parent().approx_module(0)
                self._moments = V([])
                #self.ordp = adjust_moms    #should we take min with parent().precision_cap()?
                verbose("adjust_moms %s, n %s, self.ordp %s"%(adjust_moms, n, self.ordp), level=2) 
                verbose("output: ordp %s, _moments %s"%(self.ordp, self._moments), level=2)
                verbose("\n******** end normalize ********", level=2)
                return self
            if adjust_moms == 0:
                for i in range(n):
                    self._moments[i] = self._moments[i].add_bigoh(n-i)
                verbose("adjust_moms %s, n %s, self.ordp %s"%(adjust_moms, n, self.ordp), level=2) 
                verbose("output: ordp %s, _moments %s"%(self.ordp, self._moments), level=2)
                verbose("\n******** end normalize ********", level=2)
                return self
            n -= adjust_moms
            val_diff = val_vector[n-1] - shift
            if val_diff == 0:
                V = self.parent().approx_module(n)
                self._moments = V([self._moments[i].add_bigoh(n - i) for i in range(n)])
                verbose("output: ordp %s, _moments %s"%(self.ordp, self._moments), level=2)
                verbose("\n******** end normalize ********", level=2)
                return self
            self_val = val_vector[n-1] + self_ordp
            val_vector = [val_vector[i] - shift for i in range(n)]
            self_ordp = self.ordp
            verbose("\nn: %s, self_val: %s, val_vector: %s, self.ordp: %s, self_ordp: %s"%(n, self_val, val_vector, self.ordp, self_ordp), level=2)
        #customized
    
    def moment(self, n):
        #RH: "copied" from dist.pyx
        r"""
        Returns the `n`-th moment.

        INPUT:

        - ``n`` -- an integer or slice, to be passed on to moments.

        OUTPUT:

        - the `n`-th moment, or a list of moments in the case that `n`
          is a slice.

        """
        if self.ordp == 0:
            return self._unscaled_moment(n)
        else:
            return self._unscaled_moment(n) << self.ordp
    
    def _unscaled_moment(self, n):
        #RH: "copied" from dist.pyx
        r"""
        Returns the `n`-th moment, unscaled by the overall power of p stored in self.ordp.
        """
        return self._moments[n]
    
    def reduce_precision(self, new_prec):
        #RH: adapted from dist.pyx
        if new_prec > self.precision_relative():
            raise ValueError("new_prec(=%s) must be less than relative precision of self."%(new_prec))
        moments = self._moments[:new_prec]
        ordp = self.ordp   #Should this be updated?
        return CoeffMod_OMS_element(moments, self.parent(), ordp, check=False)  #should latter be true?
    
    def reduce_precision_absolute(self, new_prec):
        #This isn't quite right, e.g. mu = (2*3^2 + 2*3^3 + O(3^4), 2*3 + O(3^3), O(3^2), 2 + O(3)) with new_prec=2
        if new_prec > self.precision_absolute():
            raise ValueError("new_prec(=%s) must be less than absolute precision of self."%(new_prec))
        ordp = self.ordp
        if new_prec - ordp <= 0:
            moments = self.parent().approx_module(0)([])
            ordp = new_prec
        else:
            moments = self._moments[:new_prec - ordp]
            moments[new_prec - ordp - 1] = moments[new_prec - ordp - 1].add_bigoh(1)
        return CoeffMod_OMS_element(moments, self.parent(), ordp, check=False)
    
    def lift(self, DD=None):
        #RH: should be good
        r"""
            This function lifts ``self`` into the space of families of overconvergent distributions ``DD``. If ``DD`` is None,
            it creates the space of families of overconvergent distributions of appropriate weight, character, etc. Otherwise,
            unlike simply coercing ``self`` into ``DD``, which would require compatibility of weights, actions, etc.,
            this function just brutally lifts the moments of ``self`` into ``DD``
        """
        if DD is None:
            DD = self.parent().lift()
            return DD(self)
        length = min(ZZ(len(self._moments)), DD.length_of_moments())
        prec = DD.length_reverse_lookup(length)
        actual_length = DD.length_of_moments(prec)
        ## actual_length might be one smaller than length because families of distributions can't achieve every possible length
        VV = DD.approx_module(prec)
        #V = self.parent().approx_module(length)
        new_moments = VV(self._moments[:actual_length])
        return DD.Element(new_moments, DD, ordp=self.ordp, check=False, var_prec=DD.precision_cap()[1])
    
    def solve_diff_eqn(self):
        r"""
            Solves the difference equation.
            
            See Theorem 4.5 and Lemma 4.4 of [PS].
            
            INPUT:
            
            - ``self`` - an overconvergent distribution `\mu` of absolute
              precision `M`
            
            OUTPUT:
            
            - an overconvergent distribution `\nu` of absolute precision
              `M - \lfloor\log_p(M)\rfloor - 1` such that
            
            .. math::
            
                \nu|\Delta = \mu,\text{ where }\Delta=\begin{pmatrix}1&1\\0&1
                \end{pmatrix} - 1.
            
            EXAMPLES::
            
                sage: D = OverconvergentDistributions(0,7,base=ZpCA(7,5))
                sage: D10 = D.change_precision(10)
                sage: mu10 = D10((O(7^10), 4 + 6*7 + 5*7^3 + 2*7^4 + 5*7^5 + O(7^9), 5 + 7^3 + 5*7^4 + 6*7^5 + 7^6 + 6*7^7 + O(7^8), 2 + 7 + 6*7^2 + 6*7^4 + 7^5 + 7^6 + O(7^7), 3*7 + 4*7^2 + 4*7^3 + 3*7^4 + 3*7^5 + O(7^6), 5 + 3*7 + 2*7^2 + 7^3 + 3*7^4 + O(7^5), 1 + 7^2 + 7^3 + O(7^4), 6*7 + 6*7^2 + O(7^3), 2 + 3*7 + O(7^2), 1 + O(7)))
                sage: nu10 = mu10.solve_diff_eqn()
                sage: MS = OverconvergentModularSymbols(14, coefficients=D)
                sage: MR = MS.source()
                sage: Id = MR.gens()[0]
                sage: nu10 * MR.gammas[Id] - nu10 - mu10
                7^8 * ()
                sage: D = OverconvergentDistributions(0,7,base=Qp(7,5))
                sage: D10 = D.change_precision(10)
                sage: mu10 = D10((O(7^10), 4 + 6*7 + 5*7^3 + 2*7^4 + 5*7^5 + O(7^9), 5 + 7^3 + 5*7^4 + 6*7^5 + 7^6 + 6*7^7 + O(7^8), 2 + 7 + 6*7^2 + 6*7^4 + 7^5 + 7^6 + O(7^7), 3*7 + 4*7^2 + 4*7^3 + 3*7^4 + 3*7^5 + O(7^6), 5 + 3*7 + 2*7^2 + 7^3 + 3*7^4 + O(7^5), 1 + 7^2 + 7^3 + O(7^4), 6*7 + 6*7^2 + O(7^3), 2 + 3*7 + O(7^2), 1 + O(7)))
                sage: nu10 = mu10.solve_diff_eqn()
                sage: MS = OverconvergentModularSymbols(14, coefficients=D);
                sage: MR = MS.source();
                sage: Id = MR.gens()[0]
                sage: nu10 * MR.gammas[Id] - nu10 - mu10
                7^8 * ()
                sage: R = ZpCA(5, 5); D = OverconvergentDistributions(0,base=R);
                sage: nu = D((R(O(5^5)), R(5 + 3*5^2 + 4*5^3 + O(5^4)), R(5 + O(5^3)), R(2*5 + O(5^2), 2 + O(5))));
                sage: nu.solve_diff_eqn()
                5 * (1 + 3*5 + O(5^2), O(5))
                
            Check input of relative precision 2::
            
                sage: from sage.modular.pollack_stevens.coeffmod_OMS_element import CoeffMod_OMS_element
                sage: R = ZpCA(3, 9)
                sage: D = OverconvergentDistributions(0, base=R, prec_cap=4)
                sage: V = D.approx_module(2)
                sage: nu = CoeffMod_OMS_element(V([R(0, 9), R(2*3^2 + 2*3^4 + 2*3^7 + 3^8 + O(3^9))]), D, ordp=0, check=False)
                sage: mu = nu.solve_diff_eqn()
                sage: mu
                3^2 * ()
        """
        #RH: see tests.sage for randomized verification that this function works correctly
        p = self.parent().prime()
        if self.is_zero():
            M = ZZ(self.precision_absolute())
            mu = self.parent()(0)
            mu.ordp = M - M.exact_log(p) - 1
            return mu
        if self._unscaled_moment(0) != 0:
            raise ValueError("Distribution must have total measure 0 to be in image of difference operator.")
        M = ZZ(len(self._moments))
        abs_prec = self.precision_absolute()
        ## RP: This should never happen -- the distribution must be 0 at this point if M==1
        if M == 1:
            return self.parent()(0)
        if M == 2:
            if p == 2:
                raise ValueError("Not enough accuracy to return anything")
            else:
                out_prec = abs_prec - abs_prec.exact_log(p) - 1
                if self.ordp >= out_prec:
                    mu = self.parent()(0)
                    mu.ordp = out_prec
                    return mu
                mu = self.parent()(self._unscaled_moment(1))
                mu.ordp = self.ordp
                return mu
        R = self.parent().base_ring()
        K = R.fraction_field()
        bern = [bernoulli(i) for i in range(0,M-1,2)]
        minhalf = ~K(-2)    #bernoulli(1)
        # bernoulli(1) = -1/2; the only nonzero odd bernoulli number
        v = [minhalf * self.moment(m) for m in range(M-1)] #(m choose m-1) * B_1 * mu[m]/m            
        for m in range(1,M):
            scalar = K(self.moment(m)) * (~K(m))
            for j in range(m-1,M-1,2):
                v[j] += binomial(j,m-1) * bern[(j-m+1)//2] * scalar
        ordp = min(a.valuation() for a in v)
        #Is this correct in ramified extensions of QQp?
        verbose("abs_prec: %s, ordp: %s"%(abs_prec, ordp), level=2)
        if ordp != 0:
            new_M = abs_prec - 1 - (abs_prec).exact_log(p) - ordp
            verbose("new_M: %s"%(new_M), level=2)
            V = self.parent().approx_module(new_M)
            v = V([R(v[i] >> ordp) for i in range(new_M)])
        else:
            new_M = abs_prec - 1 - (abs_prec).exact_log(p)
            verbose("new_M: %s"%(new_M), level=2)
            V = self.parent().approx_module(new_M)
            v = V([R(v[i]) for i in range(new_M)])
        v[new_M-1] = v[new_M-1].add_bigoh(1)  #To force normalize to deal with this properly
        mu = CoeffMod_OMS_element(v, self.parent(), ordp=ordp, check=False)
        verbose("mu.ordp: %s, mu._moments: %s"%(mu.ordp, mu._moments), level=2)
        return mu.normalize()

def test_correctness_and_precision_of_solve_diff_eqn(number=20, verbosity=1):
    """
    ``number`` is how many different random distributions to check. 
    
    Currently, avoids the prime 2.
    """
    from sage.misc.prandom import randint
    from sage.rings.arith import random_prime
    from sage.rings.padics.factory import ZpCA
    from sage.modular.pollack_stevens.coeffmod_OMS_space import OverconvergentDistributions
    from sage.structure.sage_object import dumps
    errors = []
    munus = []
    for i in range(number):
        Mspace = randint(1, 20)    #Moments of space
        M = randint(max(0, Mspace - 5), Mspace)
        p = random_prime(13, lbound=3)
        k = randint(0, 6)
        Rprec = Mspace + randint(0, 5)
        R = ZpCA(p, Rprec)
        D = OverconvergentDistributions(k, base=R, prec_cap=Mspace)
        S0 = D.action().actor()
        Delta_mat = S0([1,1,0,1])
        mu = D.random_element(M)
        mu_save = dumps(mu)#[deepcopy(mu.ordp), deepcopy(mu._moments)]
        if verbosity > 0:
            print "\nTest #{0} data (Mspace, M, p, k, Rprec) =".format(i+1), (Mspace, M, p, k, Rprec)
            print "mu =", mu
        
        nu = mu * Delta_mat - mu
        nu_save = [deepcopy(nu.ordp), deepcopy(nu._moments)]
        mu2 = nu.solve_diff_eqn()
        nu_abs_prec = nu.precision_absolute()
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
            print "     mured =", mu.reduce_precision_absolute(expected)
            print "       mu2 =", mu2
            print "        nu = ", nu
        if not agree:
            errors.append((i+1, 1))
            munus.append((mu_save, nu_save, mu2, (Mspace, M, p, k, Rprec)))
        if verbosity > 0:
            print "    Test finding mu from mu|Delta accurate: %s"%(agree)
            print "        nu_abs_prec  soln_abs_prec_expected  actual  agree"
        mu2_abs_prec = mu2.precision_absolute()
        agree = (expected == mu2_abs_prec)
        if verbosity > 0:
            print "        %s             %s                       %s      %s"%(nu_abs_prec, expected, mu2_abs_prec, agree)
        if not agree:
            errors.append((i+1, 2))
            munus.append((mu_save, nu_save, mu2, (Mspace, M, p, k, Rprec)))
        
        if mu.precision_relative() > 0:
            mu._moments[0] = R(0, mu.precision_relative())
        mu_save = [deepcopy(mu.ordp), deepcopy(mu._moments)]
        if verbosity > 0:
            print "    mu modified =", mu
        nu = mu.solve_diff_eqn()
        mu_abs_prec = mu.precision_absolute()
        expected = mu_abs_prec - mu_abs_prec.exact_log(p) - 1
        nud = nu * Delta_mat - nu
        nu_save = [deepcopy(nu.ordp), deepcopy(nu._moments)]
        agree = (nud - mu).is_zero(expected)
        if verbosity > 1:
            print "    Just so you know:"
            print "        mu =", mu
            print "     mured =", mu.reduce_precision_absolute(expected)
            print "       nud =", nud
        if not agree:
            errors.append((i+1, 3))
            munus.append((mu_save, nu_save, (Mspace, M, p, k, Rprec)))
        if verbosity > 0:
            print "    Test finding nu with nu|Delta == mu: %s"%(agree)
            print "        mu_abs_prec  soln_abs_prec_expected  actual  agree"
        nu_abs_prec = nu.precision_absolute()
        agree = (expected == nu_abs_prec)
        if verbosity > 0:
            print "        %s             %s                       %s      %s"%(mu_abs_prec, expected, nu_abs_prec, agree)
        if not agree:
            errors.append((i+1, 4))
            munus.append((mu_save, nu_save, (Mspace, M, p, k, Rprec)))
    if len(errors) == 0:
        if verbosity > 0:
            print "\nTest passed with no errors."
        return
    if verbosity > 0:
        print "\nTest failed with errors: %s\n"%(errors)
    return errors, munus
    
#def create__CoeffMod_OMS_element(moments, parent, ordp):
#    """
#    Used for unpickling.
#    """
#    return CoeffMod_OMS_element(moments, parent, ordp=ordp, check=False)
