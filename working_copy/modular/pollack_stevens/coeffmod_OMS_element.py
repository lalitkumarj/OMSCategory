# This should be cythoned once it's done.

from copy import copy   #Not necessary in cython version
from sage.misc.misc import verbose
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.arith import binomial, bernoulli
from sage.modular.pollack_stevens.coeffmod_element import CoefficientModuleElement_generic
from sage.rings.integer_ring import ZZ

maxordp = 2 ** 14   #change this to what it is in dist.pyx

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
                if len(moments) == 0 or moments == 0:   #should do something with how "zero" moments is
                    V = parent.approx_module(0)
                    moments = V([])
                    ordp = parent.precision_cap()
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
                moments = V(moments >> ordp)
        
        self._moments = moments
        #if var_prec is None:
        #    self._var_prec = parent._prec_cap[1]
        #else:
        #    self._var_prec = var_prec
        self.ordp = ordp
    
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
        return self._addsub(right, False)
    
    def _sub_(self, right):
        #RH: "copied" from dist.pyx
        return self._addsub(right, True)
    
    def _addsub(self, right, negate):
        #RH: "copied" from dist.pyx
        ans = self.parent()(0)
        aprec = min(self.ordp + len(self._moments), right.ordp + len(right._moments))
        ans.ordp = min(self.ordp, right.ordp)
        rprec = aprec - ans.ordp
        V = ans.parent().approx_module(rprec)
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
        """
        Scalar multiplication self*right.
        """
        #RH: mostly "copied" from dist.pyx but then changed
        ans = CoeffMod_OMS_element(None, self.parent(), None, False)
        #p = self.parent().prime()
        if right.is_zero():
            verbose("right is zero: %s"%(right))
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
            
            EXAMPLES::
            
                
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
                True
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
                sage: mu * 3^4 == 0
                True
                sage: mu * 3^3 == 0
                False
        """
        #RH: "copied" from dist.pyx, but then changed
        left = copy(left)
        right = copy(right)
        left.normalize()
        right.normalize()
        rprec = min(left.precision_relative(), right.precision_relative())
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
        v = a.valuation(p)
        while v >= n - i:
            i += 1
            verbose("p moment %s"%i)
            try:
                a = self._unscaled_moment(i)
            except IndexError:
                raise ValueError("self is zero")
            v = a.valuation(p)
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
            v = a.valuation(p)
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
        #RH: "copied" from dist.pyx, but then adjusted
        p = self.parent().prime()
        n = self.precision_relative()
        if n == 0:
            return self.ordp
        return self.ordp + min([self.parent().precision_cap()] + [self._unscaled_moment(a).valuation(p) for a in range(n) if not self._unscaled_moment(a).is_zero()])
    
    def normalize(self):
        #RH: "copied" from dist.pyx, but then changed (see issue #17)
        V = self._moments.parent()
        R = V.base_ring()
        n = self.precision_relative()
        self_val = self.valuation()
        shift = self_val - self.ordp
        self.ordp = self_val
        #Factor out powers of uniformizer and check precision
        m = n
        adjust_moms = 0
        verbose("n: %s; shift: %s; _mom: %s"%(n, shift, self._moments))
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
        verbose("adjust_mom: %s\n_moms: %s"%(adjust_moms, self._moments))
        if adjust_moms >= n:
            V = self.parent().approx_module(0)
            self._moments = V([])
            self._ordp = adjust_moms    #should we take min with parent().precision_cap()?
        elif adjust_moms > 0:
            n -= adjust_moms
            V = self.parent().approx_module(n)
            self._moments = V([self._moments[i].add_bigoh(n - i) for i in range(n)])
        else:
            for i in range(n):
                self._moments[i] = self._moments[i].add_bigoh(n-i)
        return self
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
            return self.parent().prime()**(self.ordp) * self._unscaled_moment(n)
    
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
        length = min(ZZ(len(self._moments)), DD.precision_cap()[0])
        VV = DD.approx_module(length)
        V = self.parent().approx_module(length)
        new_moments = VV(self._moments[:length])
        return DD.Element(new_moments, DD, ordp=self.ordp, check=False)
    
    def solve_diff_eqn(self):
        r"""
            Solves the difference equation.
            
            See Theorem 4.5 and Lemma 4.4 of [PS].
            
            OUTPUT:
            
            - a distribution v so that self = v | Delta, where Delta = [1, 1; 0, 1] - 1.
            
            EXAMPLES::
            
                sage: D = OverconvergentDistributions(0,7,base=ZpCA(7,5))
                sage: D10 = D.change_precision(10)
                sage: mu10 = D10((O(7^10), 4 + 6*7 + 5*7^3 + 2*7^4 + 5*7^5 + O(7^9), 5 + 7^3 + 5*7^4 + 6*7^5 + 7^6 + 6*7^7 + O(7^8), 2 + 7 + 6*7^2 + 6*7^4 + 7^5 + 7^6 + O(7^7), 3*7 + 4*7^2 + 4*7^3 + 3*7^4 + 3*7^5 + O(7^6), 5 + 3*7 + 2*7^2 + 7^3 + 3*7^4 + O(7^5), 1 + 7^2 + 7^3 + O(7^4), 6*7 + 6*7^2 + O(7^3), 2 + 3*7 + O(7^2), 1 + O(7)))
                sage: nu10 = mu10.solve_diff_eqn()
                sage: set_verbose(2)
                sage: MS = OverconvergentModularSymbols(14, coefficients=D);
                sage: MR = MS.source();
                sage: Id = MR.gens()[0]
                sage: nu10 * MR.gammas[Id] - nu10 - mu10
                7^9 * ()
                sage: D = OverconvergentDistributions(0,7,base=Qp(7,5))
                sage: D10 = D.change_precision(10)
                sage: mu10 = D10((O(7^10), 4 + 6*7 + 5*7^3 + 2*7^4 + 5*7^5 + O(7^9), 5 + 7^3 + 5*7^4 + 6*7^5 + 7^6 + 6*7^7 + O(7^8), 2 + 7 + 6*7^2 + 6*7^4 + 7^5 + 7^6 + O(7^7), 3*7 + 4*7^2 + 4*7^3 + 3*7^4 + 3*7^5 + O(7^6), 5 + 3*7 + 2*7^2 + 7^3 + 3*7^4 + O(7^5), 1 + 7^2 + 7^3 + O(7^4), 6*7 + 6*7^2 + O(7^3), 2 + 3*7 + O(7^2), 1 + O(7)))
                sage: nu10 = mu10.solve_diff_eqn()
                sage: MS = OverconvergentModularSymbols(14, coefficients=D);
                sage: MR = MS.source();
                sage: Id = MR.gens()[0]
                sage: nu10 * MR.gammas[Id] - nu10 - mu10
                7^9 * ()
        """
        #RH: see tests.sage for randomized verification that this function works correctly
        # assert self._moments[0][0]==0, "not total measure zero"
        # print "result accurate modulo p^",self.moment(0).valuation(self.p)
        #v=[0 for j in range(0,i)]+[binomial(j,i)*bernoulli(j-i) for j in range(i,M)]
#        if self.is_zero():
#            return self.parent()(0)
#        if self._unscaled_moment(0) != 0:
#            raise ValueError("Distribution must have total measure 0 to be in image of difference operator.")
#        M = len(self._moments)
#        if M == 1:
#            return self.parent()(0)
#        if M == 2:
#            return self.parent()(self.moment(1))
#        R = self.parent().base_ring()
#        p = self.parent().prime()
#        pm = p ** self.ordp
#        shift = 0
#        while pm <= M:   #looking in F_k(M-1)
#            pm *= p
#            shift += 1
#        self.ordp += shift
#        #K = R.fraction_field()
#        V = self.parent().approx_module(M-1)
#        bern = [bernoulli(i) for i in range(0,M-1,2)]
#        #What about p=2
#        minhalf = ~(-2)    #bernoulli(1)
#        # bernoulli(1) = -1/2; the only nonzero odd bernoulli number
#        v = [minhalf * self.moment(m) for m in range(M-1)] #(m choose m-1) * B_1 * mu[m]/m            
#        for m in range(1,M):
#            scalar = self.moment(m) / m
#            for j in range(m-1,M-1,2):
#                v[j] += binomial(j,m-1) * bern[(j-m+1)//2] * scalar
#        ordp = min(a.valuation(p) for a in v)
#        #Is this correct in ramified extensions of QQp?
#        if ordp < 0:
#            print "This should never happen!"
#            scalar = K(p) ** (-ordp)
#            v = V([R(a * scalar) for a in v])
#        elif ordp > 0:
#            scalar = p ** ordp
#            v = V([R(a // scalar) for a in v])
#        else:
#            v = V([R(a) for a in v])
#        return CoeffMod_OMS_element(v, self.parent(), ordp=ordp-shift, check=False)
        p = self.parent().prime()
        if self.is_zero():
            M = ZZ(self.precision_absolute())
            mu = self.parent()(0)
            mu.ordp = M-M.exact_log(p)-1
            return mu
        if self._unscaled_moment(0) != 0:
            raise ValueError("Distribution must have total measure 0 to be in image of difference operator.")
        M = ZZ(len(self._moments))
        ## RP: This if should never happen -- the distribution must be 0 at this point if M==1
        if M == 1:
            return self.parent()(0)
        if M == 2:
            if p == 2:
                raise ValueError("Not enough accuracy to return anything")
            else:
                mu = self.parent()(self._unscaled_moment(1))
                mu.ordp = self.ordp
                return mu
        R = self.parent().base_ring()
        K = R.fraction_field()
        V = self.parent().approx_module(M-1)
        bern = [bernoulli(i) for i in range(0,M-1,2)]
        minhalf = ~K(-2)    #bernoulli(1)
        # bernoulli(1) = -1/2; the only nonzero odd bernoulli number
        v = [minhalf * self.moment(m) for m in range(M-1)] #(m choose m-1) * B_1 * mu[m]/m            
        for m in range(1,M):
            scalar = K(self.moment(m)) * (~K(m))
            for j in range(m-1,M-1,2):
                v[j] += binomial(j,m-1) * bern[(j-m+1)//2] * scalar
        ordp = min(a.valuation(p) for a in v)
        #Is this correct in ramified extensions of QQp?
        if ordp < 0:
            scalar = K(p) ** (-ordp)
            v = V([R(a * scalar) for a in v])
        elif ordp > 0:
            scalar = K(p) ** ordp
            v = V([R(a // scalar) for a in v])
        else:
            v = V([R(a) for a in v])
        mu = CoeffMod_OMS_element(v, self.parent(), ordp=ordp, check=False)
        e = self.ordp - mu.ordp  ## this is the amount the valuation dropped
        f = M.exact_log(p)        ## this is the amount we need it to drop
        mu = mu.reduce_precision(mu.precision_relative()-(f-e))
        return mu.normalize()