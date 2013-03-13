from sage.misc.cachefunc import cached_function
from sage.structure.sequence import Sequence
from sage.rings.padics.factory import Qp
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.misc_c import prod
from sage.functions.other import factorial
from sage.matrix.constructor import Matrix

def ps_normalize(f, p, p_prec):
    """reduces all of the coefficients of the power series modulo p^N"""
    v = Sequence(f)
    v = [v[a] % (p ** p_prec) for a in range(len(v))]
    S = f.parent()
    return S(v)

def logp_fcn(p, p_prec, a):
    r"""
    INPUT:
    
    - ``p``: prime
    - ``p_prec: desired ``p``-adic precision
    - ``z``
    
    OUTPUT:
    
    - The ``p``-adic logarithm of 
    this is the *function* on Z_p^* which sends z to log_p(z) using a power series truncated at p_prec terms"""
    R = Qp(p, 2 * p_prec)
    a = a / R.teichmuller(a)
    return sum([((-1) ** (m - 1)) * ((a - 1) ** m) / m for m in range(1, p_prec)])

def logpp(p, p_prec):
    """returns the (integral) power series for log_p(1+p*z) -- extra p here!"""
    SS = PolynomialRing(QQ, 'y')
    y = SS.gen()
    return sum([((-1) ** (m - 1)) * ((p * y) ** m) / m for m in range(1, p_prec)])

def logpp_gam(p, p_prec):
    """returns the (integral) power series log_p(1+p*z)*(1/log_p(1+p)) where the denominator is computed with some accuracy"""
    L = logpp(p, p_prec)
    loggam = ZZ(logp_fcn(p, p_prec * (p ** 2), 1 + p))
    return ps_normalize(L / loggam, p, p_prec)

#@cached_function
def logpp_binom(n, p, p_prec):
    """returns the (integral) power series p^n*(log_p(1+p*z)/log_p(1+p) choose n)"""
    #prod=1+0*z
    L = logpp_gam(p, p_prec)
    ans = prod([(L - j) for j in range(n)])
    #for j in range(0,n):
    #    prod=prod*(L-j)
    ans *= (p ** n) / factorial(n)
    
    return ps_normalize(ans.truncate(p_prec), p, p_prec)

@cached_function
def automorphy_factor_vector(p, a, c, k, chi, p_prec, var_prec, R):
    """
    EXAMPLES::
        
        sage: from sage.modular.pollack_stevens.families_util import automorphy_factor_vector
        sage: automorphy_factor_vector(3, 1, 3, 0, None, 4, 3, PowerSeriesRing(ZpCA(3), 'w'))
        [1 + O(3^20), O(3^21) + (3 + 3^2 + 2*3^3 + O(3^21))*w + (3^2 + 2*3^3 + O(3^22))*w^2, O(3^22) + (3^2 + 2*3^3 + O(3^22))*w + (2*3^2 + O(3^22))*w^2, O(3^22) + (3^2 + 3^3 + O(3^22))*w + (2*3^3 + O(3^23))*w^2]
    """
    S = PolynomialRing(R, 'z')
    z = S.gens()[0]
    w = R.gen()
    aut = S(1)
    for n in range(1, var_prec):
        LB = logpp_binom(n, p, p_prec)
        ta = ZZ(Qp(p, 2 * max(p_prec, var_prec)).teichmuller(a))
        arg = (a / ta - 1) / p + c / (p * ta) * z
        aut += LB(arg) * (w ** n)
    aut *= (ta ** k)
    if not (chi is None):
        aut *= chi(a)
    aut = aut.list()
    len_aut = len(aut)
    if len_aut == p_prec:
        return aut
    elif len_aut > p_prec:
        return aut[:p_prec]
    return aut + [R.zero_element()] * (p_prec - len_aut)

def new_automorphy_factor_vector(p, a, c, k, chi, p_prec, var_prec, R):
    #if not R.is_capped_relative():
    #    Rcr =
    if p_prec > R.precision_cap():
        raise ValueError("Specified p-adic precision, p_prec, should be at most that of R.")
    S = PolynomialRing(R, 'z')
    z = S.gens()[0]
    w = R.gen()
    aut = S(R1, p_prec)
    ta = R.teichmuller(R(a, p_prec))
    
@cached_function
def automorphy_factor_matrix(p, a, c, k, chi, p_prec, var_prec, R):
    """
    EXAMPLES::
        
        sage: from sage.modular.pollack_stevens.families_util import automorphy_factor_matrix
        sage: automorphy_factor_matrix(3, 1, 3, 0, None, 4, 3, PowerSeriesRing(ZpCA(3), 'w'))
        [                                                          1 + O(3^20)                                                                     0                                                                     0                                                                     0]
        [O(3^20) + (3 + 3^2 + 2*3^3 + O(3^20))*w + (3^2 + 2*3^3 + O(3^20))*w^2                                                           1 + O(3^20)                                                                     0                                                                     0]
        [          O(3^20) + (3^2 + 2*3^3 + O(3^20))*w + (2*3^2 + O(3^20))*w^2 O(3^20) + (3 + 3^2 + 2*3^3 + O(3^20))*w + (3^2 + 2*3^3 + O(3^20))*w^2                                                           1 + O(3^20)                                                                     0]
        [            O(3^20) + (3^2 + 3^3 + O(3^20))*w + (2*3^3 + O(3^20))*w^2           O(3^20) + (3^2 + 2*3^3 + O(3^20))*w + (2*3^2 + O(3^20))*w^2 O(3^20) + (3 + 3^2 + 2*3^3 + O(3^20))*w + (3^2 + 2*3^3 + O(3^20))*w^2                                                           1 + O(3^20)]
    """
    #RH: there's a problem when p_prec = 1
    aut = automorphy_factor_vector(p, a, c, k, chi, p_prec, var_prec, R)
    M = Matrix(R, p_prec)
    for c in range(min(p_prec, len(aut))):
        for r in range(min(p_prec, len(aut) - c)):
            M[r + c, c] = aut[r]
    return M