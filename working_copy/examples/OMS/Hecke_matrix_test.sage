def test_Hecke_matrices_at_p(input_dict, rank=2, prec=5, sign=0, verbose=True):
    if verbose:
        print "\nVerifying rank %s Hecke polynomials at p (up to precision %s)"%(rank, prec)
        print "Constructing %s spaces of (sage) modular symbols."%(len(input_dict[rank]))
    MFs = {}
    for Npk in input_dict[rank]:
        N, p, k = Npk
        MF = ModularSymbols(N * p, k + 2, sign, Qp(p))
        if verbose:
            print MF
        MFs[Npk] = MF
    if verbose:
        print "Constructing bases of OMS spaces and verifying the Hecke polynomials."
    bases = {}
    polys = {}
    #problems = []
    #not_eigens = []
    count = 0
    for Npk in input_dict[rank]:
        if verbose:
            print "\nModular symbols space #%s"%(count+1)
            count += 1
            print Npk
        N, p, k = Npk
        #if k > 6 or p > 11:
        #    print "Passing this one."
        #    continue
        M = OverconvergentModularSymbols(N, k, p=p, prec_cap=prec, sign=sign)
        d = M.dimension_of_ordinary_subspace()
        if verbose:
            print "\nFinding basis (of dimension %s)."%(d)
        B = M.basis_of_ordinary_subspace(d)
        bases[Npk] = B
        if verbose:
            print "\nBasis:"
            for b in B:
                print b
            print "\nCompare Hecke polynomials."
        RP = PolynomialRing(M.base_ring().fraction_field(), 'x')
        if verbose:
            print "Find classical Hecke polynomial at %s:"%(p)
        hecke_poly_MF = ordinary_hecke_poly(RP(MFs[Npk].hecke_polynomial(p)))
        if verbose:
            print hecke_poly_MF
            print "\nFind OMS Hecke polynomial at %s:"%(p)
        hecke_poly_M = M.hecke_matrix(p, B).charpoly()
        if verbose:
            print hecke_poly_M
            print "\nAre they equal at {0}?".format(p), hecke_poly_M == hecke_poly_MF
        polys[Npk] = [hecke_poly_MF, hecke_poly_M]
    return bases, polys, MFs
    
def test_Hecke_matrices(input_dict, rank=2, max_ell=20, prec=5, sign=0, verbose=True):
    if verbose:
        print "\nVerifying rank %s Hecke polynomials up to ell = %s (up to precision %s)"%(rank, max_ell+1, prec)
        print "Constructing %s ordinary spaces of (sage) modular symbols."%(len(input_dict[rank]))
    MFs = {}
    for Npk in input_dict[rank]:
        N, p, k = Npk
        MF = ModularSymbols(N * p, k + 2, sign, Qp(p))
        if verbose:
            print MF
        MFs[Npk] = find_ordinary_subspace(MF)
    if verbose:
        print "Constructing bases of OMS spaces and verifying the Hecke polynomials."
    bases = {}
    polys = {}
    #problems = []
    #not_eigens = []
    count = 0
    for Npk in input_dict[rank]:
        if verbose:
            print "\nModular symbols space #%s"%(count+1)
            count += 1
            print Npk
        N, p, k = Npk
        #if k > 6 or p > 11:
        #    print "Passing this one."
        #    continue
        M = OverconvergentModularSymbols(N, k, p=p, prec_cap=prec, sign=sign)
        d = M.dimension_of_ordinary_subspace()
        if verbose:
            print "\nFinding basis (of dimension %s)."%(d)
        B = M.basis_of_ordinary_subspace(d)
        bases[Npk] = B
        if verbose:
            print "\nBasis:"
            for b in B:
                print b
            print "\nCompare Hecke polynomials."
        RP = M.base_ring()[x]
        for ell in prime_range(max_ell + 1):
            if verbose:
                print "Find classical Hecke polynomial at %s:"%(ell)
            hecke_poly_MF = RP(MFs[Npk].hecke_polynomial(ell))
            if verbose:
                print hecke_poly_MF
                print "\nFind OMS Hecke polynomial at %s:"%(ell)
            hecke_poly_M = M.hecke_matrix(ell, B).charpoly()
            if verbose:
                print hecke_poly_M
                print "\nAre they equal at {0}?".format(ell), hecke_poly_M == hecke_poly_MF
            polys[Npk] = [hecke_poly_MF, hecke_poly_M]
    return bases, polys, MFs

def find_ordinary_subspace(MF):
    p = MF.base_ring().prime()
    Up = MF.hecke_matrix(p)
    decomp = Up.decomposition()
    ordinary_Ws = []
    for W, _ in decomp:
        fW = Up.restrict(W).charpoly()
        slopes = fW.newton_slopes()
        if 0 in slopes:
            ordinary_Ws.append(W)
            if slopes != [0] * len(slopes):
                print "OH CRAP!", slopes
    return MF.submodule(sum(ordinary_Ws))

def ordinary_hecke_poly(f, verbose=True):
    """Given a Hecke polynomial (with coefficients in QQp), extract only the part with units roots."""
    if verbose:
        print "f =", f
    facts = f.factor()
    #print facts
    with_units = []
    for g, m in facts:
        if verbose:
            print g
        if g == f.parent().gen():
            continue
        slopes = g.newton_slopes()
        if verbose:
            print "slopes = ", slopes
        if 0 in slopes:
            with_units += [g] * m
            if slopes != [0] * len(slopes):
                print "OH CRAP!", slopes
    return prod(with_units)