def compute_dims_of_ord_subspace(Nbound=30, ps=[3,5,7,11,13], ks=None, input_dict={}, cusp=False, sign=-1):
    if len(input_dict) > 0:
        raise NotImplementedError
    else:
        Nstart = 1
    for N in range(Nstart, Nbound+1):
        verbose("Level %s"%(N))
        for p in ps:
            if N % p == 0:
                continue
            verbose("    Prime %s"%(p))
            for k in range(p-1 if ks is None else ks[p]):
                verbose("    Computing dim for weight %s"%(k))
                M = OverconvergentModularSymbols(N, k, p=p, prec_cap=5, sign=sign)
                d = M.dimension_of_ordinary_subspace(cusp=cusp)
                verbose("    Got: %s"%(d))
                if d not in input_dict:
                    input_dict[d] = []
                    verbose("        This is a new record!")
                input_dict[d].append((N, p, k))
    return input_dict

def verify_rank_one(input_dict, max_ell=20, prec=10, verbose=True):
    """
        First run compute_dims_of_cusp_ord_subspace(), then pass the output as
        the input_dict of this function.
        Or instead, run this command in sage
            dims_dict = load("working_copy/examples/data/dims_of_ord_cusp.sobj")
        and pass this as input_dict.
    """
    if verbose:
        print "\nVerifying rank one up to precision %s"%(prec)
        print "Constructing %s newforms."%(len(input_dict[1]))
    NFs = {}
    for Npk in input_dict[1]:
        N, p, k = Npk
        NFs[Npk] = Newforms(N, k + 2, names='a') if k != 0 else Newforms(N * p, k + 2, names='a')
    if verbose:
        print "Constructing OMSs and verifying they match up."
    Phis = {}
    problems = []
    not_eigens = []
    count = 0
    for Npk in input_dict[1]:
        if verbose:
            print "\nNewform #%s"%(count+1)
            count += 1
            print Npk
        N, p, k = Npk
        #if k > 6 or p > 11:
        #    print "Passing this one."
        #    continue
        f = NFs[Npk][0]
        if verbose:
            print "Finding non-Eisenstein prime"
        eis_ell = 2
        eis = True
        while eis:
            if (f[eis_ell] == 1 + eis_ell ** (k + 1)) or ((N*p) % eis_ell == 0):
                eis_ell = next_prime(eis_ell)
            else:
                eis = False
                if verbose:
                    print "Non-Eisenstein prime =", eis_ell
        M = OverconvergentModularSymbols(N, k, p=p, prec_cap=prec)
        if prec is None:
            prec = M.precision_cap()
        Phi = M.random_element()
        is_ord = False
        j = 0
        while not is_ord and j < 15:
            for i in range(prec + 2):
                if verbose:
                    iter_nb = j * (prec + 2) + i + 1
                    print "Iteration {0} of Up.".format(iter_nb)
                Phi = Phi.hecke(p)
            if verbose:
                print "Killing Eisenstein component with ell = {0} (eis={1}).".format(eis_ell, eis)
            Phi = Phi.hecke(eis_ell) - (1 + eis_ell ** (k + 1)) * Phi
            try:
                is_ord = Phi.is_ordinary()
                if verbose:
                    print "Up eigenvalue:", Phi.Tq_eigenvalue(p)
            except NotImplementedError:
                pass
            j += 1
        Phis[Npk] = Phi
        if verbose:
            print "Working with an eigensymbol with precision %s and valuation %s"%(Phi.precision_absolute(),Phi.valuation())       
        R = M.base_ring()
        for ell in prime_range(max_ell + 1):
            a_ell_f = R(f[ell], prec)
            try:
                a_ell_Phi = Phi.Tq_eigenvalue(ell)
            except ValueError:
                if verbose:
                    print "                 Further iterating Up to get a T{0} eigensymbol".format(ell)
                not_eigen = True
                j = 0
                while not_eigen and j < 15:
                    for i in range(prec + 2):
                        if verbose:
                            iter_nb += 1
                            print "Iteration {0} of Up.".format(iter_nb)
                        Phi = Phi.hecke(p)
                    if verbose:
                        print "Killing Eisenstein component with ell = {0} (eis={1}).".format(eis_ell, eis)
                    Phi = Phi.hecke(eis_ell) - (1 + eis_ell ** (k + 1)) * Phi  #Should we kill again?
                    try:
                        a_ell_Phi = Phi.Tq_eigenvalue(ell)
                        not_eigen = False
                    except ValueError:
                        pass
                    j += 1
                if not_eigen:
                    not_eigens.append([Npk, Phi])
                    continue
                if verbose:
                    print "Working with an eigensymbol with precision %s and valuation %s"%(Phi.precision_absolute(),Phi.valuation())
            if k != 0 and ell == p:
                #un-p-stablize
                a_ell_Phi = a_ell_Phi + (p ** (k + 1)) / a_ell_Phi
            checker = a_ell_f == a_ell_Phi
            if verbose:
                print "a_%s: %s ?= %s"%(ell, a_ell_f, a_ell_Phi)
                print checker, "up to precision", (a_ell_f - a_ell_Phi).valuation()
            if not checker:
                problems.append([Npk, Phi])
                if verbose:
                    print "PROBLEM at {0}".format(Npk)
    return Phis, NFs, problems, not_eigens
