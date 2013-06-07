import sys

Ns = [11, 14, 15, 17, 19, 20, 21, 24]
k = 2
M = 4
var_prec = M

max_ell = 11
ells = prime_range(max_ell + 1)

for N in Ns:
    print N
    E = EllipticCurve('%sa'%N)
    a_ells_E = dict(zip(ells, [E.ap(ell) for ell in ells]))
    for p in N.prime_factors():
        if p == 2:
            continue
        print "   ", p
        sys.stdout.flush()
        r = (k - 2) % (p - 1)
        ZmodpM = Zmod(p ** (M+1))
        
        DD = FamiliesOfOverconvergentDistributions(0, base_coeffs=ZpCA(p, M), prec_cap=[M,var_prec])
        MM = FamiliesOfOMS(N, r, coefficients=DD, sign=-1)
        print "\nGenerating random modular symbol."
        sys.stdout.flush()
        Phis = MM.random_element()
        Phis.normalize();
        print "\nProjecting to ordinary subspace."
        sys.stdout.flush()
        for i in range(M + 2):
            Phis = Phis.hecke(p)
        Phis.normalize();
        
        a_ells = {}
        for ell in ells:
            a_ells[ell] = Phis.Tq_eigenvalue(ell)
            print "\na_{0} = {1}".format(ell, a_ells[ell])
            sys.stdout.flush()
        
        R = MM.base_ring().base_ring()
        for ell in ells:
            print "\na_{0}:".format(ell)
            known = R(a_ells_E[ell])
            ours = a_ells[ell][0]
            #if ell == p:
            #    ours = ours + p**(k-1) / ours
            print "Known: {0}".format(known)
            print "Us:    {0}".format(ours)
            print "Valuation of difference: {0}".format((ours - known).valuation())
            sys.stdout.flush()
        print "\n"
    print "\n*******************************"