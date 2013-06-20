import sys

start_time = walltime()

#Initialization data
p = 11
k = 2
r = (k-2) % (p-1)
M = 8  #Number of moments
var_prec = M    #Precision on the variable
max_iter = 50   #Maximum number of iterations of Up in projection to ordinary part

#Verification data
max_ell = 20
ells = prime_range(max_ell + 1)
E = EllipticCurve('11a1')
a_ells_E = dict(zip(ells, [E.ap(ell) for ell in ells]))
f = Newforms(1, 12)[0]
a_ells_f = dict(zip(ells, [f[ell] for ell in ells]))
weights = (2, 12)
a_ells_dict = dict(zip(weights, [a_ells_E, a_ells_f]))
unstabilize = {2: False, 12: True}
ZmodpM = Zmod(p ** (M+1))
weights_in_w = dict(zip(weights, [(ZmodpM(1+p) ** (wt-2) - 1).lift() / p for wt in weights]))

print "This script computes the Hida family of X_0(11) at p = {0} with {1} moments and w-adic precision {2}.".format(p, M, var_prec)
print "It verifies the computation by specializing to weights {0} and checking Hecke eigenvalues for primes up to {1}.".format(weights, max_ell)

DD = FamiliesOfOverconvergentDistributions(0, base_coeffs=ZpCA(p, M), prec_cap=[M,var_prec])
MM = FamiliesOfOMS(p, r, coefficients=DD)

print "\nGenerating random modular symbol."
sys.stdout.flush()
before = walltime()
Phis = MM.random_element()
print "Time elapsed:", walltime() - before

print "\nNormalizing the symbol."
sys.stdout.flush()
before = walltime()
Phis.normalize()
print "Time elapsed:", walltime() - before

print "Projecting to ordinary subspace."
sys.stdout.flush()
before = walltime()
for i in range(M+2):
    print "    Iteration %s of Up"%(i+1)
    sys.stdout.flush()
    Phis = Phis.hecke(p)
#i = 0
#try:
#    is_ord = Phis.is_ordinary()
#except:
#    is_ord = False
#while i < max_iter and not is_ord:
#    print "    Iteration %s of Up"%(i+1)
#    sys.stdout.flush()
#    Phis = Phis.hecke(p)
#    if (i+1) % (M+2) == 0:
#        try:
#            is_ord = Phis.is_ordinary()
#        except:
#            pass
#    i += 1
print "Time elapsed:", walltime() - before

print "Killing Eisenstein component (the dumb way)."
sys.stdout.flush()
before = walltime()
Phis = Phis.hecke(p) - Phis
print "Time elapsed:", walltime() - before

print "\nNormalizing the symbol again."
sys.stdout.flush()
before = walltime()
Phis.normalize()
print "Time elapsed:", walltime() - before

print "\nCompute the Hecke eigenvalues:"
sys.stdout.flush()
before = walltime()
a_ells = {}
for ell in ells:
    a_ells[ell] = Phis.Tq_eigenvalue(ell)
    print "\na_{0} = {1}".format(ell, a_ells[ell])
    sys.stdout.flush()
print "Time elapsed:", walltime() - before

print "\nCompare to known weights {0}".format(weights)

R = MM.base_ring().base_ring()
for wt in weights:
    print "\nWeight {0}:".format(wt)
    for ell in ells:
        print "\na_{0}:".format(ell)
        known = R(a_ells_dict[wt][ell])
        ours = a_ells[ell].truncate()(weights_in_w[wt])
        if ell == p and unstabilize[wt]:
            ours = ours + p**(wt-1) / ours
        print "Known: {0}".format(known)
        print "Us:    {0}".format(ours)
        print "Valuation of difference: {0}".format((ours - known).valuation())
        sys.stdout.flush()

print "\nTotal time elapsed:", walltime() - start_time