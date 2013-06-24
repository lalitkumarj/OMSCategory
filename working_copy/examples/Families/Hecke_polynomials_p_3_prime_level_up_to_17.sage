import sys

max_ell = 19
ells = prime_range(max_ell + 1)
p = 3
k = 2
r = k % (p-1)
sign = -1
M = 6
var_prec = M

start = walltime()
HPs = {}
discs = {}
for ell in ells:
    MM = FamiliesOfOMS(ell, r, sign=sign, p=p, prec_cap=[M, var_prec], base_coeffs=ZpCA(p, M))
    print "Level %s has dimension %s."%(ell, MM.dimension_of_ordinary_subspace(cusp=True))
    sys.stdout.flush()
    HP = MM.hecke_polynomial_in_T_variable(p, verbose=False)
    print "\nChar. poly. of Up:"
    print HP
    sys.stdout.flush()
    HPs[(p, ell)] = HP
    disc = HP.discriminant()
    print "\nIts discriminant:"
    print disc
    print "-------------------\n"
    sys.stdout.flush()
    discs[(p, ell)] = disc

print "Time elapse:", walltime() - start