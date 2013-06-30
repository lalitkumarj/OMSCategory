import sys

start_time = walltime()

#Initialization data
p = 3
N = 19
k = 2
r = (k-2) % (p-1)
M = 6  #Number of moments
var_prec = M    #Precision on the variable
max_iter = 50   #Maximum number of iterations of Up in projection to ordinary part

#print "This script computes the Hida family of X_0(19) at p = {0} with {1} moments and w-adic precision {2}.".format(p, M, var_prec)
#print "It verifies the computation by specializing to weight {0} and checking Hecke eigenvalues for primes up to {1}.".format(weights, max_ell)

DD = FamiliesOfOverconvergentDistributions(0, base_coeffs=ZpCA(p, M), prec_cap=[M,var_prec])
MM = FamiliesOfOMS(N, r, coefficients=DD, sign=-1)

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

val = Phis.valuation()
if val > 0:
    for key in Phis._map._dict.keys():
        Phis._map._dict[key].ordp -= val

print "Projecting to ordinary subspace."
sys.stdout.flush()
before = walltime()
for i in range(M+4):
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

print "Isolating connected component of rank 2."
sys.stdout.flush()
before = walltime()
for i in range(M+4):
    Phis = Phis.hecke(p) + Phis
print "Time elapsed:", walltime() - before

print "\nNormalizing the symbol again."
sys.stdout.flush()
before = walltime()
Phis.normalize()
print "Time elapsed:", walltime() - before

B = [Phis, Phis.hecke(p)]
HM = MM.hecke_polynomial_in_T_variable(p, basis=B)

#add more code once the above works...
