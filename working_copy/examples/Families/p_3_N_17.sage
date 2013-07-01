import sys
from sage.modular.pollack_stevens.modsym_OMS_families_space import Iwasawa_invariants

start_time = walltime()

#Initialization data
p = 3
N = 17
k = 2
r = (k-2) % (p-1)
M = 10  #Number of moments
var_prec = M    #Precision on the variable
max_iter = 50   #Maximum number of iterations of Up in projection to ordinary part

#print "This script computes the Hida family of X_0(19) at p = {0} with {1} moments and w-adic precision {2}.".format(p, M, var_prec)
#print "It verifies the computation by specializing to weight {0} and checking Hecke eigenvalues for primes up to {1}.".format(weights, max_ell)

DD = FamiliesOfOverconvergentDistributions(0, base_coeffs=ZpCA(p, M), prec_cap=[M,var_prec])
MM = FamiliesOfOMS(N, r, coefficients=DD, sign=1)

print "This example looks at p=3 and N=17"
print "In weight 2, there are 3 ordinary forms (all newforms).  One is defined over Q_3 and the other two are conjugate over Q_9."
print "The Hida algebra should then be Lambda \oplus (Lambda \otimes Z_9)."

print "\nIn what follows, we compute the characteristic polynomial of U_3 projected to the two-dimensional local piece of the Hecke algebra."
print "We are using precision: %s"%([M,var_prec])

B = []

while len(B) < 2:
    print "\nGenerating a random modular symbol."
    sys.stdout.flush()
    before = walltime()
    Phis = MM.random_element()
    print "Time elapsed:", walltime() - before

    Phis.normalize()

    val = Phis.valuation()
    if val > 0:
        for key in Phis._map._dict.keys():
            Phis._map._dict[key].ordp -= val

    print "Projecting to ordinary subspace."
    sys.stdout.flush()
    before = walltime()
    for i in range(M+4):
        #print "    Iteration %s of Up"%(i+1)
        sys.stdout.flush()
        Phis = Phis.hecke(p)
    print "Time elapsed:", walltime() - before

    print "Positive valuation after projection.  Scaling away"
    val = Phis.valuation()
    if val > 0:
        for key in Phis._map._dict.keys():
            Phis._map._dict[key].ordp -= val

    print "Isolating connected component of rank 2."
    sys.stdout.flush()

    print "  Killing off Eisenstein series and newform defined over Z_3 (applying U_3 - 1)"
    sys.stdout.flush()
    before = walltime()
    for i in range(M+4):
        Phis = Phis.hecke(p) - Phis
    print "Time elapsed:", walltime() - before

    print "Positive valuation after projection.  Scaling away"
    val = Phis.valuation()
    if val > 0:
        for key in Phis._map._dict.keys():
            Phis._map._dict[key].ordp -= val

    Phis.normalize()

    if MM.is_start_of_basis(B + [Phis]):
        B = B + [Phis]
        print "Adding this new symbol to our basis"
    else:
        print "FAILED. Symbol combined with previous basis does not generate free module"

print "\nComputing U_3"
sys.stdout.flush()
before = walltime()
T3_poly = MM.hecke_polynomial_in_T_variable(3,basis=B,verbose=False)
b = T3_poly.padded_list()[1]
c = T3_poly.padded_list()[0]
print "Time elapsed:", walltime() - before

print "\nThe characteristic polynomial of U_3 is: %s"%(T3_poly)
print "The discriminant of char poly of U_3 is: %s"%(b^2-4*c)
print "The Iwasawa invariants (mu,lambda) of the discriminant are",(Iwasawa_invariants(b^2-4*c))


for q in primes(10):
    if q != p:
        print "\nComputing T_%s"%(q)
        sys.stdout.flush()
        before = walltime()
        Tq_poly = MM.hecke_polynomial_in_T_variable(q,basis=B,verbose=False)
        b = Tq_poly.padded_list()[1]
        c = Tq_poly.padded_list()[0]
        print "Time elapsed:", walltime() - before

        print "The characteristic polynomial of T_%s is: %s"%(q,Tq_poly)
        print "The discriminant of char poly of T_%s is: %s"%(q,b^2-4*c)
        print "The Iwasawa invariants (mu,lambda) of the discriminant are",Iwasawa_invariants(b^2-4*c)
        print "Time elapsed:", walltime() - before
