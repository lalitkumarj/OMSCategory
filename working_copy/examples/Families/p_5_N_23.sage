import sys
from sage.modular.pollack_stevens.modsym_OMS_families_space import Iwasawa_invariants

start_time = walltime()

#Initialization data
p = 5
N = 23
k = 2
r = (k-2) % (p-1)
M = 10  #Number of moments
var_prec = M    #Precision on the variable
max_iter = 50   #Maximum number of iterations of Up in projection to ordinary part

#print "This script computes the Hida family of X_0(19) at p = {0} with {1} moments and w-adic precision {2}.".format(p, M, var_prec)
#print "It verifies the computation by specializing to weight {0} and checking Hecke eigenvalues for primes up to {1}.".format(weights, max_ell)

DD = FamiliesOfOverconvergentDistributions(0, base_coeffs=ZpCA(p, M), prec_cap=[M,var_prec])
MM = FamiliesOfOMS(N, r, coefficients=DD, sign=-1)

print "This example looks at p=5 and N=23"
print "In weight 2, there are 7 ordinary forms all new: 1 is defined over Z_5, 2 are conjugate over Q(sqrt(5)) (and thus are congruent), 4 are conjugate over a quartic extension K/Q in which 5 is inert."
print "Moreover, there are no other congruences.  Thus, the Hida algebra has 4 local pieces: T = T_1 + T_2 + T_4 as ordered above."
print "We know that T_1 is Lambda and T_4 is Lambda \otimes O where O is a quartic unramified extension of Z_5."
print "\nIn what follows, we compute the characteristic polynomial of various Hecke operators projected to T_2."
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
    for i in range(2*M+4):
        #print "    Iteration %s of Up"%(i+1)
        sys.stdout.flush()
        Phis = Phis.hecke(p)
    print "Time elapsed:", walltime() - before

    print "Isolating T_1"
    sys.stdout.flush()

    print "  Killing off rank 1 piece (applying T_2 - 2)"
    sys.stdout.flush()
    before = walltime()
    for i in range(2*M+4):
        Phis = Phis.hecke(2) - 2 * Phis
    print "Time elapsed:", walltime() - before

    Phis.normalize()

    print "  Killing off rank 4 piece (applying U_5 - 1)"
    sys.stdout.flush()
    before = walltime()
    for i in range(2*M+4):
        Phis = Phis.hecke(5) - Phis 
    print "Time elapsed:", walltime() - before

    Phis.normalize()

    if MM.is_start_of_basis(B + [Phis]):
        B = B + [Phis]
        print "Adding this new symbol to our basis"
    else:
        print "FAILED. Symbol combined with previous basis does not generate free module"

print "\nComputing U_%s"%(p)
sys.stdout.flush()
before = walltime()
Tp_poly = MM.hecke_polynomial_in_T_variable(p,basis=B,verbose=False)
b = Tp_poly.padded_list()[1]
c = Tp_poly.padded_list()[0]
print "Time elapsed:", walltime() - before

print "The characteristic polynomial of U_%s is: %s"%(p,Tp_poly)
print "The discriminant of char poly of U_%s is: %s"%(p,b^2-4*c)
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



