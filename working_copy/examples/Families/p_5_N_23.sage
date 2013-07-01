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
print "In weight 2, there are 9 ordinary forms, 2 old and 7 new.  The oldforms are conjugate over Q(sqrt(5)) and thus congruent."
print "For the newforms: 1 is defined over Z_5 and is congruent to the 2 oldforms, 2 are conjugate over Q(sqrt(5)) (and thus are congruent), 4 are conjugate over a quartic extension K/Q in which 5 is inert."
print "Moreover, there are no other congruences.  Thus, the Hida algebra has 3 local pieces: T = T_2 + T_3 + T_4 labeled by rank."
print "We know that T_4 is Lambda \otimes O where O is a quartic unramified extension of Z_5."
print "\nIn what follows, we compute the characteristic polynomial of various Hecke operators projected to T_2 and T_3."
print "We are using precision: %s"%([M,var_prec])

B2 = []
B3 = []

while (len(B2) < 2) or (len(B3) < 3):
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
    for i in range(M+3):
        #print "    Iteration %s of Up"%(i+1)
        sys.stdout.flush()
        Phis = Phis.hecke(p)
    print "Time elapsed:", walltime() - before

    print "Positive valuation after projection.  Scaling away"
    val = Phis.valuation()
    if val > 0:
        for key in Phis._map._dict.keys():
            Phis._map._dict[key].ordp -= val

    print "  Killing off rank 4 piece (applying U_5 - 1)"
    sys.stdout.flush()
    before = walltime()
    for i in range(M+3):
        Phis = Phis.hecke(5) - Phis
    print "Time elapsed:", walltime() - before

    print "Positive valuation after projection.  Scaling away"
    val = Phis.valuation()
    if val > 0:
        for key in Phis._map._dict.keys():
            Phis._map._dict[key].ordp -= val

    if len(B2) < 2:
        print "Isolating T_2"
        Phis2 = Phis
        sys.stdout.flush()

        print "  Killing off rank 3 piece (applying T_2 - 2)"
        sys.stdout.flush()
        before = walltime()
        for i in range(M+3):
            Phis2 = Phis2.hecke(2) - 2 * Phis2
        print "Time elapsed:", walltime() - before

        Phis2.normalize()

        print "Positive valuation after projection.  Scaling away"
        val = Phis2.valuation()
        if val > 0:
            for key in Phis2._map._dict.keys():
                Phis2._map._dict[key].ordp -= val

        if MM.is_start_of_basis(B2 + [Phis2]):
            B2 = B2 + [Phis2]
            print "Adding this new symbol to our basis"
        else:
            print "FAILED. Symbol combined with previous basis does not generate free module"

    if len(B3) < 3:
        print "Isolating T_3"
        Phis3 = Phis
        sys.stdout.flush()

        print "  Killing off rank 2 piece (applying T_3 + 1)"
        sys.stdout.flush()
        before = walltime()
        for i in range(M+3):
            Phis3 = Phis3.hecke(3) + Phis3
        print "Time elapsed:", walltime() - before

        Phis3.normalize()

        print "Positive valuation after projection.  Scaling away"
        val = Phis3.valuation()
        if val > 0:
            for key in Phis3._map._dict.keys():
                Phis3._map._dict[key].ordp -= val

        if MM.is_start_of_basis(B3 + [Phis3]):
            B3 = B3 + [Phis3]
            print "Adding this new symbol to our basis"
        else:
            print "FAILED. Symbol combined with previous basis does not generate free module"

print "\nComputing U_%s on rank 2 piece:"%(p)
sys.stdout.flush()
before = walltime()
Tp_poly = MM.hecke_polynomial_in_T_variable(p,basis=B2,verbose=False)
b = Tp_poly.padded_list()[1]
c = Tp_poly.padded_list()[0]
print "Time elapsed:", walltime() - before

print "The characteristic polynomial of U_%s is: %s"%(p,Tp_poly)
print "The discriminant of char poly of U_%s is: %s"%(p,b^2-4*c)
print "The Iwasawa invariants (mu,lambda) of the discriminant are",(Iwasawa_invariants(b^2-4*c))


for q in primes(10):
    if q != p:
        print "\nComputing T_%s on rank 2 piece:"%(q)
        sys.stdout.flush()
        before = walltime()
        Tq_poly = MM.hecke_polynomial_in_T_variable(q,basis=B2,verbose=False)
        b = Tq_poly.padded_list()[1]
        c = Tq_poly.padded_list()[0]
        print "Time elapsed:", walltime() - before

        print "The characteristic polynomial of T_%s is: %s"%(q,Tq_poly)
        print "The discriminant of char poly of T_%s is: %s"%(q,b^2-4*c)
        print "The Iwasawa invariants (mu,lambda) of the discriminant are",Iwasawa_invariants(b^2-4*c)
        print "Time elapsed:", walltime() - before

print "\nComputing U_%s on rank 3 piece:"%(p)
sys.stdout.flush()
before = walltime()
Tp_poly = MM.hecke_polynomial_in_T_variable(p,basis=B3,verbose=False)
b = Tp_poly.padded_list()[2]
c = Tp_poly.padded_list()[1]
d = Tp_poly.padded_list()[0]
print "Time elapsed:", walltime() - before

print "The characteristic polynomial of U_%s is: %s"%(p,Tp_poly)
print "The discriminant of char poly of U_%s is: %s"%(p,b^2*c^2 - 4*c^3 - 4*b^3*d - 27*d^2 + 18*b*c*d)
print "The Iwasawa invariants (mu,lambda) of the discriminant are",(Iwasawa_invariants(b^2*c^2 - 4*c^3 - 4*b^3*d - 27*d^2 + 18*b*c*d))


for q in primes(10):
    if q != p:
        print "\nComputing T_%s on rank 3 piece:"%(q)
        sys.stdout.flush()
        before = walltime()
        Tq_poly = MM.hecke_polynomial_in_T_variable(q,basis=B3,verbose=False)
        b = Tq_poly.padded_list()[2]
        c = Tq_poly.padded_list()[1]
        d = Tq_poly.padded_list()[0]
        print "Time elapsed:", walltime() - before

        print "The characteristic polynomial of T_%s is: %s"%(q,Tq_poly)
        print "The discriminant of char poly of T_%s is: %s"%(q,b^2*c^2 - 4*c^3 - 4*b^3*d - 27*d^2 + 18*b*c*d)
        print "The Iwasawa invariants (mu,lambda) of the discriminant are",Iwasawa_invariants(b^2*c^2 - 4*c^3 - 4*b^3*d - 27*d^2 + 18*b*c*d)
        print "Time elapsed:", walltime() - before





