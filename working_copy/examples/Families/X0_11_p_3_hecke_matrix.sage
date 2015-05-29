import sys

start_time = walltime()

#Initialization data
p = 3
N = 11
k = 2
r = (k-2) % (p-1)
M = 12  #Number of moments
var_prec = M    #Precision on the variable

print "This script computes the Up matrix for tame level {0}, p = {1} with {2} moments and w-adic precision {3}.".format(N, p, M, var_prec)
print "It also prints out the Hecke polynomial at p = {0}, as well as its discriminant.".format(p)

DD = FamiliesOfOverconvergentDistributions(r, base_coeffs=ZpCA(p, M), prec_cap=[M,var_prec])
MM = FamiliesOfOMS(N, r, coefficients=DD, sign=-1)  #With this sign, there is no Eisenstein component

print "\nFinding basis of ordinary subspace."
sys.stdout.flush()
before = walltime()
B = MM.basis_of_ordinary_subspace()
print "Basis found."
print "Time elapsed:", walltime() - before

print "\nFinding Hecke matrix at p."
sys.stdout.flush()
before = walltime()
H = MM.hecke_matrix(p, B)
print H
print "Time elapsed:", walltime() - before

print "\nFinding Hecke polynomial at p and its discriminant."
sys.stdout.flush()
before = walltime()
HP = H.charpoly()
print "\nChar poly:\n", HP
disc = HP[1]^2 - 4 * HP[0] * HP[2]
w = HP.parent().base_ring().gen()
RT.<T> = PowerSeriesRing(HP.parent().base_ring().base_ring())
disc_in_T = disc.subs({w : T/p})
print "\nDiscriminant:\n", disc_in_T
print "\nTime elapsed:", walltime() - before