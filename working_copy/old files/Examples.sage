from sage.modular.pollack_stevens.coeffmod_OMS_families_space import FamiliesOfOverconvergentDistributions
DD = FamiliesOfOverconvergentDistributions(0, base=PowerSeriesRing(ZpCA(3, 10), 'w', default_prec=5))
print DD
d = DD([1]*10)
print "d =", d
print " "
S0 = DD.action().actor()
g = S0([1,1,3,1])
print "d*g =", d*g, "\n"
h = S0([1,0,0,1])
d2 = d*h
print "d2 =", d2
print "This should equal d. Does it?", d == d2
print "2*d =", 2*d

from sage.modular.pollack_stevens.coeffmod_OMS_space import OverconvergentDistributions
D = OverconvergentDistributions(0, base=ZpCA(3, 10))
print D
d = D([1]*10)
print "d =", d
print " "
S0 = D.action().actor()
g = S0([1,1,3,1])
print "d*g =", d*g, "\n"
h = S0([1,0,0,1])
d2 = d*h
print "d2 =", d2
print "This should equal d. Does it?", d == d2
print "2*d =", 2*d

from sage.modular.pollack_stevens.modsym_OMS_families_space import FamiliesOfOMS
from sage.modular.pollack_stevens.coeffmod_OMS_families_space import FamiliesOfOverconvergentDistributions
DD = FamiliesOfOverconvergentDistributions(0, base=PowerSeriesRing(ZpCA(3, 10), 'w', default_prec=5))
M = FamiliesOfOMS(11, 0, coefficients=DD)

from sage.modular.pollack_stevens.modsym_OMS_space import OverconvergentModularSymbols
from sage.modular.pollack_stevens.coeffmod_OMS_space import OverconvergentDistributions
p = 11
D = OverconvergentDistributions(0, base=ZpCA(p, 10))
M = OverconvergentModularSymbols(11, 0, coefficients=D)
MR = M.source()
Phi = M.random_element()
Psi = Phi.hecke(p)
for i in range(10):
    print i
    Psi = Psi.hecke(p)
Psi = Psi.hecke(3) - (1+3)*Psi
E = EllipticCurve('11a')
R = M.base_ring()
for q in prime_range(30):
    print q, R(E.ap(q)), Psi.Tq_eigenvalue(q), R(E.ap(q)) == Psi.Tq_eigenvalue(q)

set_verbose(2)
from sage.modular.pollack_stevens.modsym_OMS_families_space import FamiliesOfOMS
from sage.modular.pollack_stevens.coeffmod_OMS_families_space import FamiliesOfOverconvergentDistributions
p = 11
DD = FamiliesOfOverconvergentDistributions(0, base_coeffs=ZpCA(p, 6), prec_cap=[6,4])
MM = FamiliesOfOMS(11, 0, coefficients=DD)
MMR = MM.source()
Phis = MM.random_element()
Psis = Phis.hecke(p)
for i in range(8):
    print "*************  Iteration %s of Up"%(i+1)
    Psis = Psis.hecke(p)
Psis = Psis.hecke(p) - Psis

E = EllipticCurve('11a')
R = M.base_ring()
for q in prime_range(30):
    print q, R(E.ap(q)), Psi.Tq_eigenvalue(q), R(E.ap(q)) == Psi.Tq_eigenvalue(q)