from sage.modular.pollack_stevens.coeffmod_OMS_families_space import FamiliesOfOverconvergentDistributions
DD = FamiliesOfOverconvergentDistributions(0, base=PowerSeriesRing(ZpCA(3, 10), 'w', default_prec=5));
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
