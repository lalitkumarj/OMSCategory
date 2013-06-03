from sage.rings.integer_ring import ZZ
from sage.rings.padics.factory import Qp
from sage.functions.other import ceil
from sage.functions.log import log
from sage.modular.arithgroup.all import Gamma0
from sage.matrix.constructor import Matrix
from sage.modular.pollack_stevens.modsym_space import ModularSymbolSpace_generic
from sage.modular.pollack_stevens.coeffmod_OMS_families_space import FamiliesOfOverconvergentDistributions
from sage.modular.pollack_stevens.modsym_OMS_families_element import ModSym_OMS_Families_element
from sage.interfaces.gp import gp
from sage.modular.pollack_stevens.manin_map import ManinMap, M2Z
D = Distributions(0, 11, 10)
MR = ManinRelations(11)
data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
f = ManinMap(D, MR, data)

# for p in primes(100):
# 	G = Gamma0(p)
# 	srcs = ManinRelations(G.level())
# 	for k in range(0):
# 		_nice_gamma(p,k)
	
