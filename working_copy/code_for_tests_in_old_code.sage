# CODE BLOCK 1:
# Generates data for random testing in function test_Dist_solve_vs_old_code().

from sage.matrix.matrix_integer_2x2 import Matrix_integer_2x2 as mi2x2
import sage.matrix.all as matrix
M2Z = matrix.MatrixSpace(ZZ,2)
Delta_mat = mi2x2(M2Z, [1, 1, 0, 1], copy = False, coerce = False)

results_dict = {}
musnus_dict = {}
#results_str = "{"
ps = [3,5,7,13]
M = 7
for p in ps:
    print "p =", p
    R = Qp(p)
    V = R ** M
    for k in range(p):
        print "    k =", k
        mu = random_dist(p, k, M)
        mu.moments[0] = 0
        print "mu =", mu
        nu = mu.solve_diff_eqn()
        Vmu = V(mu.moments)
        Vnu = V(nu.moments)
	results_dict[(p, k)] = (Vmu, Vnu)
        musnus_dict[(p, k)] = (mu, nu)
	#results_str += repr((p, k)) + ": (Vs[(p, k)](" + repr(Vmu) + "), VQs[(p, k)][-1](" + repr(Vnu) + ")),"
        print nu.act_right(Delta_mat) - nu == mu
    print " "
#results_str = results_str[:-1] + "}"
pks = results_dict.keys()
pks.sort()
results_str = ""
for pk in pks:
    results_str += "\n    random_dict[{0}] = [Vs[{0}]({1}), VQs[{0}][-1]({2})]".format(pk, results_dict[pk][0], results_dict[pk][1])

print results_str

# Checking results of code block 1
for pk in pks:
    mu, nu = musnus_dict[pk]
    print nu.act_right(Delta_mat) - nu == mu


#############################

# CODE BLOCK 2:
# 

from sage.matrix.matrix_integer_2x2 import Matrix_integer_2x2 as mi2x2
import sage.matrix.all as matrix
M2Z = matrix.MatrixSpace(ZZ,2)
Delta_mat = mi2x2(M2Z, [1, 1, 0, 1], copy = False, coerce = False)

results_dict = {}
musnus_dict = {}
ps = [3,5,7,13]
M = 7
w = PowerSeriesRing(QQ, 'w').gen()
for p in ps:
    print "p =", p
    R = PowerSeriesRing(Qp(p), 'w')
    V = R ** M
    for k in range(p):
        print "    k =", k
        mus = random_dist_fam(p, M, M, k, w)
        mus.moments[0] = 0*w
        print "mus =", mus
        nus = mus.solve_diff_eqn()
        Vmu = V(mus.moments)
        Vnu = V(nus.moments)
	results_dict[(p, k)] = (Vmu, Vnu)
        musnus_dict[(p, k)] = (mus, nus)
	#results_str += repr((p, k)) + ": (Vs[(p, k)](" + repr(Vmu) + "), VQs[(p, k)][-1](" + repr(Vnu) + ")),"
        print nus.act_right(Delta_mat) - nus == mus
    print " "
#results_str = results_str[:-1] + "}"
pks = results_dict.keys()
pks.sort()
results_str = ""
for pk in pks:
    results_str += "\n    w = VQs[{0}][-1].base_ring().gen()".format(pk)
    results_str += "\n    random_dict[{0}] = [Vs[{0}]({1}), VQs[{0}][-1]({2})]".format(pk, results_dict[pk][0], results_dict[pk][1])

print results_str