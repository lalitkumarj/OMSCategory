#D = Distributions(0, 11, 10)
#MR = ManinRelations(11)
#data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
#f = ManinMap(D, MR, data)

for p in prime_range(3,100):
    srcs = ManinRelations(p)
    for k in range(0, 20, 1):
	print "\n", p, k
        ret = srcs._nice_gamma(p, k % (p-1))
        print "\n", ret[0]
        print "\n", ret[1]
        val = ret[1].matrix()[1,0].valuation(p)
        if val > 1:
            print "WHOA!", val

for p in prime_range(3,7):
    MR = ManinRelations(p)
    for k in range(0, p, 2):
        print MR._nice_gamma(p, k % (p-1))