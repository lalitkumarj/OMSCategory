from sage.structure.sage_object import SageObject

#################################################################################################################
##  A family of distributions -- i.e. an element of D \hat{\otimes} A(W_r) -- is represented by a vector whose 
##  i-th entry is the i-th moment of the distribution (which is a power series in w).
#################################################################################################################

class FamiliesElement(SageObject):
    def __init__(self,moments,M,d,parent):
            """
        Initializes a familiy of distributions

        INPUT:
            - p -- prime number
            - deg -- one more than the maximal degree of the power series in w (i.e. working modulo w^deg)
            - r -- an integer from 0 to p-2 marking our specified disc in weight space (whose tame char is omega^r)
            - moments -- a vector of power series in w which are the moments of our distribution
            - char -- optional argument which is a Dirichlet character (denoting the nebentype character)
            - M = p-adic precision
            - d the w-adic precision

        OUTPUT:
        
        A family of distributions with data as specified by the inputs.

            """
        self.parent = parent
	self.moments=moments
        self.M = self.parent.M_max
        self.d = self.parent.d_max

            
    def __repr__(self):
        return repr(self.moments)

    def moment(self,n):
        return self.moments[n]

    def num_moments(self):
        return len(self.moments)

    def truncate(self):
        v=self.moments
        w=[v[j].truncate(self.d) for j in range(self.num_moments())]
        return FamiliesElement(vector(w),self.M,self.d,self.parent)

    def change_d(self,new_d):
        assert new_d<=self.d, "can only lower degree"
        v=[self.moments[a].truncate(new_d) for a in range(self.num_moments())]
        return FamiliesElement(vector(v),self.M,new_d,self.parent)
    
    def __add__(self,right):
#        assert self.num_moments()==right.num_moments(), "the accuracies are different"
#        assert self.deg==right.deg, "the degrees in w are different"
        return FamiliesElement(self.moments+right.moments,self.M,self.d,self.parent)

    def scale(self,left):
        return FamiliesElement(left*self.moments,self.M,self.d,self.parent)

    def __sub__(self,right):
        return self+right.scale(-1)

    def __cmp__(self,right):
        return cmp((self.p,self.deg,self.disc(),self.moments,self.char()),(right.p,right.deg,self.disc(),right.moments,right.char()))

    def zero(self):
        return dist_fam(vector([self.moment(0)*0 for i in range(0,len(self.moments))]),self.M,self.d,self.parent)
    
    def is_zero(self):
        """Return true if all of self's moments are zero."""
        return self.moments.is_zero()

    def specialize(self,k):
        """evaluates at ((1+p)^k-1)/p"""
#        assert k % (self.parent.p - 1) == self.parent._disc, "Wrong component of weight space"
#        w=self.parent.gen()
#        v=[]
#        for j in range(0,self.num_moments()):
#            v=v+[Rational(self.moment(j).substitute(w=((1+self.parent.p)^k-1)/self.parent.p))]
#        return DistributionElement(self.parent.p,k,vector(v),self.char())
        raise NotImplementedError

    def valuation(self):
        return min([val(a,self.parent.p) for a in self.moment()])

    def normalize(self):
        """
        foo
        """
        N=self.num_moments()
        v=[]
        p = self.parent.p
        for j in range(N):
            m = self.moment(j)
            l = m.list()
            M = ceil((N-j)*(p-2)/(p-1))
            l = [a%(p^M) for a in l]
            S = m.parent()
            v = v+[S(v)]
        return FamiliesElement(vector(v),self.M,self.d,self.parent)
    
    def change_precision(self,M):
        """only hangs onto M moments"""
        assert M<=self.num_moments(),"not enough moments"
        v=[self.moment(i) for i in range(M)]
        mu=FamiliesElement(vector(v),self.M,self.d,self.parent)
        return mu

    def act_by_ps_fam(self,F):
        gam=form_acting_matrix_on_dist_fam(F)
        v=(Matrix(self.moments)*gam)[0]
        w=[v[j].truncate(self.deg) for j in range(self.num_moments())]
        return dist_fam(self.p,self.deg,self.disc(),vector(w),self.char())

    def act_right_weight_zero(self,gam):
        a=gam[0,0]
        b=gam[0,1]
        c=gam[1,0]
        d=gam[1,1]
        G=form_acting_matrix_on_dist(self.p,self.num_moments(),0,a,b,c,d)
        return FamiliesElement((Matrix(self.moments)*G)[0],self.M,self.d,self.parent)

    def act_right(self,gam):
        w=self.gen()
        K=aut(self.p,self.deg,self.num_moments(),gam[0,0],gam[1,0],self.disc(),self.char(),w)
        return self.act_by_ps_fam(K).act_right_weight_zero(gam)
        
    def solve_diff_eqn(self):
        w=self.gen()
        mus=self.zero()
        for j in range(1,self.num_moments()):        
            v=[Rational(0) for i in range(self.num_moments())]
            v[j]=Rational(1)
            mu=dist(self.p,0,v,self.char())
            nu=mu.solve_diff_eqn()
            mus=mus+nu.lift_to_dist_fam(self.deg,self.disc(),w).scale(self.moment(j))
        return mus

#@cached_function
def form_acting_matrix_on_dist_fam(F):
    """first row is just F, then F shifted over 1, etc."""
    list=copy(F)
    v=[]
    for j in range(0,len(list)):
        v=v+[copy(list)]
        list.insert(0,0)
        list.pop()
    return Matrix(v).transpose()

def val(F,p):
    v=F.list()
    if v==[]:
        return Infinity
    else:
        return min([valuation(v[j],p) for j in range(len(v))])

## produces a random distribution with values in S_w  (i.e. the coefficient of w^j has valuation at least j*c_p
def random_dist_fam(p,M,deg,r,w,char=None):
    v = []
    R = w.parent()
    pM = p ** M
    pjs = []
    comp_pjs = True
    cur_pow = -1    #for computing pjs
    for a in range(M):
        flist = []
        for j in range(deg):
            if comp_pjs:
                test_pow = ceil(j * (p-2)/(p-1))
                if test_pow > cur_pow:
                    cur_pow = test_pow
                    pjs.append(p ** cur_pow)
                else:
                    pjs.append(pjs[-1])
            flist.append(pjs[j] * ZZ(floor(random() * pM)))
        comp_pjs = False
        v.append(R(flist))
    return dist_fam(p,deg,r,vector(v),char).normalize()

def random_dist_fam_old(p,M,deg,r,w,char=None):
    v = []
    for a in range(M):
        f = 0*w
        for j in range(deg):
            f = f + p^(ceil(j*(p-2)/(p-1))) * ZZ(floor(random()*p^M))*w^j
        v = v + [f]
    return dist_fam(p,deg,r,vector(v),char).normalize()

