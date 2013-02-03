from sage.structure.sage_object import SageObject
from sage.modular.dirichlet import DirichletGroup
from sage.rings.all import ZZ, QQ
from sage.all import *
import modsym_coefficient_module_category
from families4 import WeightKAction_fam
#################################################################################################################
##  A family of distributions -- i.e. an element of D \hat{\otimes} A(W_r) -- is represented by a vector whose 
##  i-th entry is the i-th moment of the distribution (which is a power series in w).
#################################################################################################################

class FamiliesSpace(Parent):
	def __init__(self,p,M_max,d_max,r,prec,char=None):
	        """
		Initializes a familiy of distributions

		INPUT:
			- p -- prime number
			- deg -- one more than the maximal degree of the power series in w (i.e. working modulo w^deg)
			- r -- an integer from 0 to p-2 marking our specified disc in weight space (whose tame char is omega^r)
			- moments -- a vector of power series in w which are the moments of our distribution
			- char -- optional argument which is a Dirichlet character (denoting the nebentype character)
            - M the number of moments
            - d the maximum allowed w-adic precision

        OUTPUT:
        
        A family of distributions with data as specified by the inputs.

        	"""
		#
		Parent.__init__(self,category=MSCoefficientModule())
		Element = FamiliesElement
		
		self.A = PowerSeriesRing(Zp(p,prec),'w')
		self.p=p
		self.w = self.A.gen()
		self.M_max = M_max
		self.d_max=d_max
		self._disc=r     ##  this r is an integer from 0 to p-2 which represents which disc in weight space we are working on
		
		if char != None:
			self._char = char
		else:
			self._char = DirichletGroup(1,QQ)[0]
        act = WeightKAction_fam(self,self.char(),None,False)
        self.register_action(act)
			
	def char(self):
		return self._char

	def disc(self):
		return self._disc


	def zero(self):
		return FamiliesElement(vector([0 for i in range(self.M_max)]),self.M_max,self.d_max,self)
	
	def gen(self):
		""" Returns the variable of the moments of the distribution (i.e. w)"""
		return self.w

	def solve_diff_eqn(self,D):
		mus=self.zero()
		for j in range(1,D.num_moments()):		
			v=[Rational(0) for i in range(D.num_moments())]
			v[j]=Rational(1)
			#???????????????????CAN OF MEGA WORMSSSSSS
			mu=DistributionElementBase(v,self.M_max,self.d)
			nu=mu.solve_diff_eqn()
			mus=mus+nu.lift_to_dist_fam(self.deg,self.disc(),w).scale(D.moment(j))
		return mus


	def random_element(self,prec_cap):
		v = []
		pM = self.p ** self.M_max
		pjs = []
		comp_pjs = True
		cur_pow = -1	#for computing pjs
		for a in range(self.M_max):
			flist = []
			for j in range(self.d_max):
				if comp_pjs:
					test_pow = ceil(j * (self.p-2)/(self.p-1))
					if test_pow > cur_pow:
						cur_pow = test_pow
						pjs.append(self.p ** cur_pow)
					else:
						pjs.append(pjs[-1])
				flist.append(pjs[j] * ZZ(floor(random() * pM)))
			comp_pjs = False
			v.append((self.A)(flist))
		return FamiliesElement(vector(v),self.M_max,self.d_max,self)




#@cached_function
# def form_acting_matrix_on_dist_fam(F):
# 	"""first row is just F, then F shifted over 1, etc."""
# 	list=copy(F)
# 	v=[]
# 	for j in range(0,len(list)):
# 		v=v+[copy(list)]
# 		list.insert(0,0)
# 		list.pop()
# 	return Matrix(v).transpose()

#DONT DELETE THIS
def normalize(F,p,r,N):
	v=F.list()
	M=ceil((N-r)*(p-2)/(p-1))
	v=[v[a]%(p^M) for a in range(len(v))]
	S=F.parent()
	return S(v)


## produces a random distribution with values in S_w  (i.e. the coefficient of w^j has valuation at least j*c_p
