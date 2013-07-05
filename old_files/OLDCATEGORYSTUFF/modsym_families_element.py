#IMPLEMENT Hecke
class ModsymFamiliesElement:
	
	

	def is_Tq_eigensymbol(self,q,verbose=False):
		p = self.parent.prime()
		M = self.num_moments()
		R = self.parent.A
		T = PowerSeriesRing(QQ,'y')
		
		a = 0
		done = false
		for gen in self._map.gens():
			if _map[gen].moment(0)!=0:	
				
				Phiq = self.hecke(q)
				aq = R(Phiq._map[gen].moment(0)/self._map[gen].moment(0))
				if (self*aq - Phiq).normalize().is_zero():
					v=Sequence(f)
					v=[v[a]%(p^(M-self.valuation())) for a in range(len(v))]
					return (True, A(v), (M-self.valuation(), self.deg()))
				if verbose:
					print "The difference Phi | T_q - (potential eval) * Phi is:",(self.scale(aq) - Phiq).normalize()
				return (False, None, None)

		print "All of the total measures are zero!"
		return (False, None, None)


	

	def change_d(self,new_d):
		v=[self.data[r].change_d(new_d) for r in range(self.ngens())]
		return ModsymFamiliesElement(self.level,v,self.manin)
    
