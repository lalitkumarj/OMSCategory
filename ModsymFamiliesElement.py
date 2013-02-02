class ModsymFamiliesElement(Element):
	def d(self):
		return self.data[0].d

	## This function returns a number between 0 and p-2 which represents which 	
	## disc in weight the family is supported on -- the integer i corresponds to the
	## disc with characters with finite order part omega^i
	def disc(self):
		return self.data[0].parent._disc

	def change_d(self,new_d):
		v=[self.data[r].change_d(new_d) for r in range(self.ngens())]
        return ModsymFamiliesElement(self.level,v,self.manin)
    
