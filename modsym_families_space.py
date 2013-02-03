import fund_domain
import modsym_space

class ModsymFamiliesSpace(Parent):

    def __init__(self,p,N,M,d_max,r,char=None):
        #initialize the manin relations
        #do other stuff??
        Parent.__init__(self,category=ModsymSpace)
        Element = ModsymFamiliesElement
        self.manin_relations = ManinRelations(N*p)#????????
        self.coefficient_module = FamiliesSpace(p,M,d,r,prec,char)
        
        
    def ms(self):
        """demotes to a regular modular symbol"""
        return modsym(self.level,self.data,self.manin)

    def num_moments(self):
        return self.data[0].num_moments()


    ## This function returns a number between 0 and p-2 which represents which     
    ## disc in weight the family is supported on -- the integer i corresponds to the
    ## disc with characters with finite order part omega^i
    def disc(self):
        return self.data[0].disc()


    
    ## This procedure tries to find a power series c(w) such that 
    ##      self | T_q = c(w) self
    ## Returns a triple consisting of a boolean and if true c(w) and the precision (if false, None and None)
#######################################################################################################################
##  This function produces a random family of OMSs.
##
##  p -- prime
##  N -- tame level
##  char -- character of conductor dividing N (or maybe Np?)
##  M -- number of moments
##  r -- integer between 0 and p-2 indicating which disc on weight space we are working
##  w -- variable which is just carried around because I had trouble with symbols being defined over different rings
#####################################################################################################################
    def random_element(TestTorsion=False):
        if self.char == None:
            self.char = DirichletGroup(1,QQ)[0]
        manin = self.manin_relations
        t2 = 0
        t3 = 0
        if (self.p == 2) and (len(manin.gens()) > 0):
            t2 = 1
        if (self.p == 3) and (len(manin.gens()) > 0):
            t3 = 1
        MM = M + t2 + t3 + 1 + floor(log(self.M)/log(self.p))

        v = []
        ## this loop runs through each generator (different from D_infty) and picks a random value for that generator
        for g in manin.gens():

            mus = coefficient_module.random_element()
            if not (g in manin.reps_with_two_tor) and not (g in manin.reps_with_three_tor):
                        v = v + [mus]
            elif (g in manin.reps):
                ## Case of two torsion (See [PS] section 4.1)
                v = v + [1/2*(mus - mus*g)]
            else:
                ## Case of three torsion (See [PS] section 4.1)    
                v = v + [2*mus - mus*g - 1/3*mus*gam**2]

        t = coefficient_module.zero()
        ## This loops adds up around the boundary of fundamental domain except the two verticle lines
        for g in manin.gens():
            if  not (g in manin.reps_with_two_tor) and not (g in manin.reps_with_three_tor):
                t = t + v[j-1]*g - v[j-1]
            else:
                t = t - v[j-1]

        ## t now should be sum Phi(D_i) | (gamma_i - 1) - sum Phi(D'_i) - sum Phi(D''_i)
        ## (Here I'm using the opposite sign convention of [PS1] regarding D'_i and D''_i)

        ## We now need to make some adjustment of Phi(D_i) to make t have total measure 0



        for g in manin.gens():
            gam=g
            if not (g in manin.reps_with_two_tor) and not (g in manin.reps_with_three_tor):
                break


        a = gam[0][0]
        c = gam[1][0]
            #WHERE IN THE WORLD IS THIS FUNCTION!!! 
        K = aut(p,deg,M,a,c,r,char,w)
        t0 = t.moment(0)

        err = coefficient_module.zero()

        if r != 0:
                ## The following code simply computes -t0/(K0-1)
                temp = coefficient_module.A(-t0/(K[0]))
                temp = temp.truncate(deg)
                err.moments[0] = temp
        else:
            ## The following code simply computes -t0/(K1)
                temp = coefficient_module.A(-t0/(K[1]))
                temp = temp.truncate(deg)
                temp=temp.truncate(deg)
                err.moments[1] = temp

        if (g in manin.reps_with_two_tor) or (g in manin.reps_with_three_tor):
            print "All generators are two and three torsion"
            v[j-1] = v[j-1] + err*gam - err
            t = t + err*gam - err
        else:
            v[j-1] =v[j-1] + err
            t = t + err*gam - err

        if TestTorsion:
            print "Total Measure", t.moment(0),t.moment(0).valuation(p)
            print "Is t zero?", t.is_zero()
        
            #print "is two/three tor?", manin.twotor.count(rj) > 0, manin.threetor.count(rj) >0
        #v[j-1] = v[j-1] + err
        #t = t + err.act_right(gam)-err
            
        mus = coefficient_module.solve_diff_eqn(t)

        v = [mus*(-1)] + v

        Phis = ModsymFamiliesElement(N*p,v,manin)
        Psis = Phis.change_precision(M)
        e = max(-Psis.valuation(),0)
        Phis = (p^e*Phis).change_precision(M)

        return Phis


