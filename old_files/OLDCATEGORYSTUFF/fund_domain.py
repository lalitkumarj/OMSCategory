"""
Manin Relations

Code to create the Manin Relations class, which solves the "Manin
relations".  That is, a description of `Div^0(P^1(Q))` as a `\ZZ[\Gamma_0(N)]`-module
in terms of generators and relations is found.  The method used is
geometric, constructing a nice fundamental domain for `\Gamma_0(N)` and
reading the relevant Manin relations off of that picture.  The
algorithm follows the paper of Pollack and Stevens "Overconvergent
modular symbols and p-adic L-functions".
"""

######################################################################
##  Copyright (c) 2012, Rob Pollack and Jonathan Hanke
##      <rpollack@math.bu.edu>
##      <jonhanke@gmail.com>
##
##  Released under the GNU Public License, 2012.
##
######################################################################

from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
from sage.modular.modsym.all import P1List
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
from sage.modules.free_module_element import zero_vector
from copy import deepcopy
from sage.misc.cachefunc import cached_method
from sage.rings.arith import convergents,xgcd,gcd


M2ZSpace = MatrixSpace_ZZ_2x2()
def M2Z(x):
    """
    Creates an immutable 2x2 integer matrix

    INPUT:

    - ``x`` -- a list of four integers.

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.fund_domain import M2Z
        sage: A = M2Z([1,2,3,4])
        sage: hash(A)
        8
    """
    a = M2ZSpace(x)
    a.set_immutable()
    return a

Id = M2Z([1,0,0,1])
sig = M2Z([0,1,-1,0])
tau = M2Z([0,-1,1,-1])
minone_inf_path = M2Z([1,1,-1,0])

# We store these so that we don't have to constantly create them.
t00 = (0,0)
t10 = (1,0)
t01 = (0,1)
t11 = (1,1)

class PSModularSymbolsDomain(SageObject):
    def __init__(self, N, reps, indices, rels, equiv_ind):
        """
        INPUT:

        - `N` -- positive integer
        - ``reps`` -- TODO
        - ``indices`` -- TODO
        - ``rels`` -- TODO
        - ``equiv_ind`` -- TODO

        EXAMPLES::

        TODO: some good examples

        The level N must be an integer::

            sage: from sage.modular.pollack_stevens.fund_domain import PSModularSymbolsDomain
            sage: PSModularSymbolsDomain(1/2, None, None, None, None)
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: PSModularSymbolsDomain(Gamma0(11), None, None, None, None)
            Traceback (most recent call last):
            ...
            TypeError: unable to coerce <class 'sage.modular.arithgroup.congroup_gamma0.Gamma0_class_with_category'> to an integer
        """
        ## Store the level
        self._N = ZZ(N)

        ## Coset representatives of Gamma_0(N) coming from the geometric
        ## fundamental domain algorithm
        self._reps = reps

        ## This is a list of indices of the (geometric) coset representatives
        ## whose values (on the associated degree zero divisors) determine the
        ## modular symbol.
        self._indices = sorted(indices)

        self._gens = [reps[i] for i in self._indices]
        self._ngens = len(indices)
        self._rels = rels
        self._rel_dict = {}
        for j, L in enumerate(rels):
            self._rel_dict[reps[j]] = [(d, A, reps[i]) for (d, A, i) in L]
        ## A list of lists of triples (d, A, i), one for each coset
        ## representative of Gamma_0(N) (ordered to correspond to the
        ## representatives of self.reps) expressing the value of a
        ## modular symbol on the associated unimodular path as a sum of terms
        ##    d * (value on the i-th coset rep) | A
        ## where the index i must appear in self.gens_index, and the slash gives the
        ##  matrix action.

        self._equiv_ind = equiv_ind
        self._equiv_rep = {}
        for ky in equiv_ind:
            self._equiv_rep[ky] = reps[equiv_ind[ky]]

    def __len__(self):
        """
        Returns the number of coset representatives.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: len(A)
            12
        """
        return len(self._reps)

    def __getitem__(self, i):
        """
        Returns the `i`-th coset rep.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A[4]
            [-1 -2]
            [ 2  3]
        """
        return self._reps[i]

    def __iter__(self):
        """
        Returns an iterator over all coset representatives.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: for rep in A:
            ...       if rep[1,0] == 1:
            ...           print rep
            [ 0 -1]
            [ 1  3]
            [ 0 -1]
            [ 1  2]
            [ 0 -1]
            [ 1  1]
        """
        return iter(self._reps)

    def gens(self):
        """
        Returns the list of coset representatives chosen as generators.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.gens()
            [
            [1 0]  [ 0 -1]  [-1 -1]
            [0 1], [ 1  3], [ 3  2]
            ]
        """
        return self._gens

    def gen(self, n=0):
        """
        Returns the `n`-th generator.

        INPUT:

        - ``n`` -- integer (default: 0), which generator is desired

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(137)
            sage: A.gen(17)
            [-4 -1]
            [ 9  2]
        """
        return self._gens[n]

    def ngens(self):
        """
        Returns the number of generators.

        OUTPUT:

        - the number of coset representatives from which a modular
          symbol's value on any coset can be derived.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(1137)
            sage: A.ngens()
            255
        """
        return len(self._gens)

    def level(self):
        r"""
        Returns the level `N` of `\Gamma_0(N)` that we work with.

        OUTPUT:

        - The integer `N` of the group `\Gamma_0(N)` for which the
          Manin Relations are being computed.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.level()
            11
        """
        return self._N

    def indices(self, n=None):
        r"""
        Returns the indices of coset reps which were chosen as our
        generators.

        In particular, the divisors associated to these coset reps
        generate all divisors over `\ZZ[\Gamma_0(N)]`, and thus a modular
        symbol is uniquely determined by its values on these divisors.

        INPUT:

        - ``n`` -- integer (default: None)

        OUTPUT:

        - The list of indices in self.reps() of our generating
          set.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.indices()
            [0, 2, 3]
            sage: A.indices(2)
            3
            sage: A = ManinRelations(13)
            sage: A.indices()
            [0, 2, 3, 4, 5]
            sage: A = ManinRelations(101)
            sage: A.indices()
            [0, 2, 3, 4, 5, 6, 8, 9, 11, 13, 14, 16, 17, 19, 20, 23, 24, 26, 28]
        """
        if n is None:
            return self._indices
        else:
            return self._indices[n]

    def reps(self, n=None):
        r"""
        Returns the n-th coset rep associated with our fundamental
        domain or all coset reps if n is not specified.

        INPUT:

        - ``n`` -- integer (default: None)

        OUTPUT:

        - If n is given then the n-th coset representative is returned
          and otherwise all coset reps are returned.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.reps(0)
            [1 0]
            [0 1]
            sage: A.reps(1)
            [ 1  1]
            [-1  0]
            sage: A.reps(2)
            [ 0 -1]
            [ 1  3]
            sage: A.reps()
            [
            [1 0]  [ 1  1]  [ 0 -1]  [-1 -1]  [-1 -2]  [-2 -1]  [ 0 -1]  [ 1  0]
            [0 1], [-1  0], [ 1  3], [ 3  2], [ 2  3], [ 3  1], [ 1  2], [-2  1],
            <BLANKLINE>
            [ 0 -1]  [ 1  0]  [-1 -1]  [ 1 -1]
            [ 1  1], [-1  1], [ 2  1], [-1  2]
            ]
        """
        if n is None:
            return self._reps
        else:
            return self._reps[n]

    def relations(self, A=None, indices=False):
        r"""
        Expresses the divisor attached to the coset rep of A in terms
        of our chosen generators.

        INPUT:

        - ``A`` -- None, integer or a coset rep (default: None)

        - ``indices`` -- boolean (default: False), determines output
          type when ``A`` is None.

        OUTPUT:

        - A `\ZZ[\Gamma_0(N)]`-relation expressing the divisor
          attached to one (or all) coset rep(s) in terms of our
          generating set.  The type of the return value depends on
          ``A`` and ``indices``.

          - If ``A`` is a 2x2 matrix that is among the coset
            represetatives, returns a list of triples `(d, B, C)` such
            that the divisor attached to ``A`` equals the sum over
            these triples of:

              `d * B^{-1} * (divisor attached to C)`

            Here `C` will be one of the chosen generating coset reps.

          - If ``A`` is an integer, returns a list of triples `(d, B,
            i)` such that the divisor attached to the `i`-th coset rep
            equals the sum over these triples of:

              `d * B^{-1} * (divisor attached to i-th coset rep)`

            Here each index `i` must appear in ``self.indices()``.

          - If ``A`` is None and ``indices`` is False, returns a
            dictionary whose keys are the cosets reps and the values
            are the lists of triples `(d, B, C)` described above.

          - If ``A`` is None and ``indices`` is True, returns a list
            of triples `(d, B, i)` in the same order as the coset reps
            to which they correspond.

        .. NOTE::

            These relations allow us to recover the value of a modular
            symbol on any coset rep in terms of its values on our
            generating set.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.indices()
            [0, 2, 3]
            sage: A.relations(0)
            [(1, [1 0]
            [0 1], 0)]
            sage: A.relations(2)
            [(1, [1 0]
            [0 1], 2)]
            sage: A.relations(3)
            [(1, [1 0]
            [0 1], 3)]
            sage: A.relations(4)
            [(-1, [-3 -2]
            [11  7], 2)]
            sage: B=A.relations(4)[0][1]; B
            [-3 -2]
            [11  7]
            sage: B^(-1)*A.reps(2)
            [ 2 -1]
            [-3  2]
            sage: A.reps(4)
            [-1 -2]
            [ 2  3]
            sage: from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
            sage: M2Z = MatrixSpace_ZZ_2x2()
            sage: sig = M2Z([0,1,-1,0])
            sage: B^(-1)*A.reps(2) == A.reps(4)*sig
            True
        """
        if A is None:
            if indices:
                return self._rels
            else:
                return self._rel_dict
        if isinstance(A, (int, Integer, slice)):
            return self._rels[A]
        else:
            return self._rel_dict[A]


### Normalize elements of P^1(Z/N) for N arbitrary in ZZ (no overflows)
def p1_normalize_arbitrary(N, u, v,compute_s = False):
    r"""
    p1_normalize_arbitrary(N, u, v):
    
    Computes the canonical representative of
    `\mathbb{P}^1(\ZZ/N\ZZ)` equivalent to `(u,v)` along
    with a transforming scalar 's' (if compute_s is 1).
    
    INPUT:
    
    
    -  ``N`` - an integer (the modulus or level)
    
    -  ``u`` - an integer (the first coordinate of (u:v))
    
    -  ``v`` - an integer (the second coordinate of (u:v))
    
    -  ``compute_s`` - a boolean (int)
    
    
    OUTPUT: If gcd(u,v,N) = 1, then returns
   
    
    -  ``uu`` - an integer
    
    -  ``vv`` - an integer
    
    - ``ss`` - an integer such that `(ss*uu, ss*vv)` is equivalent to `(u,v)` mod `N`;

       if `\gcd(u,v,N) \not= 1`, returns 0, 0, 0.

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.fund_domain import p1_normalize_arbitrary
        sage: p1_normalize_arbitrary(90,7,77)
        (1, 11)
        sage: p1_normalize_arbitrary(90,7,77, compute_s=True)
        (1, 11, 7)
        sage: p1_normalize_arbitrary(90,7,78, True)
        (1, 24, 7)
        sage: (7*24-78*1) % 90
        0
        sage: (7*24) % 90
        78
    """
    if N == 1:
        if compute_s == True:
            return 0,0,1
        else:
            return 0,0

    u = u % N
    v = v % N
    if u<0: u += N
    if v<0: v += N
    if u == 0:
        uu = 0
        if gcd(v,N) == 1:
            vv = 1
        else:
            vv = 0
        ss = v
        if compute_s == True:
            return uu,vv,ss
        else:
            return uu,vv

    g,s,t = xgcd(u, N)
    s = s % N
    t = t % N
    if s<0: s += N
    if gcd(g, v) != 1:
        if compute_s == True:
            return 0, 0, 0
        else:
            return 0, 0

    # Now g = s*u + t*N, so s is a "pseudo-inverse" of u mod N
    # Adjust s modulo N/g so it is coprime to N.
    if g!=1:
        d = N/g
        while gcd(s,N) != 1:
            s = (s+d) % N

    # Multiply [u,v] by s; then [s*u,s*v] = [g,s*v] (mod N)
    u = g
    # v = (s*v) % N
    v = (s*v) % N

    min_v = v; min_t = 1
    if g!=1:
        Ng = N/g
        vNg = ZZ((v*Ng) % N)
        t = 1
        for k in range(2,g+1):
            v = (v + vNg) % N
            t = (t + Ng) % N
            if v<min_v and gcd(t,N)==1:
                min_v = v; min_t = t
    v = min_v
    if u<0: u = u+N
    if v<0: v = v+N
    uu = u
    vv = v
    if compute_s:
        ss = (Zmod(N)(s*min_t)**(-1)).lift()
        ss = ss % N
        return uu,vv,ss
    else:
        return uu,vv

######################################
##  Define the Manin Relation Class ##
######################################

class ManinRelations(PSModularSymbolsDomain):
    """
    This class gives a description of Div^0(P^1(QQ)) as a
    `\ZZ[\Gamma_0(N)]`-module.

    INPUT:

    - ``N`` -- a positive integer

    EXAMPLES::



    ``MR.reps_with_two_torsion`` is a list of coset reps whose
    associated unimodular path contains a point fixed by a
    `\Gamma_0(N)` element of order 2 (where the order is computed in
    `PSL_2(Z)`).

    ``MR.indices_with_two_torsion`` gives the corresponding indices in
    ``MR.reps()``::

        sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
        sage: MR = ManinRelations(11)
        sage: MR.indices_with_two_torsion
        []
        sage: MR = ManinRelations(13)
        sage: MR.indices_with_two_torsion
        [3, 4]
        sage: MR.reps_with_two_torsion
        [
        [-1 -1]  [-1 -2]
        [ 3  2], [ 2  3]
        ]
        sage: MR.reps_with_two_torsion[0] == MR.reps(3)
        True
        sage: MR = ManinRelations(17)
        sage: MR.reps_with_two_torsion
        [
        [-3 -2]  [-3 -1]
        [ 5  3], [ 4  1]
        ]
        sage: path = MR.reps_with_two_torsion[0]; path
        [-3 -2]
        [ 5  3]

    The corresponding matrices of order two are contained in
    ``MR.two_torsion``::

        sage: A = MR.two_torsion[path]; A
        [ 21  13]
        [-34 -21]
        sage: A^2
        [-1  0]
        [ 0 -1]

    You can see that multiplication by A just interchanges the
    rational cusps determined by the columns of the matrix ``path``::

        sage: A * path
        [ 2 -3]
        [-3  5]

        sage: sorted(ManinRelations(13).two_torsion.values())
        [
        [  5   2]  [  8   5]
        [-13  -5], [-13  -8]
        ]

    ``MR.reps_with_three_torsion`` is a list of coset reps whose
    associated unimodular path contains a point fixed by a
    `\Gamma_0(N)` element of order 3 in the ideal triangle directly
    below that path.  Here the order is again computed in `PSL_2(Z)`).

    ``MR.indices_with_three_torsion`` gives the corresponding indices
    in ``MR.reps()``::

        sage: MR = ManinRelations(11)
        sage: MR.indices_with_three_torsion
        []
        sage: MR = ManinRelations(13)
        sage: MR.indices_with_three_torsion
        [2, 5]
        sage: MR.reps_with_three_torsion
        [
        [ 0 -1]  [-2 -1]
        [ 1  3], [ 3  1]
        ]
        sage: MR.reps_with_three_torsion[0] == MR.reps(2)
        True
        sage: MR = ManinRelations(17)
        sage: MR.reps_with_three_torsion
        []
        sage: MR = ManinRelations(103)
        sage: MR.reps_with_three_torsion
        [
        [-4 -1]  [-1 -5]
        [ 9  2], [ 2  9]
        ]
        sage: path = MR.reps_with_three_torsion[0]; path
        [-4 -1]
        [ 9  2]

    The corresponding matrices of order three are contained in
    ``MR.three_torsion``::

        sage: A = MR.three_torsion[path]; A
        [-47 -21]
        [103  46]
        sage: A^3
        [1 0]
        [0 1]

    You can see that the columns of path, A*path and A^2*path give the
    same rational cusps::

        sage: A * path
        [ -1   5]
        [  2 -11]
        sage: A^2*path
        [  5  -4]
        [-11   9]

        sage: sorted(ManinRelations(13).three_torsion.values())
        [
        [-10  -7]  [-4 -1]
        [ 13   9], [13  3]
        ]
    """
    def __init__(self, N):
        r"""
        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: ManinRelations(11)
            Manin Relations of level 11
            sage: type(ManinRelations(30))
            <class 'sage.modular.pollack_stevens.fund_domain.ManinRelations'>
            sage: ManinRelations(1)
            Manin Relations of level 1

        Error checking::

            sage: ManinRelations(0)
            Traceback (most recent call last):
            ...
            ValueError: N must be a positive integer
            sage: ManinRelations(-5)
            Traceback (most recent call last):
            ...
            ValueError: N must be a positive integer        

        Implementation limits::
        
            sage: ManinRelations(2^20)
            Traceback (most recent call last):
            ...
            OverflowError: Modulus is too large (must be < 46340)
        """
        N = ZZ(N)
        if N <= 0:
            raise ValueError, "N must be a positive integer"
        
        self._N = N

        ## Creates and stores the Sage representation of P^1(Z/NZ)
        P = P1List(N)
        self._P = P

        ## Creates a fundamental domain for Gamma_0(N) whose boundary is a union
        ## of unimodular paths (except in the case of 3-torsion).
        ## We will call the intersection of this domain with the real axis the
        ## collection of cusps (even if some are Gamma_0(N) equivalent to one another).
        cusps = self.form_list_of_cusps()

        ## Takes the boundary of this fundamental domain and finds SL_2(Z) matrices whose
        ## associated unimodular path gives this boundary.  These matrices form the
        ## beginning of our collection of coset reps for Gamma_0(N) / SL_2(Z).
        coset_reps = self.fd_boundary(cusps)

        ## Takes the bottom row of each of our current coset reps,
        ## thinking of them as distinct elements of P^1(Z/NZ)
        p1s = [(coset_reps[j])[1] for j in range(len(coset_reps))]

        ## Initializes relevant Manin data
        gens_index = []
        twotor_index = []
        twotorrels = []
        threetor_index = []
        threetorrels = []
        rels = [0] * len(coset_reps)

        ## the list rels (above) will give Z[Gamma_0(N)] relations between
        ## the associated divisor of each coset representatives in terms
        ## of our chosen set of generators.
        ## entries of rel will be lists of elements of the form (c,A,r)
        ## with c a constant, A a Gamma_0(N) matrix, and r the index of a
        ## generator.  The meaning is that the divisor associated to the
        ## j-th coset rep will equal the sum of:
        ##
        ##   c * A^(-1) * (divisor associated to r-th coset rep)
        ##
        ## as one varies over all (c,A,r) in rel[j].
        ## (Here r must be in self.generator_indices().)
        ##
        ## This will be used for modular symbols as then the value of a
        ## modular symbol phi on the (associated divisor) of the j-th
        ## element of coset_reps will be the sum of c * phi (r-th genetator) | A
        ## as one varies over the tuples in rel[j]

        boundary_checked = [False] * len(coset_reps)

        ## The list boundary_checked keeps track of which boundary pieces of the
        ## fundamental domain have been already used as we are picking
        ## our generators

        ## The following loop will choose our generators by picking one edge
        ## out of each pair of edges that are glued to each other and picking
        ## each edge glued to itself (arising from two-torsion)
        ## ------------------------------------------------------------------
        for r in range(len(coset_reps)):
            if boundary_checked[r] == False:

                ## We now check if this boundary edge is glued to itself by
                ## Gamma_0(N)

                if P.normalize(p1s[r][0],p1s[r][1]) == P.normalize(-p1s[r][1],p1s[r][0]):
                    ## This edge is glued to itself and so coset_reps[r]
                    ## needs to be added to our generator list.

                    ## this relation expresses the fact that
                    ## coset_reps[r] is one of our basic generators
                    rels[r] = [(1,Id,r)]

                    ## the index r is adding to our list
                    ## of indexes of generators
                    gens_index.append(r)

                    ## the index r is adding to our list of indexes of
                    ## generators which satisfy a 2-torsion relation
                    twotor_index.append(r)

                    gam = coset_reps[r] * sig * coset_reps[r]._invert_unit()
                    ## gam is 2-torsion matrix and in Gamma_0(N).
                    ## if D is the divisor associated to coset_reps[r]
                    ## then gam * D = - D and so (1+gam)D=0.

                    ## This gives a restriction to the possible values of
                    ## modular symbols on D

                    ## The 2-torsion matrix gam is recorded in our list of
                    ## 2-torsion relations.
                    twotorrels.append(gam)

                    ## We have now finished with this edge.
                    boundary_checked[r] = True

                else:
                    c = coset_reps[r][t10]
                    d = coset_reps[r][t11]

                    ## In the following case the ideal triangle below
                    ## the unimodular path described by coset_reps[r]
                    ## contains a point fixed by a 3-torsion element.
                    if (c**2+d**2+c*d)%N == 0:

                        ## the index r is adding to our list of indexes
                        ## of generators
                        gens_index.append(r)

                        ## this relation expresses the fact that coset_reps[r]
                        ## is one of our basic generators
                        rels[r] = [(1,Id,r)]

                        ## the index r is adding to our list of indexes of
                        ##generators which satisfy a 3-torsion relation
                        threetor_index.append(r)

                        gam = coset_reps[r] * tau * coset_reps[r]._invert_unit()
                        ## gam is 3-torsion matrix and in Gamma_0(N).
                        ## if D is the divisor associated to coset_reps[r]
                        ## then (1+gam+gam^2)D=0.
                        ## This gives a restriction to the possible values of
                        ## modular symbols on D

                        ## The 3-torsion matrix gam is recorded in our list of
                        ## 3-torsion relations.
                        threetorrels.append(gam)

                        ## The reverse of the unimodular path associated to
                        ## coset_reps[r] is not Gamma_0(N) equivalent to it, so
                        ## we need to include it in our list of coset
                        ## representatives and record the relevant relations.

                        a = coset_reps[r][t00]
                        b = coset_reps[r][t01]
                        A = M2Z([-b,a,-d,c])
                        coset_reps.append(A)
                        ## A (representing the reversed edge) is included in
                        ## our list of coset reps

                        rels.append([(-1,Id,r)])
                        ## This relation means that phi on the reversed edge
                        ## equals -phi on original edge

                        boundary_checked[r] = True
                        ## We have now finished with this edge.

                    else:

                        ## This is the generic case where neither 2 or
                        ## 3-torsion intervenes.
                        ## The below loop searches through the remaining edges
                        ## and finds which one is equivalent to the reverse of
                        ## coset_reps[r]
                        ## ---------------------------------------------------
                        for s in range(r+1, len(coset_reps)):
                            if boundary_checked[s]:
                                continue
                            if P.normalize(p1s[s][0],p1s[s][1]) == P.normalize(-p1s[r][1],p1s[r][0]):
                                ## the reverse of coset_reps[r] is
                                ## Gamma_0(N)-equivalent to coset_reps[s]
                                ## coset_reps[r] will now be made a generator
                                ## and we need to express phi(coset_reps[s])
                                ## in terms of phi(coset_reps[r])

                                gens_index.append(r)
                                ## the index r is adding to our list of
                                ## indexes of generators

                                rels[r] = [(1,Id,r)]
                                ## this relation expresses the fact that
                                ## coset_reps[r] is one of our basic generators

                                A = coset_reps[s] * sig
                                ## A corresponds to reversing the orientation
                                ## of the edge corr. to coset_reps[r]

                                gam = coset_reps[r] * A._invert_unit()
                                ## gam is in Gamma_0(N) (by assumption of
                                ## ending up here in this if statement)

                                rels[s] = [(-1,gam,r)]
                                ## this relation means that phi evaluated on
                                ## coset_reps[s] equals -phi(coset_reps[r])|gam
                                ## To see this, let D_r be the divisor
                                ## associated to coset_reps[r] and D_s to
                                ## coset_reps[s]. Then gam D_s = -D_r and so
                                ## phi(gam D_s) = - phi(D_r) and thus
                                ## phi(D_s) = -phi(D_r)|gam
                                ## since gam is in Gamma_0(N)

                                boundary_checked[r] = True
                                boundary_checked[s] = True
                                break

        ## We now need to complete our list of coset representatives by
        ## finding all unimodular paths in the interior of the fundamental
        ## domain, as well as express these paths in terms of our chosen set
        ## of generators.
        ## -------------------------------------------------------------------

        for r in range(len(cusps)-2):

        ## r is the index of the cusp on the left of the path.  We only run
        ## thru to the number of cusps - 2 since you can't start an interior
        ## path on either of the last two cusps

            for s in range(r+2,len(cusps)):
            ## s is in the index of the cusp on the the right of the path
                cusp1 = cusps[r]
                cusp2 = cusps[s]
                if self.is_unimodular_path(cusp1,cusp2):
                    A,B = self.unimod_to_matrices(cusp1,cusp2)
                    ## A and B are the matrices whose associated paths
                    ## connect cusp1 to cusp2 and cusp2 to cusp1 (respectively)
                    coset_reps.extend([A,B])
                    ## A and B are added to our coset reps
                    vA = []
                    vB = []

                    ## This loop now encodes the relation between the
                    ## unimodular path A and our generators.  This is done
                    ## simply by accounting for all of the edges that lie
                    ## below the path attached to A (as they form a triangle)
                    ## Similarly, this is also done for B.

                    for rel in rels[r+2:s+2]:
                    ## Running between the cusps between cusp1 and cusp2
                        ## Add edge relation
                        vA.append(rel[0])
                        ## Add negative of edge relation
                        vB.append((-rel[0][0], rel[0][1], rel[0][2]))
                    ## Add relations for A and B to relations list
                    rels.extend([vA,vB])

        ## Make the translation table between the Sage and Geometric
        ## descriptions of P^1
        equiv_ind = {}
        for i, rep in enumerate(coset_reps):
            ky = P.normalize(rep[t10],rep[t11])
            equiv_ind[ky] = i

        PSModularSymbolsDomain.__init__(self, N, coset_reps, gens_index, rels, equiv_ind)

        ## A list of indices of the (geometric) coset representatives whose
        ## paths are identified by some 2-torsion element (which switches the
        ## path orientation)
        self.indices_with_two_torsion = twotor_index
        self.reps_with_two_torsion = [coset_reps[i] for i in twotor_index]

        ## A dictionary of (2-torsion in PSL_2(Z)) matrices in Gamma_0(N) that give
        ## the orientation identification in the paths listed in twotor_index above!
        self.two_torsion = {}
        for j, tor_elt in zip(twotor_index, twotorrels):
            self.two_torsion[coset_reps[j]] = tor_elt

        ## A list of indices of the (geometric) coset representatives that
        ## form one side of an ideal triangle with an interior fixed point of
        ## a 3-torsion element of Gamma_0(N)
        self.indices_with_three_torsion = threetor_index
        self.reps_with_three_torsion = [coset_reps[i] for i in threetor_index]

        ## A dictionary of (3-torsion in PSL_2(Z)) matrices in Gamma_0(N) that give
        ## the interior fixed point described in threetor_index above!
        self.three_torsion = {}
        for j, tor_elt in zip(threetor_index, threetorrels):
            self.three_torsion[coset_reps[j]] = tor_elt

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: ManinRelations(11)._repr_()
            'Manin Relations of level 11'
        """
        return "Manin Relations of level %s"%self._N

    def equivalent_index(self, A):
        r"""
        Returns the index of the rep equivalent to A.

        Here by equivalent we mean the unique coset rep whose bottom
        row is equivalent to the bottom row of A in `P^1(\ZZ/N\ZZ)`.

        INPUT:

        - ``A`` -- an element of `SL_2(\ZZ)`

        OUTPUT:

        - The unique integer j satisfying that the bottom row of
          self.reps(j) is equivalent to the bottom row of A.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: MR = ManinRelations(11)
            sage: A = matrix(ZZ,2,2,[3,5,16,27])
            sage: j = MR.equivalent_index(A); j
            8
            sage: MR.reps(8)
            [ 0 -1]
            [ 1  1]
            sage: MR.equivalent_rep(A)
            [ 0 -1]
            [ 1  1]
            sage: MR.P1().normalize(16,27)
            (1, 1)
        """
        try:
            ky = self._P.normalize(A[t10],A[t11])
        except OverflowError:
            ky = p1_normalize_arbitrary(self._P.N(),A[t10],A[t11])

        return self._equiv_ind[ky]

    def equivalent_rep(self, A):
        """
        Returns a coset representative that is equivalent to A modulo `\Gamma_0(N)`.

        INPUT:

        - ``A`` -- a matrix in `SL_2(\ZZ)`

        OUTPUT:

        - a matrix in `SL_2(\ZZ)` congruent to ``A`` modulo `\Gamma_0(N)`.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
            sage: M2Z = MatrixSpace_ZZ_2x2()
            sage: A = M2Z([5,3,38,23])
            sage: ManinRelations(60).equivalent_rep(A)
            [-7 -3]
            [26 11]
        """
        try:
            ky = self._P.normalize(A[t10],A[t11])
        except OverflowError:
            ky = p1_normalize_arbitrary(self._P.N(),A[t10],A[t11])
        return self._equiv_rep[ky]

    def P1(self):
        r"""
        Returns the Sage representation of `P^1(\ZZ/N\ZZZ)`.

        OUTPUT:

        - `P^1(Z/NZ)` where N is the level of the relations.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.P1()
            The projective line over the integers modulo 11
        """
        return self._P

    def form_list_of_cusps(self):
        r"""
        Returns the intersection of a fundamental domain for
        `\Gamma_0(N)` with the real axis.

        The construction of this fundamental domain follows the
        arguments of [PS] Section 2.  The boundary of this fundamental
        domain consists entirely of unimodular paths when
        `\Gamma_0(N)` has no elements of order 3.  (See [PS] Section
        2.5 for the case when there are elements of order 3.)

        OUTPUT:

        - A sorted list of rational numbers marking the intersection
          of a fundamental domain for `\Gamma_0(N)` with the real
          axis.


        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.form_list_of_cusps()
            [-1, -2/3, -1/2, -1/3, 0]
            sage: A = ManinRelations(13)
            sage: A.form_list_of_cusps()
            [-1, -2/3, -1/2, -1/3, 0]
            sage: A = ManinRelations(101)
            sage: A.form_list_of_cusps()
            [-1, -6/7, -5/6, -4/5, -7/9, -3/4, -11/15, -8/11, -5/7, -7/10, -9/13, -2/3, -5/8, -13/21, -8/13, -3/5, -7/12, -11/19, -4/7, -1/2, -4/9, -3/7, -5/12, -7/17, -2/5, -3/8, -4/11, -1/3, -2/7, -3/11, -1/4, -2/9, -1/5, -1/6, 0]
        """
        ## Get the level
        N = self.level()

        ## Checks that the level N is > 1
        # TODO: I'm commenting this out; I see no reason not to allow level 1, except
        # possibly the bug here that I fixed: http://trac.sagemath.org/sage_trac/ticket/12772
        #if not (N > 1):
        #    raise TypeError, "Error in form_list_of_cusps: level should be > 1"

        ## Some convenient shortcuts
        P = self.P1()
        sP = len(P.list())   ## Size of P^1(Z/NZ)

        ## Initialize some lists

        C = [QQ(-1),"?",QQ(0)]

        ## Initialize the list of cusps at the bottom of the fund. domain.
        ## The ? denotes that it has not yet been checked if more cusps need
        ## to be added between the surrounding cusps.

        full_domain = False     ## Says that we're not done yet!

        v = [False for r in range(sP)]
        ## This initializes a list indexed by P^1(Z/NZ) which keeps track of
        ## which right coset representatives we've found for Gamma_0(N)/SL_2(Z)
        ## thru the construction of a fundamental domain

        ## Includeds the coset repns formed by the original ideal triangle
        ## (with corners at -1, 0, infty)

        v[P.index(0,1)] = True
        v[P.index(1,-1)] = True
        v[P.index(-1,0)] = True


        ## Main Loop -- Ideal Triangle Flipping
        ## ====================================
        while (not full_domain):
            full_domain = True

            ## This loop runs through the current set of cusps
            ## and checks to see if more cusps should be added
            ## -----------------------------------------------
            for s in range(1, len(C), 2):  ## range over odd indices in the
                                           ## final list C
                if C[s] == "?":

                    ## Single out our two cusps (path from cusp2 to cusp1)
                    cusp1 = C[s-1]
                    cusp2 = C[s+1]

                    ## Makes the unimodular transform for the path from cusp2
                    ## to cusp1

                    b1 = cusp1.denominator()
                    b2 = cusp2.denominator()

                    ## This is the point where it is determined whether
                    ## or not the adjacent triangle should be added
                    ## ------------------------------------------------
                    pos = P.index(b2,b1)   ## The Sage index of the bottom
                                                 ## row of our unimodular
                                           ## transformation gam

                    ## Check if we need to flip (since this P1 element has not
                    ## yet been accounted for!)
                    if v[pos] == False:
                        v[pos] = True      ## Say this P1 element now occurs
                        v[P.index(b1,-(b1+b2))] = True ## Say that the other
                                                       ## two ideal triangle
                                                       ## edges also occur!
                        v[P.index(-(b1+b2),b2)] = True

                        ## Check to see if this triangle contains a fixed
                        ## point by an element of Gamma_0(N).  If such an
                        ## element is present, the fundamental domain can be
                        ## extended no further.

                        if (b1**2 + b2**2 + b1*b2)%N != 0:

                        ## this congruence is exactly equivalent to
                        ## gam * [0 -1; 1 -1] * gam^(-1) is in Gamma_0(N)
                        ## where gam is the matrix corresponding to the
                        ## unimodular path connecting cusp1 to cusp2

                            C[s] = "i"  ## The '?' is changed to an 'i'
                                            ## indicating that a new cusp needs to
                                        ##  be inserted here
                            full_domain = False
                        else:
                            C[s] = "x"  ## The '?' is changed to an 'x' and no
                                               ##  more checking below is needed! =)
                    else:
                        C[s] = "x"  ## The '?' is changed to an 'x' and no more
                                           ## checking below is needed! =)


            ## Now insert the missing cusps (where there is an 'i' in the
            ## final list)
            ## This will keep the fundamental domain as flat as possible!
            ## ---------------------------------------------------------------

            s=1
            while s < len(C):    ## range over odd indices in the final list C
                if C[s] == "i":
                    C[s]="?"

                    ## Single out our two cusps (path from cusp2 to cusp1)
                    cusp1 = C[s-1]
                    cusp2 = C[s+1]

                    ## Makes the unimodular transform for the path from cusp2
                    ## to cusp1
                    a1 = cusp1.numerator()
                    b1 = cusp1.denominator()
                    a2 = cusp2.numerator()
                    b2 = cusp2.denominator()

                    ## Inserts the Farey center of these two cusps!
                    a = a1 + a2
                    b = b1 + b2
                    C.insert(s+1, a/b)
                    C.insert(s+2, "?")
                    s = s+2
                s = s+2

        ## Remove the (now superfluous) extra string characters that appear
        ## in the odd list entries
        C = [QQ(C[s]) for s in range(0,len(C),2)]
        return C

    def is_unimodular_path(self, r1, r2):
        r"""
        Determines whether two (non-infinite) cusps are connected by a
        unimodular path.

        INPUT:

        - ``r1, r2`` -- rational numbers

        OUTPUT:

        - A boolean expressing whether or not a unimodular path
          connects r1 to r2.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.is_unimodular_path(0,1/3)
            True
            sage: A.is_unimodular_path(1/3,0)
            True
            sage: A.is_unimodular_path(0,2/3)
            False
            sage: A.is_unimodular_path(2/3,0)
            False
        """
        a = r1.numerator()
        b = r2.numerator()
        c = r1.denominator()
        d = r2.denominator()
        return (a*d - b*c)**2 == 1


    def unimod_to_matrices(self, r1, r2):
        r"""
        Returns the two matrices whose associated unimodular paths
        connect `r1 -> r2` and `r2 -> r1`, respectively.

        INPUT:

        - ``r1, r2`` -- rational numbers (that are assumed to be
          related by a unimodular path)

        OUTPUT:

        - a pair of `2 x 2` matrices of determinant 1

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.unimod_to_matrices(0,1/3)
            (
            [ 0  1]  [1 0]
            [-1  3], [3 1]
            )

        """
        a = r1.numerator()
        b = r2.numerator()
        c = r1.denominator()
        d = r2.denominator()
        if (a*d-b*c)==1:
            ans = M2Z([a,b,c,d]), M2Z([-b,a,-d,c])
        else:
            ans = M2Z([-a,b,-c,d]), M2Z([b,a,d,c])
        ans[0].set_immutable()
        ans[1].set_immutable()
        return ans

    def fd_boundary(self,C):
        r"""
        Finds matrices whose associated unimodular paths give the
        boundary of a fundamental domain.

        Here the fundamental domain is for `\Gamma_0(N)`.  (In the
        case when `\Gamma_0(N)` has elements of order three the shape
        cut out by these unimodular matrices is a little smaller than
        a fundamental domain.  See `\S2.5` of Pollack-Stevens.)

        INPUT:

        - a list of rational numbers coming from
          self.form_list_of_cusps()

        OUTPUT:

        - a list of `2 x 2` integer matrices of determinant 1 whose
          associated unimodular paths give the boundary of a
          fundamental domain for `Gamma_0(N)` (or nearly so in the
          case of 3-torsion).

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: C = A.form_list_of_cusps(); C
            [-1, -2/3, -1/2, -1/3, 0]
            sage: A.fd_boundary(C)
            [
            [1 0]  [ 1  1]  [ 0 -1]  [-1 -1]  [-1 -2]  [-2 -1]
            [0 1], [-1  0], [ 1  3], [ 3  2], [ 2  3], [ 3  1]
            ]
            sage: A = ManinRelations(13)
            sage: C = A.form_list_of_cusps(); C
            [-1, -2/3, -1/2, -1/3, 0]
            sage: A.fd_boundary(C)
            [
            [1 0]  [ 1  1]  [ 0 -1]  [-1 -1]  [-1 -2]  [-2 -1]
            [0 1], [-1  0], [ 1  3], [ 3  2], [ 2  3], [ 3  1]
            ]
            sage: A = ManinRelations(101)
            sage: C = A.form_list_of_cusps(); C
            [-1, -6/7, -5/6, -4/5, -7/9, -3/4, -11/15, -8/11, -5/7, -7/10, -9/13, -2/3, -5/8, -13/21, -8/13, -3/5, -7/12, -11/19, -4/7, -1/2, -4/9, -3/7, -5/12, -7/17, -2/5, -3/8, -4/11, -1/3, -2/7, -3/11, -1/4, -2/9, -1/5, -1/6, 0]
            sage: A.fd_boundary(C)
            [
            [1 0]  [ 1  1]  [ 0 -1]  [-1 -1]  [-1 -2]  [-2 -1]  [-1 -3]  [-3 -2]
            [0 1], [-1  0], [ 1  6], [ 6  5], [ 5  9], [ 9  4], [ 4 11], [11  7],
            <BLANKLINE>
            [-2 -1]  [-1 -4]  [-4 -3]  [-3 -2]  [-2 -7]  [-7 -5]  [-5 -3]  [-3 -4]
            [ 7  3], [ 3 11], [11  8], [ 8  5], [ 5 17], [17 12], [12  7], [ 7  9],
            <BLANKLINE>
            [-4 -1]  [-1 -4]  [ -4 -11]  [-11  -7]  [-7 -3]  [-3 -8]  [ -8 -13]
            [ 9  2], [ 2  7], [  7  19], [ 19  12], [12  5], [ 5 13], [ 13  21],
            <BLANKLINE>
            [-13  -5]  [-5 -2]  [-2 -9]  [-9 -7]  [-7 -5]  [-5 -8]  [ -8 -11]
            [ 21   8], [ 8  3], [ 3 13], [13 10], [10  7], [ 7 11], [ 11  15],
            <BLANKLINE>
            [-11  -3]  [-3 -7]  [-7 -4]  [-4 -5]  [-5 -6]  [-6 -1]
            [ 15   4], [ 4  9], [ 9  5], [ 5  6], [ 6  7], [ 7  1]
            ]
        """

        C.reverse() ## Reverse here to get clockwise orientation of boundary

        ## These matrices correspond to the paths from infty to 0 and -1 to infty
        mats = [Id, minone_inf_path]

        ## Now find SL_2(Z) matrices whose associated unimodular paths connect
        ## the cusps listed in C.
        ## --------------------------------------------------------
        for j in range(len(C)-1):
            a = C[j].numerator()
            b = C[j+1].numerator()
            c = C[j].denominator()
            d = C[j+1].denominator()
            new_mat = M2Z([a,b,c,d])
            mats.append(new_mat)

        return mats

    @cached_method
    def prep_hecke_on_gen(self, ell, gen):
        """
        This function does some precomputations needed to compute T_ell.

        In particular, if phi is a modular symbol and D_m is the divisor associated to the generator ``gen``,
        to compute (phi|T_ell)(D_m) one needs to compute phi(gam_a D_m)|gam_a where
        gam_a runs thru the ell+1 matrices defining T_ell.  One then takes gam_a D_m and writes it
        as a sum of unimodular divisors.  For each such unimodular divisor, say [M] where M is a
        SL_2 matrix, we then write M=gam*h where gam is in Gamma_0(N) and h is one of our
        chosen coset representatives.  Then phi([M]) = phi([h]) | gam^(-1).  Thus, one has

            (phi | gam_a)(D_m) = sum_h sum_j phi([h]) | gam_{hj}^(-1) * gam_a

        as h runs over all coset representatives and j simply runs over however many
        times M_h appears in the above computation.

        Finally, the output of this function is a dictionary D whose keys are the coset representatives
        in ``self.reps()`` where each value is a list of matrices, and the entries of D
        satisfy:

            D[h][j] = gam_{hj} * gam_a

        INPUT:

        - ``ell`` -- a prime
        - ``gen`` -- a generator

        OUTPUT:

        A list of lists (see above).

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
        sage: phi = ps_modsym_from_elliptic_curve(E)
        sage: phi.values()
        [-1/5, 3/2, -1/2]
        sage: M = phi.parent().source()
        sage: M.prep_hecke_on_gen(2, M.gens()[0])
        {[ 1  0]
        [-1  1]: [], [1 0]
        [0 1]: [[1 0]
        [0 2], [1 1]
        [0 2], [2 0]
        [0 1]], [ 1 -1]
        [-1  2]: [[ 1 -1]
        [ 0  2]], [ 1  0]
        [-2  1]: [], [ 0 -1]
        [ 1  1]: [], [-1 -2]
        [ 2  3]: [], [ 0 -1]
        [ 1  3]: [], [-1 -1]
        [ 2  1]: [], [ 0 -1]
        [ 1  2]: [], [-2 -1]
        [ 3  1]: [], [ 1  1]
        [-1  0]: [], [-1 -1]
        [ 3  2]: []}

        This was the output when the output was still a list::
        [[[1 0]
        [0 2], [1 1]
        [0 2], [2 0]
        [0 1]], [], [], [], [], [], [], [], [], [], [], [[ 1 -1]
        [ 0  2]]]

        The output the original version of this file claimed is the
        following, but this disagrees with what we get, and with the
        .sage version (which agree with each other)::
        [[[1 0]
        [0 2], [1 1]
        [0 2], [2 0]
        [0 1]], [], [], [], [], [], [[ 1 -1]
        [ 0  2]], [], [], [], [], []]

        """
        N = self.level()

        ans = {}
        for h in self:
            ans[h] = []
        # this will be the dictionary D above enumerated by coset reps

        #  This loop will run thru the ell+1 (or ell) matrices
        #  defining T_ell of the form [1, a, 0, ell] and carry out the
        #  computation described above.
        #  -------------------------------------
        for a in range(ell + 1):
           if (a < ell) or (N % ell != 0):
               # if the level is not prime to ell the matrix [ell, 0, 0, 1] is avoided.
               gamma = basic_hecke_matrix(a, ell)
               t = gamma * gen
               #  In the notation above this is gam_a * D_m
               v = unimod_matrices_from_infty(t[0, 0], t[1, 0]) + unimod_matrices_to_infty(t[0, 1], t[1, 1])
               #  This expresses t as a sum of unimodular divisors

               # This loop runs over each such unimodular divisor
               # ------------------------------------------------
               for A in v:
                   #  B is the coset rep equivalent to A
                   B = self.equivalent_rep(A)
                   #  C equals A^(-1).
                   C = A._invert_unit()
                   #  gaminv = B*A^(-1)
                   gaminv = B * C
                   #  The matrix gaminv * gamma is added to our list in the j-th slot
                   #  (as described above)
                   tmp = M2Z(gaminv * gamma)
                   tmp.set_immutable()
                   ans[B].append(tmp)

        return ans

def basic_hecke_matrix(a, ell):
    """
    Returns the matrix [1, a, 0, ell] (if a<ell) and [ell, 0, 0, 1] if a>=ell

    INPUT:

    - `a` -- an integer or Infinity
    - ``ell`` -- a prime

    OUTPUT:

    - a 2 x 2 matrix of determinant ell

    EXAMPLES:

        sage: sage.modular.pollack_stevens.manin_map.basic_hecke_matrix(0, 7)
        [1 0]
        [0 7]
        sage: sage.modular.pollack_stevens.manin_map.basic_hecke_matrix(5, 7)
        [1 5]
        [0 7]
        sage: sage.modular.pollack_stevens.manin_map.basic_hecke_matrix(7, 7)
        [7 0]
        [0 1]
        sage: sage.modular.pollack_stevens.manin_map.basic_hecke_matrix(19, 7)
        [7 0]
        [0 1]
    """
    # TODO: probably a bottleneck.
    if a < ell:
        return M2Z([1, a, 0, ell])
    else:
        return M2Z([ell, 0, 0, 1])

def unimod_matrices_to_infty(r, s):
    """
    Returns a list of matrices whose associated unimodular paths
    connect 0 to r/s.  This is Manin's continued fraction trick, which
    gives an expression {0,r/s} = {0,oo} + ... + {a,b} + ... + {*,r/s},
    where each {a,b} is the image of {0,oo} under a matrix in SL_2(ZZ).

    INPUT:

    - `r`, `s` -- rational numbers

    OUTPUT:

    - a list of matrices in `SL_2(\ZZ)`

    EXAMPLES::

        sage: v = sage.modular.pollack_stevens.manin_map.unimod_matrices_to_infty(19,23); v
        [
        [1 0]  [ 0  1]  [1 4]  [-4  5]  [ 5 19]
        [0 1], [-1  1], [1 5], [-5  6], [ 6 23]
        ]
        sage: [a.det() for a in v]
        [1, 1, 1, 1, 1]
    """
    if s == 0:
        return []
    # the function contfrac_q in
    # https://github.com/williamstein/psage/blob/master/psage/modform/rational/modular_symbol_map.pyx
    # is very, very relevant to massively optimizing this.
    L = convergents(r / s)
    # Computes the continued fraction convergents of r/s
    v = [M2Z([1, L[0].numerator(), 0, L[0].denominator()])]
    # Initializes the list of matrices
    for j in range(0, len(L)-1):
        a = L[j].numerator()
        c = L[j].denominator()
        b = L[j + 1].numerator()
        d = L[j + 1].denominator()
        v.append(M2Z([(-1)**(j + 1) * a, b, (-1)**(j + 1) * c, d]))
        # The matrix connecting two consecutive convergents is added on
    return v


def unimod_matrices_from_infty(r, s):
    """
    Returns a list of matrices whose associated unimodular paths
    connect 0 to r/s.  This is Manin's continued fraction trick, which
    gives an expression {oo,r/s} = {oo,0} + ... + {a,b} + ... + {*,r/s},
    where each {a,b} is the image of {0,oo} under a matrix in SL_2(ZZ).

    INPUT:

    - `r`, `s` -- rational numbers

    OUTPUT:

    - a list of SL_2(Z) matrices

    EXAMPLES:

        sage: v = sage.modular.pollack_stevens.manin_map.unimod_matrices_from_infty(19,23); v
        [
        [ 0  1]  [-1  0]  [-4  1]  [-5 -4]  [-19   5]
        [-1  0], [-1 -1], [-5  1], [-6 -5], [-23   6]
        ]
        sage: [a.det() for a in v]
        [1, 1, 1, 1, 1]
    """
    if s != 0:
        L = convergents(r / s)
        # Computes the continued fraction convergents of r/s
        v = [M2Z([-L[0].numerator(), 1, -L[0].denominator(), 0])]
        # Initializes the list of matrices
        # the function contfrac_q in https://github.com/williamstein/psage/blob/master/psage/modform/rational/modular_symbol_map.pyx
        # is very, very relevant to massively optimizing this.
        for j in range(0, len(L) - 1):
            a = L[j].numerator()
            c = L[j].denominator()
            b = L[j + 1].numerator()
            d = L[j + 1].denominator()
            v.append(M2Z([-b, (-1)**(j + 1) * a, -d, (-1)**(j + 1) * c]))
            # The matrix connecting two consecutive convergents is added on
        return v
    else:
        return []
