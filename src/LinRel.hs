{-# LANGUAGE ConstraintKinds, ScopedTypeVariables, TypeApplications, AllowAmbiguousTypes, NoImplicitPrelude
#-}
module LinRel where

import Numeric.LinearAlgebra
import Prelude hiding ((<>))
type BEnum a = (Enum a, Bounded a) 
enumAll :: (BEnum a) => [a]
enumAll = [minBound .. maxBound]

card :: forall a. (BEnum a) => Int
card = (fromEnum (maxBound @a)) - (fromEnum (minBound @a)) + 1

-- HLinRel holds A x = b constraint
data HLinRel a b = HLinRel (Matrix Double) (Vector Double)

-- x = A l + b. Generator constraint. 
data VLinRel a b = VLinRel (Matrix Double) (Vector Double)


h2v :: HLinRel a b -> VLinRel a b
h2v (HLinRel a b) = VLinRel a' b' where
        b' = a <\> b -- leasty squares solution
        a' = nullspace a

-- is x = A l + b, then a' . x = a' . a l + a' b = a' b because a' . a = 0
v2h :: VLinRel a b -> HLinRel a b
v2h (VLinRel a b) = HLinRel a' b' where
        b' = a' #> b -- matrix multiply
        a' = nullspace (tr a) -- orthogonal space to range of a.

lid :: forall a. BEnum a => HLinRel a a
lid =  HLinRel (i ||| (- i)) (konst 0 (2 * s)) where 
                            s = card @a
                            i = ident s

vzero = konst 0

hcompose :: forall a b c. (BEnum a, BEnum b, BEnum c) => HLinRel b c -> HLinRel a b -> HLinRel a c
hcompose (HLinRel a b) (HLinRel a' b') = let a'' = fromBlocks [[ a , vzero ( sb , cc)] , [vzero ( sb' , ca)  , a' ]  ] in
                                         let b'' = vjoin [b, b'] in 
                                         let (VLinRel q p) = h2v (HLinRel a'' b'') in -- kind of a misuse
                                         let q' = (takeRows ca q) === (flipud (takeRows cc (flipud q))) in
                                         let [x,y,z] =  takesV [ca,cb,cc] p in
                                         let p'=  vjoin [x,z] in
                                         v2h (VLinRel q' p') 
                                       where 
                                           ca = card @a
                                           cb = card @b 
                                           cc = card @c
                                           sb   = size b
                                           sb'  = size b'

hmeet :: HLinRel a b -> HLinRel a b -> HLinRel a b
hmeet (HLinRel a b) (HLinRel a' b') = HLinRel (a === a') (vjoin [b,b'])
{- If they don't meet are we still ok? -}

hjoin :: HLinRel a b -> HLinRel a b -> HLinRel a b
hjoin v w = v2h $ vjoin' (h2v v) (h2v w)
vjoin' :: VLinRel a b -> VLinRel a b -> VLinRel a b
vjoin' (VLinRel a b) (VLinRel a' b') = VLinRel (a ||| a' ||| (asColumn (b - b'))) b

-- no constraints, everything
htop :: forall a b. (BEnum a, BEnum b) => HLinRel a b 
htop = HLinRel (vzero (0,ca + cb)) (konst 0 0) where 
                                      ca = card @a
                                      cb = card @b 
hbottom :: forall a b. (BEnum a, BEnum b) => HLinRel a b 
hbottom = HLinRel (ident (ca + cb)) (konst 0 (ca + cb)) where 
                                    ca = card @a
                                    cb = card @b  
                              

hconverse :: forall a b. (BEnum a, BEnum b) => HLinRel a b -> HLinRel b a 
hconverse (HLinRel a b) = HLinRel ( (dropColumns ca a) |||  (takeColumns ca a)) b where 
    ca = card @a
    cb = card @b  

    -- this is numerically unacceptable
-- forall l. A' ( A l + b) == b'
vhsub :: VLinRel a b -> HLinRel a b -> Bool
vhsub (VLinRel a b) (HLinRel a' b') = ((a' <> a) == 0) && (a' #> b == b')

hsub :: HLinRel a b -> HLinRel a b -> Bool
hsub h1 h2 = vhsub (h2v h1) h2

-- I can't do this right?
-- hcomplement :: HLinRel a b -> HLinRel a b
-- hcomplement  

hpar :: HLinRel a b -> HLinRel c d -> HLinRel (Either a c) (Either b d)
hpar (HLinRel mab v) (HLinRel mcd v') = HLinRel (fromBlocks [ [mab, 0], [0 , mcd]]) (vjoin [v, v']) where

hleft :: forall a b. (BEnum a, BEnum b) => HLinRel a (Either a b)
hleft = HLinRel ( i ||| (- i) ||| (konst 0 (ca,cb))) (konst 0 ca) where 
    ca = card @a
    cb = card @b  
    i = ident ca

hright :: forall a b. (BEnum a, BEnum b) => HLinRel b (Either a b)
hright = HLinRel ( i ||| (konst 0 (cb,ca)) ||| (- i) ) (konst 0 cb) where 
    ca = card @a
    cb = card @b  
    i = ident cb



-- smart constructors
hLinRel :: forall a b. (BEnum a, BEnum b) => Matrix Double -> Vector Double -> Maybe (HLinRel a b) 
hLinRel m v | cols m == (ca + cb) && (size v == rows m)  = Just (HLinRel m v)
            |  otherwise = Nothing  where 
                 ca = card @a
                 cb = card @b  



                 {-
                 
                 Is there a reasonable intepretation of kron?

                 -}
{-
everything can be definedc inefficiently via v2s and h2v functions

right division

-}

    {-
Call them affine relations

Join and meet aren't union and intersection.
They are the affine closure of union and intersection.




Linear has some niceness.
Homgeonous coordinates usually do.
For clarity and familiaryt I have chosebn not to do it this way
Or maybe I will do it?

par


-}


{-

import numpy as np


def meet(a,b):
    pass
def compose(a,b): # a after b
    assert(a.inN == b.outN)

    combo = np.block([[a.constraints, np.zeros((a.constraints.shape[0] , b.inN) )     ],
                      [np.zeros((b.constraints.shape[0] , a.outN)) ,      b.constraints]])
    print(combo)
    gens = LinRel(a.inN + a.outN, b.outN, gens=combo).gens
    print("gens",gens) 
    gens = np.vstack((gens[:a.outN, :], gens[-b.inN: , :]) )
    print(gens)
    return LinRel(a.outN, b.inN, gens=gens)
def top(outN, inN):
    return LinRel(outN, inN, constraints = np.array([[]]))
def bottom(outN, inN):
    return LinRel(outN, inN, gens = np.array([[]]))

def converse(a):
    return LinRel(a.inN, a.outN, np.hstack((a.constraints[:, a.inN:], a.constraints[:, :a.inN])))
def complement(a):
    return LinRel(a.outN, a.inN, constraints = a.gens.T.conj())

def right_div(a,b):
    pass # return complement( compose(a, complement(b)) ) # something like this
def inclusion(a,b):
    
    s = a.constraints @ b.gens
    np.all(a.constraints @ b.gens <= tol )

    if rcond is None:
        rcond = np.finfo(s.dtype).eps * max(max(a.shape), max(b.shape))
        tol = max(np.amax(a), np.amax(b))
        tol = np.amax(s) * rcond
def fromMat(mat):
    (outN, inN) = mat.shape
    return LinRel(inN, outN,constraints = np.hstack((mat,-np.eye(outN))))
def id(N):
    return fromMat(np.eye(N))
# make 0,1 first index for in/out? Oh, but then in out have to be same size.
# A[0, ...] @ x + A[1, ...] @ y = 0   
# Then I can form the kron of linear relations. 
# store sperate A B matrices? A @in + B @ out
class LinRel():
    def __init__(self, outN, inN, constraints = None, gens = None, rcond=None):
        #assert(inN <= constraints.shape[1])
        self.inN = inN
        self.outN = outN
        if constraints is not None: #baiscally scipy.linalg.null_space
            u, s, vh = np.linalg.svd(constraints, full_matrices=True)
            M, N = u.shape[0], vh.shape[1]
            if rcond is None:
                rcond = np.finfo(s.dtype).eps * max(M, N)
            tol = np.amax(s) * rcond
            num = np.sum(s > tol, dtype=int)
            self.gens = vh[num:,:].T.conj()
            self.constraints = vh[:num,:]
        if gens is not None: #basically scipy.linalg.orth
            u, s, vh = np.linalg.svd(gens, full_matrices=True)
            M, N = u.shape[0], vh.shape[1]
            if rcond is None:
                rcond = np.finfo(s.dtype).eps * max(M, N)
            tol = np.amax(s) * rcond
            num = np.sum(s > tol, dtype=int)
            self.gens = u[:, :num]
            self.constraints = u[:, num:].T.conj()

    def shape(self):
        return (self.outN, self.inN)
    def size(self):
        return self.outN + self.inN
    # operator overloadings
    def __matmul__(a,b):
        return compose(a,b)
    def __invert__(a): # ~
        return complement(a)
    def __or__(a,b): # an argument could be made for + and *
        return join(a,b)
    def __and__(a,b):
        return meet(a,b)
    def __sub__(a,b):
        return  (a) & (-b)
    def __le__(a,b): # Are the others automatic?
        return inclusion(a,b)
    def __str__(self):
        return " Constraints: \n%s, \nGens:\n %s\n" % (str(self.constraints), str(self.gens))




ex = LinRel(1,2, np.array([[3,4,0]]))
e2 = LinRel(2,1, np.array([[3,4,0]]))

assert(np.all(LinRel(1,2, ex.constraints).constraints == ex.constraints) )
print( ex.constraints @ ex.gens)
assert(np.all( np.abs(ex.constraints @ ex.gens) <= 1e-15) )
print(ex @ e2)
print(e2 @ ex)
print(e2 @ id(e2.inN))

'''
Quadratic optimization can be mixed in.
Quad && LinRel

AffineRel = maintain homgenous coord, or always insert 1 -1 keeping homoeg coord
+ discrete? Maintain a couple copies of 


'''


Linear relations
hrep - Ax = 0
or
vrep - y = sum x_i

hrep <-> vrep = row echelon

in and out variables. In and out subspaces.
in : [] - list of indices
out : []

in is a vrep of input space.

in = projection/injection matrix n x d
out = projection/injection matrix. (d-n) x d
in * out = 0. orthogonal

auxiliary variables allowed

compose relations
in1 out1
in2 out2

in = d x n
stack out1 and in2 into matrix. with them equal.

np.hstack(A1, out1 - in2, A2)
in = no.hatck(in, zeores)
in = no.hatck(zeores, out)

drect sum as monoidal product
block([ A, 0  ],
      [ 0, A  ])
in = [in1, in2]
out = [out1, out2]

converse = flip in out

meet = combine the two constraint matrices
join = convert? combine in out?

1-d as unit object
<= is subspace ordering

negation = orthogonalization
division => 


a linear problem is a linear relation of 1-d.

use in1 and out2 as new in/out

fan 
snd(30,10) = project bottom 10 = idetnity matrix stacked.
fst(30, 10) project top 10
id(20)
id(n) = LinRel(np.zeros(0,2*n), [zeros, eye], [eye, zero] ) 
id(n) = LinRel((0,n), eye(n), eye(n)  )
        LinRel [I, -I], [I, 0], [0,I]

"internal" space
class LinRel():
    def __init__(A, in, out):

svd(in * A , smallest )
1 - A*A

Ax + By = 0

vstack adds constraints
hstack adds variables

def idRel(n):
  return LinRel(sparse.lil_matrix((0,n)), n, n)
def mapRel(A):
  (r,c) = A.shape
  newA = sparse.hstack(A, - sparse.eye(r))
  return LinRel(newA, r, c) 

class LinRelV():
class LinRelH():

class LinRel():
  def init(self,A, in, out):
    self.A = A
    self.in = in
    self.out = out
  def compose(self, b):
    assert(self.out == b.in, "Shapes must match")
    
    i = sparse.eye(b.in)
    cons = sparse.hstack([0, i, -i, 0])
    ina = self.A[:,:self.in]
    auxa = self.A[:,self.in:-self.out]
    outa = self.A[:,-self.out:]
    inb = b.A[:,:b.in]
    auxb = b.A[:,b.in:-b.out]
    outb = b.A[:,-b.out:]

    newA = sparse.bmat(  [[ina, auxa, outa, 0, 0,    0],
                          [0  , 0   , i,   -i, 0 ,   0]
                          [0,  0,   0,    inb, auxb,outb]])
    LinRel(newA, self.in, b.out)
  def meet(self,b):
    #hmm. I suppose we acutlaly should split the thing apart again
    assert(self.in == b.in)
    assert(self.out == b.out)
    assert()
    newA = sparse.vstack([self.A, b.A])
    return ()
  def complement(self):
    linalg.svd(self.A)
    return LinRel(get_nonzeroeigs)
  def __negate__(self):
    return self.complement()
  def rdiv(self):
  def transpose():
    self.converse()
  def T():
    self.converse()
  
  # the svd gives you most of what you need?
  def inclusion():
    x = linalg.nullspace(self.A)
    return b.A @ x == 0 #check that every generator is in. Makes sense. Except numrically is trash.
  def __leq__(self):
    self.inclusion(b)

  

  def converse(rel):
      newA = np.hstack( [ rel.A[:,-rel.out: ] , rel.A[:,rel.in:-rel.out], rel.A[:,:rel.in ] ])
      return LinRel(newA, rel.out, rel.in)


compose(A)
linalg.nullspace(A)
range?



(cons1,d1) = self.A.shape
(cons2,d2) = b.A.shape
constrain = np.hstack 
newA = sparse.vstack 


using bigM I can encode the complement of a H-polytope 
but then what?
I do it again I guess?


complementation of relation ->
At least one must be inverted.

polytope inclusion
-> search for point not in B.
Or, do sandardinni encoding

Really It should generate new variables for an instantiation.
Snd should not reuse the same variables every time.
class PolyRel()
  invars = []
  constraints = [] # store >= values, not full constraints?
  outvars = []
  def __init__():
    all fresh variables
  def compose(self,b):
    PolyRel(self.constraints + self., invars = outvars )
  def complement():
    zs = []
    for c in constraints:
      z, constraints = reify(c)
      # z = cvx.Variable(1, boolean=True) # actually make same shape as c
      # c += c + M * z
    sum(zs) >= 1 # one or more constraints is disatisfied.
  def rdiv():

yeah. We should use a dsl, compile it, then encode it.
data Rel a b where
  Compose ::
  Complement ::
  Converse ::


relu

l1 >= 0
l2 >= 0
l1 <= M * z
l2 <= M * (1 - z)
x = lambda1 - lambda2
y = lambda1

Maybe insetado f subspaces, we should be thinking ellipses and svd.



-}