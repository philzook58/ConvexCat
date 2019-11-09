{-# LANGUAGE ConstraintKinds, ScopedTypeVariables, TypeApplications, AllowAmbiguousTypes, NoImplicitPrelude,
    GeneralizedNewtypeDeriving, FlexibleContexts
#-}
module LinRel where

import Numeric.LinearAlgebra
import Prelude hiding ((<>))
import Debug.Trace
type BEnum a = (Enum a, Bounded a) 
enumAll :: (BEnum a) => [a] -- What about Void?
enumAll = [minBound .. maxBound]




data Void = Void' | Void'' -- SORRY MOM.
instance Enum Void where
  fromEnum Void' =  1
  fromEnum Void'' = 0
  toEnum 1 = Void'
  toEnum 0 = Void''
instance Bounded Void where
  maxBound = Void''
  minBound = Void'

instance (Enum a, Enum b, Bounded a) => Enum (Either a b) where
  fromEnum (Left a) = fromEnum a
  fromEnum (Right b) = fromEnum b + fromEnum (maxBound @a) + 1
  toEnum n | n <= ((fromEnum (maxBound @a)) - (fromEnum (minBound @a))) = Left (toEnum n)
           | otherwise = Right (toEnum (n - fromEnum (maxBound @a))) 


instance (Bounded a, Bounded b) => Bounded (Either a b) where
  maxBound = Right maxBound
  minBound = Left minBound

card :: forall a. (BEnum a) => Int
card = (fromEnum (maxBound @a)) - (fromEnum (minBound @a)) + 1

-- HLinRel holds A x = b constraint
data HLinRel a b = HLinRel (Matrix Double) (Vector Double) deriving Show

-- x = A l + b. Generator constraint. 
data VLinRel a b = VLinRel (Matrix Double) (Vector Double) deriving Show


-- f(x) = xQx + bx
data QuadOp a b = QuadOp (Matrix Double) (Vector Double)
{-
qid :: QuadOp a a
qid = ? -- a = a ? Need Relations too.

qcompose :: QuadOp b c -> QuadOp a b -> QuadOp a c
qcompose (QuadOp q c) (QuadOp q' c') = QuadOp q'' c'' where
        where 
                ca = card @a
                cb = card @b 
                cc = card @c
                a = subMatrix (cb, cb) q
                b = subMatrix (cb,cc) q
                c = subMatrix (cc, cb) q
                d = subMatrix (cc,cc) q
                a' = subMatrix (ca, ca) q'
                b' = subMatrix (ca, cb) q'
                c' = subMatrix (cb, ca  q'
                d' = subMatrix (cb, cb) q'
                m = - (a + d')
                q'' = fromBlocks [[a' + b' <> m </> c'   , b' <> m </> c ],  -- can memoize some of this
                                  [b <> m </> c' , d - ]]
                
                a'' =  a' - b' <> m </> c'
                [v1', v2'] = takesV [ca,cb] c'
                [v1,  v2] = takesV [cb, cc] c
                v3 = (v2' + v1) 
                c'' = vJoin [v1' + m </> c #>  v3, ]

        -}
    {-
    
    Any cost outside the constraint space is irrelevant.
    Q = V m V where m is fullrank.
    c = Vc also

    minimal sets
    -}    
        
        -- break into pieces, form schur complement.

-- if A x = b then x is in the nullspace + a vector b' solves the equation
h2v :: HLinRel a b -> VLinRel a b
h2v (HLinRel a b) = VLinRel a' b' where
        b' = a <\> b -- least squares solution
        a' = nullspace a

-- if x = A l + b, then A' . x = A' A l + A' b = A' b because A' A = 0
v2h :: VLinRel a b -> HLinRel a b
v2h (VLinRel a' b') = HLinRel a b where
        b = a #> b' -- matrix multiply
        a = tr $ nullspace (tr a') -- orthogonal space to range of a.

hid :: forall a. BEnum a => HLinRel a a
hid =  HLinRel (i ||| (- i)) (vzero s) where 
                            s = card @a
                            i = ident s

vzero :: Konst Double d c => d -> c Double
vzero = konst 0

hcompose :: forall a b c. (BEnum a, BEnum b, BEnum c) => HLinRel b c -> HLinRel a b -> HLinRel a c
hcompose (HLinRel m b) (HLinRel m' b') = let a'' = fromBlocks [[       ma',           mb' ,    0       ],
                                                               [         0 ,    mb,        mc          ]] in
                                         let b'' = vjoin [b', b] in 
                                         let (VLinRel q p) = h2v (HLinRel a'' b'') in -- kind of a misuse
                                         let q' = (takeRows ca q)  -- drop rows belonging to @b
                                                       === 
                                                  (dropRows (ca + cb) q) in
                                         let [x,y,z] =  takesV [ca,cb,cc] p in
                                         let p'=  vjoin [x,z] in -- rebuild without rows for @b
                                         v2h (VLinRel q' p') -- reconstruct HLinRel
                                       where 
                                           ca = card @a
                                           cb = card @b 
                                           cc = card @c
                                           sb = size b -- number of constraints in first relation
                                           sb' = size b' -- number of constraints in second relation
                                           ma' = takeColumns ca m'
                                           mb' = dropColumns ca m'
                                           mb = takeColumns cb m
                                           mc = dropColumns cb m

(<<<) :: forall a b c. (BEnum a, BEnum b, BEnum c) => HLinRel b c -> HLinRel a b -> HLinRel a c
(<<<) = hcompose
-- stack the constraints
hmeet :: HLinRel a b -> HLinRel a b -> HLinRel a b
hmeet (HLinRel a b) (HLinRel a' b') = HLinRel (a === a') (vjoin [b,b'])




{- If they don't meet are we still ok? 

I am not sure. Might be weird corner cases?

-}

hjoin :: HLinRel a b -> HLinRel a b -> HLinRel a b
hjoin v w = v2h $ vjoin' (h2v v) (h2v w)

-- hmatrix took vjoin from me :(
-- joining means combining generators and adding a new generator
-- Closed under affine combination l * x1 + (1 - l) * x2 
vjoin' :: VLinRel a b -> VLinRel a b -> VLinRel a b
vjoin' (VLinRel a b) (VLinRel a' b') = VLinRel (a ||| a' ||| (asColumn (b - b'))) b

-- no constraints, everything
-- trivially true
htop :: forall a b. (BEnum a, BEnum b) => HLinRel a b 
htop = HLinRel (vzero (1,ca + cb)) (konst 0 1) where 
                                      ca = card @a
                                      cb = card @b 
{-
hbottom :: forall a b. (BEnum a, BEnum b) => HLinRel a b 
hbottom = HLinRel (vzero (1,ca + cb)) (konst 1 1) where 
                                      ca = card @a
                                      cb = card @b       
                                      -}                                
-- all the constraints! Only the origin.
-- no. it should be the empty set. Impossible to satisfy.
-- 0 x = 1 is impossible
-- not gonna play nice.
{-
hbottom :: forall a b. (BEnum a, BEnum b) => HLinRel a b 
hbottom = HLinRel (ident (ca + cb)) (konst 0 (ca + cb)) where 
                                    ca = card @a
                                    cb = card @b  
  -}                            

hconverse :: forall a b. (BEnum a, BEnum b) => HLinRel a b -> HLinRel b a 
hconverse (HLinRel a b) = HLinRel ( (dropColumns ca a) |||  (takeColumns ca a)) b where 
    ca = card @a
    cb = card @b  

    -- this is numerically unacceptable
-- forall l. A' ( A l + b) == b'
vhsub :: VLinRel a b -> HLinRel a b -> Bool
vhsub (VLinRel a b) (HLinRel a' b') = (naa' <=  1e-10 * (norm_2 a') * (norm_2 a)  ) && ((norm_2 ((a' #> b) - b')) <= 1e-10 * (norm_2 b')  ) where
          naa' = norm_2 (a' <> a)

hsub :: HLinRel a b -> HLinRel a b -> Bool
hsub h1 h2 = vhsub (h2v h1) h2

heq :: HLinRel a b -> HLinRel a b -> Bool
heq a b = (hsub a b) && (hsub b a)


instance Ord (HLinRel a b) where
  (<=) = hsub
  (>=) = flip hsub 




instance Eq (HLinRel a b) where
  (==) = heq


-- I can't do this right?
-- hcomplement :: HLinRel a b -> HLinRel a b
-- hcomplement  

hpar :: HLinRel a b -> HLinRel c d -> HLinRel (Either a c) (Either b d)
hpar (HLinRel mab v) (HLinRel mcd v') = HLinRel (fromBlocks [ [mab, 0], [0 , mcd]]) (vjoin [v, v']) where



{-
        -- Void is unit for Either.
-- void has no inhabitants.... This is a bad boy.


hcup :: HLinRel Void (Either a a)
hcap :: HLinRel (Either a a) Void




-}
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


hfan ::  forall a b c. BEnum a => HLinRel a b -> HLinRel a c -> HLinRel a (Either b c)
hfan (HLinRel m v) (HLinRel m' v') = HLinRel (fromBlocks [ [ma, mb, 0], [ma', 0, mc']]) (vjoin [v,v']) where
        ca = card @a
        ma = takeColumns ca m 
        mb = dropColumns ca m 
        ma' = takeColumns ca m' 
        mc' = dropColumns ca m' 


hdump :: HLinRel a Void
hdump = HLinRel 0 0

hlabsorb :: HLinRel a b -> HLinRel (Either Void a) b
hlabsorb (HLinRel m v) = (HLinRel m v)

htrans :: HLinRel a (Either b c) -> HLinRel (Either a b) c 
htrans (HLinRel m v) = HLinRel m v

hswap :: forall a b. (BEnum a, BEnum b) => HLinRel (Either a b) (Either b a)
hswap = HLinRel (fromBlocks [[ia ,0,0 ,-ia], [0, ib,-ib,0]]) (konst 0 (ca + cb)) where 
        ca = card @a
        cb = card @b  
        ia = ident ca
        ib = ident cb


hsum :: forall a. BEnum a => HLinRel (Either a a) a
hsum = HLinRel ( i ||| i ||| - i ) (konst 0 ca)  where 
        ca = card @a 
        i= ident ca

hdup :: forall a. BEnum a => HLinRel a (Either a a)
hdup = HLinRel (fromBlocks [[i, -i,0 ], [i, 0, -i]]) (konst 0 (ca + ca))  where 
        ca = card @a 
        i= ident ca

-- hcup :: forall a. BEnum a => HLinRel Void (Either a a)
-- hcup = -- or mainulate hid

-- hcap :: forall a. BEnum a => HLinRel (Either a a) Void

-- smart constructors
hLinRel :: forall a b. (BEnum a, BEnum b) => Matrix Double -> Vector Double -> Maybe (HLinRel a b) 
hLinRel m v | cols m == (ca + cb) && (size v == rows m)  = Just (HLinRel m v)
            |  otherwise = Nothing  where 
                 ca = card @a
                 cb = card @b  

-- a 2d space at every wire or current and voltage.
data IV = I | V deriving (Show, Enum, Bounded, Eq, Ord)


resistor :: Double -> HLinRel IV IV
resistor r = HLinRel ( (2><4)  [ 1,0,-1,   0,
                                 r, 1, 0, -1]) (konst 0 2)  

bridge :: Double -> HLinRel (Either IV IV) (Either IV IV)
bridge r = HLinRel (  (4><8) [ 1,0, 1,  0, -1, 0, -1,  0, -- current conservation
                               0, 1, 0, 0, 0, -1 , 0,  0, --voltage maintained
                               0, 0, 0, 1, 0,  0,  0, -1, -- voltage maintained
                               r, 1, 0,-1, -r,  0,  0, 0  ]) (konst 0 4)  

{- Legendre transformations in thermo are for open systems. SOmething to that -}
{- Dependent sources. Well, these are goddamn cheating. -}

{-

A wire in a circuit has a potential (questionably) and a current running through it
So our wires should have both of these variables.

-}
newtype VProbe = VProbe () deriving (Enum, Bounded, Show, Eq, Ord)
vprobe :: HLinRel IV VProbe
vprobe = HLinRel ( (2><3)  [1,0,0,
                            0,1,-1]) (konst 0 2)                

vsource :: Double -> HLinRel IV IV
vsource v = HLinRel ( (2><4) [ 1,0,-1,   0,
                               0, 1, 0, -1]) (fromList [0,v])  

isource :: Double -> HLinRel IV IV
isource i = HLinRel ( (2><4) [  1,0, -1,   0 , -- current conservation
                                 1, 0, 0,  0]) (fromList [0,i])  


-- the currents add, but the voltages dup. sum and dup are dual
-- Or should it be |--|   a parallel short?
-- Ad then we could open circuit one of them and absorb the Void
-- to derive this
cmerge :: HLinRel (Either IV IV) IV
cmerge = HLinRel ( (3><4)  [1, 0, 1, 0, -1, 0,
                            0,1,0,0,0 ,  -1  ,
                            0,0,0,1, 0, -1])  (konst 0 3)

open :: HLinRel IV Void
open = HLinRel ( (1><2) [1,0]) (konst 0 1)


cap :: HLinRel  (Either IV IV) Void
cap  = hcompose open cmerge

cup :: HLinRel Void (Either IV IV)
cup = hconverse cap

ground :: HLinRel IV Void
ground = HLinRel ( (1><2) [ 0 , 1 ]) (vzero 1) 

-- resistors in parallel.

ex1 = hcompose (bridge 10) (bridge 10)
ex2 = hcompose (resistor 10) (resistor 30) -- resistors in series.
r20 :: HLinRel IV IV
r20 = resistor 20

divider :: Double -> Double -> HLinRel (Either IV IV) (Either IV IV)
divider r1 r2 = hcompose (bridge r2) (hpar (resistor r1) hid) 



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