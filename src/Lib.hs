{-# LANGUAGE NoImplicitPrelude, TypeSynonymInstances, RankNTypes, StandaloneDeriving,
ScopedTypeVariables, TypeApplications, GADTs, TypeOperators, DeriveTraversable,
FlexibleInstances, AllowAmbiguousTypes, UndecidableInstances, TypeFamilies, LambdaCase #-}
module Lib
    where
import Numeric.LinearAlgebra
import Control.Category
import Prelude hiding ((.), id)
import Control.Arrow ((***))
someFunc :: IO ()
someFunc = putStrLn "someFunc"



{-
type V a = a -> Double
type D a = a -> Double
D (V a)
type V a  = ContT () Q a -- ContT b Q a = (a ->Q b) -> Q b. Has some way of dealing with a index. Some kind of
SelT -- takes a diag, unsummed?
-- (a -> Q b) -> Q a

newtype Dag a = Dag (a -> Q ())

class Dagger a where
    type Dual :: * -> * -- *?
    dag :: a -> Dual a
instance Dagger (D (Vec a)) where
    type Dual = Vec
    dag (D v) = join $ ContT ()
-- dagger needs to embed a search procedure?

instance Dagger (D (Vec Bool)) where
    type Dual = Vec Bool
    dag (D v) = join $ ContT ()
instance Dagger a b where
    dag :: a -> b
    dag' :: b -> a
instance Dagger k where
    type family Dual :: * -> *
    dag :: k a (Dual a) -- dag'?

-- also a kind of dot product
dag1 :: ((a -> Q ()) -> Q ()) -> (((a -> Q ()) -> Q ()) -> Q ()
dag1 f = 

-- it's a kind of dot product.
lift :: f Double -> (a -> Double) -> Double
lift w f = fold (+) $ (mapWithIndex \k v -> (f k) * v) w

-- tensor is a_i b_i, mutiply out without summing
tensor :: f a ->  f a  -> f a
tensor = intersectWith (*)

type V' a = (a -> Q ()) -> Q ()
tensor :: V' a -> V' a -> V' a

dag :: (V' Bool -> Q ()) -> V' Bool
dag f = \g -> True False

dag :: (Bool -> Q ()) -> V' Bool
dag :: (Bool -> Q ()) -> (Bool -> Q ()) -> Q ()
dag f g = (f True) (g True) + (f False) (g False)



dag :: (Bool -> Q ()) -> Q Bool
dag f = (f True) >>= true + (f False) >>= false

dag :: Q Bool -> (Bool -> Q ())
dag xs b = lookup b xs 
dag :: Eq a => Q a -> (a -> Q ()) -- maybe Ord
dag xs b = lookup b xs 

-- If we mixed Hedges select sat solver with the Q monad?
-- (Bool -> Q Bool) -> Q Bool ? the unsummed trace? diag? Hedgest hing would come out weirder.


true :: () -> Q Bool
true _ = pure True
false :: () -> Q Bool
false _ = pure False

-- the seelect function.
-- (b -> m ()) -> m b

-}


type M = Matrix Double

{-

data BMatrix a b = BMatrix M M M M | Id -- A B C D

instance Category BMatrix where
    id = Id
    Id . x = x
    x . Id = x
    (BMatrix a b c d) . (BMatrix e f g h) = BMatrix <\> where
            q = d + e
            qlu = lupacked q
            qinv = lusolve qlu
            a' = a - (b <> (qinv c))
            b' =  b <> (qinv f)
            c' = g <> (qinv c)
            d' = h - g <> (qinv f)

par :: BMatrix a b -> BMatrix c d -> BMatrix (a,c) (b,d)
par (Matrix a b c d) (BMatrix e f g h) = BMatrix (diagBlock [a,e]) (diagBlock [b,f]) (diagBlock [c,g])  (diagBlock [d,h]) 
dup = BMatrix a (a,a)
dup = BMatrix (diagBlock [ident,ident]) (ident ||| ident) (ident === ident) (diagBlock [ident,ident])
id = BMatrix ident ident ident ident
avg = 
-}
data QuadFun a b = QuadFun {mat :: M, vec :: Vector Double, c :: Double } deriving Show


--data QuadFun a b where
--   QuadFun :: (Bounded a, Bounded b, Enum a, Enum b) => M -> Vector Double -> Double -> QuadFun a b } 

-- This is kron product.
instance (Enum a, Enum b, Bounded a) => Enum (a,b) where
  toEnum x =  let (a, b) = quotRem x (fromEnum (maxBound :: a)) in
                (toEnum a   , toEnum b)
  fromEnum (a,b) =  (fromEnum (maxBound :: a)) * (fromEnum a) + (fromEnum b)
  {-
instance (Enum a, Enum b, Bounded a) => Enum (Either a b) where
    toEnum x | x > (fromEnum (maxBound :: a)) = Left (toEnum x)
             | otherwise = Right $ toEnum (x - 1 - (fromEnum (maxBound :: a)))
    fromEnum (Left a) =  fromEnum a
    fromEnum (Right b) =  (fromEnum (maxBound :: a)) + 1 + (fromEnum b) 
-}
-- deriving instance (Enum a, Enum b) => Enum (a,b)


count :: forall a. (Enum a, Bounded a) => Int
count = (fromEnum (maxBound @a)) - (fromEnum (minBound @a)) + 1

id ::forall a. (Enum a, Bounded a) => QuadFun a a
id = let s = 2 * (count @ a) in QuadFun (ident s) (konst 0 s) 0

--dup :: QuadFun a (a,a)
--dup  = QuadFun (ident

class (Bounded a, Enum a) => BEnum a where
instance (Bounded a, Enum a) => BEnum a where

t7 :: QuadFun () ()
t7 = QuadFun (konst 4 (2,2)) (konst 7 2) 3

t8 = (mat t7) <\> (vec t7)
t4 = par t7 t7
t9 = (mat t4) <\> (vec t4)

compose :: forall a b c. BEnum b => QuadFun b c -> QuadFun a b -> QuadFun a c
compose (QuadFun m' w' c') (QuadFun m w c) = QuadFun m'' w'' c'' where
   n = count @b
   k = rows m -- assuming square
   l = rows m'
   i = ident n
   q = konst 0 (k-n,n)  === i  -- I underneath zeros
   q' = -i === konst 0 (l-n,n)   -- -i above zeros
   m'' = fromBlocks [[m,    q, 0], 
                     [tr q ,0, tr q'], 
                     [0,    q', m']]
   w'' = vjoin [w, konst 0 n, w']
   c'' = c + c'

identOut :: forall a b. BEnum b => QuadFun a b -> M
identOut (QuadFun m w c) = konst 0 (k-n,n)  === i  where
    n = count @b
    k = rows m
    i = ident n
identIn :: forall a b. BEnum a => QuadFun a b -> M
identIn (QuadFun m w c) = -i === konst 0 (k-n,n)   where
    n = count @a
    k = rows m
    i = ident n
type a :+: b = Either a b
-- I kind of feel like our sign convention is flipped
par :: forall a b c d. (BEnum a, BEnum b, BEnum c, BEnum d) => QuadFun a c -> QuadFun b d -> QuadFun (a :+: b) (c :+: d)
par x@(QuadFun m w c) y@(QuadFun m' w' c') = QuadFun m'' w'' c'' where
    ia = ident (count @a)
    ib = ident (count @b)
    ia' = identIn x
    ib' = identIn y
    ia't = tr ia'
    ib't = tr ib'
    ic = - ident (count @c)
    id = - ident (count @d)
    ic' = identOut x
    id' = identOut y
    ic't = tr ic'
    id't = tr id'
    n = (count @a) + (count @b)
    n' = (count @c) + (count @d)
    --iab = fromBlocks [[0,0,ia,0], [0,0,0,ib], [ia, 0,0,0], [0,ib,0,0]]
    --iab' = fromBlocks [[tr (identIn x),0],  [0 , tr (identIn y)]]
    --zab = konst 0 (n,n)
    -- should be 2 ins + 2 langrange + 2 matr + 2lagarne + 2 outs = 10x10 block matrix
    m'' = fromBlocks [[0 ,0 , ia, 0, 0, 0, 0 ,0,0,0], 
                    [0, 0,  0, ib, 0, 0, 0, 0,0,0 ] , 
                    [ia, 0, 0, 0,  ia't, 0,0,0,0, 0],
                    [0, ib, 0, 0,  0,ib't,0, 0,0 ,0],  
                    [0, 0, ia',0,  m, 0 , ic',0,0,0],
                    [0, 0,  0, ib',0, m', 0,id',0,0], 
                    [0, 0,  0, 0, ic't,0 ,0,0,ic,0], 
                    [0, 0,  0, 0,  0,id't,0,0,0,id],
                    [0, 0,  0, 0,  0, 0,  ic,0,0,0],
                    [0, 0,  0, 0,  0, 0,  0,id,0,0]]
    w'' = vjoin [konst 0 (2*n), w,w', konst 0 (2*n')]
    c'' = c + c'

dup :: forall a. BEnum a => QuadFun a (a,a)
dup = QuadFun m 0 0 where -- 1 input, 2 lagrange mutipliers, and 2 outputs.
    ia = ident (count @a)
    m = fromBlocks [ [0, ia,ia, 0,0],
                     [ia, 0, 0, -ia,0],
                     [0, ia, 0, 0, -ia],
                     [0,-ia, 0, 0, 0 ],
                     [0, 0, -ia, 0,0]]

-- fst...?
-- for snd, we could just leave it alone
-- fuse?
-- swap.
data Void
fuse :: forall a. BEnum a => QuadFun (a,a) Void
fuse = QuadFun m 0 0 where -- 2 inputs, 1 lagrange multiplier.
    ia = ident (count @a)
    m = fromBlocks [[0,0,    ia],
                    [0 ,0 , -ia],
                    [ia, -ia, 0]]

-- The analog for this slicing for convex sets may be to slice the 
-- dimensions into in and out dimensions. Then 



data Cell a = Cell {phi :: Double, j :: Double, next :: a} deriving (Show, Traversable, Foldable, Functor) -- composition of cell is a statically sized vector. Derive applicative?



instance Applicative Cell where
    pure x = Cell 0 0 x
    (Cell phi1 j1 f) <*> (Cell phi2 j2 x) = Cell (phi1 + phi2) (j1 + j2) (f x) -- ? I doubt this makes any sense.
{-
data Cell a = Cell { phi :: a, j :: a}
type f :+: g = Product f g
gaussK :: Lens ((Cell :+: Cell) a) (Cell a)

gaussK :: -> Lens (Cell :*: f a) -- no Cell is very not this.

type Row = (Cell Cell Cell Cell Cell)


Lens (a, b, other) (b, other)


gaussK :: Lens 

data Cell2 f a b = Cell2 {phi :: a , j :: a, next :: f b} 
data Cell2 f a = Cell2 {phi :: f Double a , j ::f Double, next :: f a} 

-- zipA :: f a -> f b -> f (a,b) 
-- zipA x y = (,) <$> x <*> y 
-- monProd

fmap2 = fmap . fmap
fmap4 = fmap2 . fmap2
fmap8 = fmap4 . fmap4

zipA x y = (,) <$> x <*> y 
zip2 = zipA . (fmap zipA) -- hmm. Maybe we should be using Compose. we're going to need famp2. Compose will get us these instances for free.
zip4 = zip2 . (fmap2 zip2)
zip8 = zip4 . (fmap4 zip4)


parK :: Lens (f b) b -> Lens (g a) a -> Lens (f g a) a
parK :: forall b. Lens (f b) (f' b) -> Lens (g a) (g' a) -> Lens (f g a) (f' (g' a))
parK :: forall b. Lens (f g a) (f' g a) -> Lens (g a) (g' a) -> Lens (f g a) (f' (g' a))
parK = compose (fmap l2) l1




Do the y direction inductively

ydir :: Lens (Cell a, Cell a) (a,a)
ydir = \(Cell phi j x, Cell phi2 j2 y) -> ((x,y), \(x',y') -> (Cell phi j x, Cell phi2 j2 y)      )


ydir :: Lens (Cell a, Cell a) (a,a)
ydir = \(a,b) -> ((Cell phi j x) , f) = (gaussK a) in ((Cell phi j y) , g) = (gaussK b) in 

Lens ( (X (), (X (), a) (X (), a)
ydir = \(r1, (r2, z) ->  zip8    )    -- = gaussK r1 in  = gaussK r2 in 

Lens (f a) a -> Lens (f a) a -> Lens (f a, f a) a
Lens (Cell a) a -> Lens 

Lens x y x' y'

-}

{-
can i just use regular lenses? I need to set in kind of a weird way though.
-}

data Lens a b = Lens (a -> (b, b -> a)) 
-- newtype Lens a b = Lens forall r. (a -> (b -> (b -> a) -> r) -> r) -- cpsify for speed? van Laarhoven?

-- SLens s a b = SLens (a -> ST s (b, b -> ST s a)) mutable update
-- MLens a b = MLens (a -> m (b , b -> m a)) -- monadic lens
-- KLens a b = KLens (a -> (b, k b a)) ---- this is sort of what Conal wrote.
comp :: Lens b c -> Lens a b -> Lens a c
comp (Lens f) (Lens g) = Lens $ \x -> let (y, g') = g x in
                                      let (z, f') = f y in
                                      (z, g' . f')

-- Gauss Seidel of the standard 1 -2 1 matrix
-- j1 is the effective source incliuding the influence from phi values up in the stack
-- Double Cell construction is rather wonky
-- gaussK :: Vec (S (S n)) Double -> Vec (S n)
-- we could also go for the lagrange mutiplier interface.
-- can also push the rho value needed for stability up and down.

-- -2 phi1 + 1 phi2 = ~j1
-- -1 phi1 + -2 phi2 + ...? = j2

{-
We are moving the lower diagonal to the right hand side as the splitting.
That ends up to claculating an effective j based on previous values.
We then triangular solve the upper diagonal, which we are able to find the new diagonal element in terms of the lower values.

-}
-- it's kind of weird that we mutate j on the way down but restore it on the way up.
-- But we do need a way to access j. (Phi (Phi a), J a) (Phi )

gaussK :: Lens (Cell (Cell a)) (Cell a) -- we need the context Cell
gaussK = Lens $ \case (Cell phi1 j1 (Cell phi2 j2 y)) -> (Cell phi2 (j2 - phi1) y , \case (Cell phi2' j2' z) -> let j1' = j1 - phi2' in  -- moving the triangular upsolve to the right hand side
                                                                                                             Cell (- j1' / 2) j1 (Cell phi2' j2 z))

-- interface in , internals, interface out. was the previous way of talking about it. But then the internals grow, which is fine
-- compose :: Lens in internal out -> Lens in' internal' out' -> Lens in (internal, out, in', internal') out
-- some ind of 2 category? Enriched? The morphisms have this internal structure.
-- this is more a containing relationship. The larger context can be converted into the smaller context.                                                                                                             
-- 

g2 :: Lens (Cell (Cell (Cell a))) (Cell a)
g2 = gaussK `comp` gaussK
g4 = g2 `comp` g2
g8 = g4 `comp` g4
g16 = g8 `comp` g8
g32 = g16 `comp` g16

capZero :: Lens (Cell ()) () -- not even sure i really need this? does runGauss (f `comp` capZero) ~ runGauss f
capZero = Lens $ \case (Cell phi j _) -> ((), \_ -> Cell (- j / 2) j ())

-- runGauss ::  Lens a b -> a -> a -- removes the open ended nature of the thing.
runGauss ::  Lens a () -> a -> a -- this might make more sense as what you really want. This is some kind of cap operation.
runGauss (Lens f) x = let (y , trisolve)  = f x in trisolve y
startingVal :: Cell (Cell (Cell ()))
startingVal = (Cell 0 0 (Cell 0 1 (Cell 0 0 ())))

iters = iterate (runGauss (capZero `comp` g2)) startingVal
--iters' = iterate (runGauss g2) startingVal -- They are different. 

sol = ((3><3) [-2,1,0,
               1,-2,1,
               0,1,-2])  <\> (vector [0,1,0])

-- can do inner block solves also (via an iterative method perhaps)
-- this is reminsecent of a multiscale attack.
parC :: Lens a a' -> Lens b b' -> Lens (a,b) (a',b')
parC (Lens f) (Lens g) = Lens $ \(x,y) -> let (x', f') = f x in
                                          let (y', g') = g y in
                                          ((x',y') , f' *** g')

-- profucntor optics paper
data FunList a b t = Done t | More a (FunList a b (b -> t))
-- contents and fill
-- contents (s -> a^n)
-- fill s -> b^n -> t
-- s -> exists n. (a^n, a^n -> s) is a traversal. Could use a Vec. But we already have a Vec. Cell is a Vec
-- s -> exists n, (Vec n a, Vec n a -> s)
-- maybe we do need an applicative. We sort of need to zip together two rows to start going 2d.

-- huh. A block wise schur is kind of the monoidal ppoduct here. No. monoidal product is pure dsum.
-- actually composition? is kind of what plays the game of block wise stuff.
{-

if we want this to be not the case,
comp :: Lens (a,a') a'


schur :: Lens a a' -> (a -> b) -> (b -> a) -> Lens b b' -> Lens (a,b) (a',b')
schur (Lens a) b c (Lens d) = Lens $ \(x,y) -> let (x', a') = a x in
                                               let (y', d') = d x in
                                               (       )
-}
-- Lens (f (g a)) (g a)

                                               -- we can change this all into a continuation form
-- 

{-
par :: forall a b c d. (BEnum a, BEnum b, BEnum c, BEnum d) => BMatrix a b -> 
                         BMatrix c d -> BMatrix (a,c) (b,d)
par q@(QuadFun m w r) q'@(QuadFun m' w' r') = QuadFun q'' w'' r'' where
     a' = diagBlock [sliceA q, sliceA q']
     b' = diagBlock [sliceB q, sliceB q']
     c' = diagBlock [sliceC q, sliceC q']
     d' = diagBlock [sliceD q, sliceD q']
     q'' = fromBlocks [[a',b'], [c',d']]
     u' = vjoin [sliceU v, sliceU v']
     v' = vjoin [sliceV v, sliceV v']
     w'' = vjoin [u', v']
     r'' = r + r'





sliceA :: forall a b. (BEnum a, BEnum b) => QuadFun a b -> M
sliceA (QuadFun m v c) = subMatrix (0,0) (count @a, count @a)  m
sliceB :: forall a b. (BEnum a, BEnum b) => QuadFun a b -> M
sliceB (QuadFun m v c) = subMatrix (0,count @a) (count @a, count @b) m
sliceC :: forall a b. (BEnum a, BEnum b) => QuadFun a b -> M
sliceC (QuadFun m v c) = subMatrix (count @a,0) (count @b, count @a) m
sliceD :: forall a b. (BEnum a, BEnum b) => QuadFun a b -> M
sliceD (QuadFun m v c) = subMatrix (count @a,count @a) (count @b, count @b) m

sliceU :: forall a b. (BEnum a, BEnum b) => QuadFun a b -> Vector Double
sliceU (QuadFun m v c) = subVector 0 (count @a) v
sliceV :: forall a b. (BEnum a, BEnum b) => QuadFun a b -> Vector Double
sliceV (QuadFun m v c) = subVector (count @a) (count @b) v

n = count @a
m = count @b
q = rows m
p = q - m
a = subMatrix (0,0) (count @a, count @a)  m


a = subMatrix
--Other options
-- Keep a matrix with the end and the beginning are the 

-- Keep Constraints too?

-- what in the world
-- stack build --ghc-options /usr/lib/libiconv.dylib

-- size
-- rows
-- cols

compose :: QuadFun b c -> QuadFun a b -> QuadFun a c
compose (QuadFun m w c) (QuadFun m' w' c') = QuadFun m'' w'' c'' where
    n = (count @b)
    corner =  fromBlocks [[ident n, 0],[0, -(ident n)]]
    corner' = fromBlocks [[0,0],[corner,0]]
    m'' = fromBlocks [[m ,corner'], [tr corner', m']]  
    w'' = vjoin [w, const 0 n, w']
    c'' = c + c

par (QuadFun m w c) (QuadFun m' w' c') = 


type a :+: b = Either a b
data QuadFun x = QuadFun


((a,b,c),
 (d,e,f),
 (g,h,i)) = 


 data BMatrix' = BMatrix' M M M M M M M M M 
 compose (BMatrix' a' b' c' d' e' f' g' h' i') (BMatrix' a b c d e f g h i) =
    BMatrix' a'' b'' c'' d'' e'' f'' g'' h'' i'' where
    a'' = a
    b'' = b ||| c ||| 0
    c'' = 0
    d'' = d === g === 0
    e'' = fromBlocks [[e,f,0], [h, i + a', b'], [0, d', e']]
    f'' = 0 === c' === f'
    g'' =  0
    h'' = 0 ||| g' ||| h'
    i'' = i'
compose' (BVec u' v' w') (BVec u v w) = BVec u (vjoin [v ,w + u', v']) w'


par :: QuadFun a b -> QuadFun c d -> QuadFun (Either a c) (Either b d)
par (BMatrix' a' b' c' d' e' f' g' h' i') (BMatrix' a b c d e f g h i) =
   BMatrix' a'' b'' c'' d'' e'' f'' g'' h'' i'' where
   a'' = diagBlock [a',a]
   b'' = diagBlock [b',b]
   c'' = diagBlock [c',c]
   d'' = diagBlock [d',d]
   e'' = diagBlock [e',e]
   f'' = diagBlock [f',f]
   g'' = diagBlock [g',g]
   h'' = diagBlock [h',h]
   i'' = diagBlock [i',i]

-- If we're just recording the entire program, we might as well just record it in a DSL
-- rather than directly building the matrix
-- data QuadFunDSL = Par QuadFun
-- FreeCat + (Lit HMatrix)

-- iterative matrix solution proximal method?
-- 
-- gauss seidel
type VBlock = [Vector Double]
-- (V, V)
-- Maybe schur solve works. But with possible splitting
-- Maybe unsolvable schur
-- could damp a little if the diagonal isn't dominant enough 
-- (lusolve a (w - V <#| v', ) -- WHat I'm suggesting here is block Jacobi? No
-- No it isn't. Because we're passing back v' it is guass seidel.
-- in order to do block gauss seidel
-- type L = VBlock -> (VBlock, VBlock -> VBlock)

-- Writing in the monoidal categroy style makes parallelism more apparent.


-- [vjoin blocks] -> ([lesser blocks] , update [lesserblocks] -> [blocks]
-- \x:xs -> (a*x + xs, \b -> b )  

-- storage in the lens?
-- In my AD lens, I had all the weights in the input. This was ungainly
-- Maybe a choice monad? Everyone could have access to a global state
-- forall s. (ActualLens s a, a -> (b, db -> da))
-- compsoe (ActualLens s a, a -> (b, db -> da))
-- par 
-- ActualLens are functional references. Do they allow us to refiy sharing?
-- (ActualLens acc a,   )

-- alternative a -> State s (b, db -> State ds da)

-- Mutiple directions of composition. Horizontsl and vertical?
-- a 2-category?
-- ActualLens s a -> ActualLens s' s

-- Jules Hedges:
-- Lens sigma a b -> Lens sigm' b c -> Lens (sig,sig') a c
-- Lens :: Sigma -> (sigma,  ) -> Lens a b -- hides sigma, we'll never be able to get at it 

-}