{-# LANGUAGE NoImplicitPrelude, TypeSynonymInstances, RankNTypes, StandaloneDeriving,
ScopedTypeVariables, TypeApplications, GADTs, TypeOperators, 
FlexibleInstances, AllowAmbiguousTypes, UndecidableInstances #-}
module Lib
    where
import Numeric.LinearAlgebra
import Control.Category
import Prelude hiding ((.), id)
someFunc :: IO ()
someFunc = putStrLn "someFunc"

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
  
instance (Enum a, Enum b, Bounded a) => Enum (Either a b) where
    toEnum x | x > (fromEnum (maxBound :: a)) = Left (toEnum x)
             | otherwise = Right $ toEnum (x - 1 - (fromEnum (maxBound :: a)))
    fromEnum (Left a) =  fromEnum a
    fromEnum (Right b) =  (fromEnum (maxBound :: a)) + 1 + (fromEnum b) 

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