module Mixed where
import Linear.V4
import Data.Foldable
{-



branch and bound - store a best yet
Have a hierarchy of relaxed problems - converting

V4 a b c d = V4 a b c d

V4' a b = V4' a b a

-- all possible choices of a and b
-- somehow connect these functor together
V1' a b = V1 a
V1' a b = V1 b

-- this also perhaps let us carry down tips? So it's not all bad.
V4 (Either Int Double)
data VarType = Int | Double
data V4' = (V4 VarType, Vector Double, Vector Int) -- not intrisinacally safe, but makes sense.
data Vec a (Either Int Double) = [(a, Either Int Double)] -- free vector style.

1. How to check 


for performance?
V4 (Either Int Double) -> (Vector Int, Vector Double)

maximum :: f (Either a b ) -> Either a b

really needs to return a position

(Traversable f, Representable f) 
ifoldMapRep? maybe

bbsearch :: 

allbinary :: (Applicative f, Traversable f) => [f Bool]
allbinary = sequenceA $ pure [True, False]



-}

-- could generalize this to any bounded enum
-- the ord (f b) is ugly as hell
bruteforce :: (Applicative f, Traversable f, Ord a, Ord (f b), Bounded b, Enum b) => (f b -> a) -> (a, f b)
bruteforce f = maximum $ map (\x -> (f x, x)) allbinary where allbinary = sequenceA $ pure  [minBound .. maxBound] -- [True, False]

mosttrue :: (Int, V4 Bool)
mosttrue = bruteforce $ sum . (fmap fromEnum)

-- constrained maximization using a filtering function
bruteforce' :: (Applicative f, Traversable f, Ord a, Ord (f b), Bounded b, Enum b) => (f b -> a) -> (f b -> Bool) ->  (a, f b)
bruteforce' f constraint = maximum $ map (\x -> (f x, x)) $ filter constraint $ allbinary where allbinary = sequenceA $ pure  [minBound .. maxBound]

leasttrue :: (Int, V4 Bool)
leasttrue = bruteforce $ sum . (fmap (negate . fromEnum))

thismanytrue :: Int -> (Int, V4 Bool)
thismanytrue n =  bruteforce $ (\n' -> negate $ abs (n - n')) . sum . (fmap fromEnum)

partialbruteforce :: (Applicative f, Traversable f, Ord a, Ord (f b), Bounded b, Enum b) => (f b -> a) -> f (Maybe b) -> (a, f b)
partialbruteforce obj pf = maximum $ map (\x -> (obj x, x)) $ traverse (maybe [minBound .. maxBound] pure) pf

constrainedmosttrue = partialbruteforce (sum . (fmap fromEnum)) (V4 (Just False) Nothing Nothing (Just True))


{-

pruning search.
a pruning function, tells us whetehr we can possibly beat the current best score
(a -> f (Maybe b) -> Bool)
perhaps the function should also select some good suggestions
a -> f (Maybe b) -> Maybe (f (Maybe b))

a -> f (Maybe b) -> [ f (Maybe b) ] -- return refined possibilities that can possibly beat a.
prunesearch :: (a -> f (Maybe b) -> [ f (Maybe b) ]) -> (f b -> a) -> f (Maybe b) -> (a, f b)
prunesearch prune obj partialf = 



    a relaxation function. A bound. 
    f (Maybe b) -> (a, f (Either b c))
    (f Double) -> (Double, f Double)
    objective :: (f (Either b c) -> a)

using Ordering via objective function.    
instance Galois (f Int) (f Double)
data Galois Int Double = {}

objective ::
f (Either a b) -> c

Left/Right patterns give different domains. The most restrictive is all right, the least restrivie is all left.
Incomparable are incomarable.

relaxedsolve -- does not change LR pattern
prunesolve -- does change LR pattern. 

There is an ordering relations for possible solutions given an LR pattern. The ordering is given by the induced ordering of obj.
we have abstract/concrete pair
concrete :: Pattern -> f (Either a b) -> f (Either a b) 

coherence of type class instances for the LR pattern

instance Galois V2 V2
   concrete = id
   abstract = id
instance  Galois f f',   => V2 f 

-- could carry around the objective function implicitly.
newtype Obj1 f = Obj1 f
instance Ord Obj1 where
    compare x y = compare (obj x) (obj y)

newtype ConeInEq f = 

instance => Lattice (ConeInEq f Double) =
    \/

instance Lattice (V4 a) where
    -- pointwise?


-- not entirely dissimilar from A*
prunesearch bound obj partial

partial solution is path already found. bound is underestimator of path to go. That is A*



    -}