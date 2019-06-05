
{-# LANGUAGE FlexibleContexts #-}
module Limpy where

import Data.Functor.Rep
-- import Cone
import qualified Numeric.Limp.Rep as R
import Numeric.Limp.Program
import Numeric.Limp.Canon
import Numeric.Limp.Solvers.Cbc
import Linear.V3

-- not happy. Let's try glpk


listulate :: (Bounded (Rep f), Enum (Rep f), Representable f, Num a, Eq a) => f a -> [(Rep f, a)]
listulate x =  [ (i, index x i)  | i <- allIndices, index x i /= 0] where allIndices = [minBound .. maxBound]


type HRep f a = [f a]

-- findPoint :: HRep f Double -> f Double
findPoint hrep = solve $ minimise (r1 minBound) constraints [] where
    constraints = (foldr1 (:&&) (map (\f -> let table = [(Right x, R.R n) | (x, n) <- listulate f] in (LR table (0)) :>= conR 0) hrep)) :&& (r1 minBound :== con 1) 

-- Rep V3 is not Ord. Sigh.
ex1 :: HRep V3 Double
ex1 = [V3 0 0 1, V3 1 0 0, V3 1 1 1]

-- unlistulate result = tabulate (\j -> lookup 0 j result)