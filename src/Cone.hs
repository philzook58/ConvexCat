{-# LANGUAGE DerivingVia, GeneralizedNewtypeDeriving, DeriveTraversable, ScopedTypeVariables, TypeApplications  #-}

module Cone where

import Linear.Metric
import Linear.Vector
import Linear.V4
import Linear.Epsilon
import Data.Coerce
import Linear.V1
import Data.List (unfoldr)



-- look at Linear.Affine for Point, which is what we're modelling off of. I took some typeclasses out that I don't know what they do.
-- There is a slot for the shape of the space f (_^7 7 dimensional or whatever)  and a slot for the underlying number a.  
newtype Ray f a = Ray (f a) deriving (Eq, Ord, Show, Read, Monad, Functor, Applicative, Foldable, Traversable, Additive, Metric, Fractional , Num, Epsilon)
newtype Dual f a = Dual (f a) deriving (Eq, Ord, Show, Read, Monad, Functor, Applicative, Foldable, Traversable, Additive, Metric, Fractional , Num, Epsilon)
type HalfSpace f a = Dual (Ray f) a -- should I be using Compose?
absorbdual :: Dual (Dual f) a -> f a
absorbdual = coerce

injectdual :: f a -> Dual (Dual f) a
injectdual = coerce

polar :: HalfSpace f a -> Ray f a
polar = coerce

dpolar :: Ray f a -> HalfSpace f a 
dpolar = coerce

-- ConvSet returns Nothing if (f a) is in the set and (Just Dual) giving a hyperplane in which the entire set is contained for which the original argument is not in the set.
-- In other words, it tells you if in set, or gives simple "proof" that you are not in the set. And a clue for possible future inquiries.
type ConvCone f a = f a -> Maybe (Dual f a)
-- Or could call objective
newtype Max a = Max (V1 a) deriving (Eq, Ord, Show, Read, Monad, Functor, Applicative, Foldable, Traversable, Additive, Metric, Fractional , Num, Epsilon)
newtype Domain f a = Domain (f a) deriving (Eq, Ord, Show, Read, Monad, Functor, Applicative, Foldable, Traversable, Additive, Metric, Fractional , Num, Epsilon)

-- newtype ConvexProgSet f a = Max :*: (Domain f) a 
-- newtype
type VRep f a = [Ray f a]
type HRep f a = [HalfSpace f a]

-- an hrep for the convex cone of a single ray. Consists of the polar of the ray to cut off the nagetive side, and a projection of a complete basis. Not linearly independent. Use orthogonalize and prune if you want that.
hrep :: (Metric f, Traversable f, Fractional a) => Ray f a -> HRep f a 
hrep r = pplane : (fmap (projectoff pplane) (basisFor pplane)) where pplane = dpolar r

vrep :: forall f a. (Metric f, Traversable f, Fractional a) => HalfSpace f a -> VRep f a 
vrep = coerce @(Ray f a -> HRep f a) @(HalfSpace f a -> VRep f a) hrep

projectoff :: (Metric v, Fractional a) => v a -> v a -> v a -- subtracts the first argument off of the second 
projectoff b u =  u ^-^ (project b u)

reflectover :: (Metric v, Fractional a) => v a -> v a -> v a -- reflect the second argument about the first
reflectover b u =  u ^-^ (2 *^ (project b u))

elemRH :: (Metric f, Ord a, Num a) => Ray f a -> HalfSpace f a -> Bool
elemRH v h = (polar h) `dot` v >= 0





halfcone :: (Metric f, Ord a, Num a) => HalfSpace f a -> ConvCone (Ray f) a
halfcone h r | r `elemRH` h = Nothing
             | otherwise  = Just h




intersectHH :: HRep f a -> HRep f a -> HRep f a
intersectHH = (++)

hullVV :: VRep f a -> VRep f a -> VRep f a 
hullVV = (++)

-- intersection of VV, hull of HH, and hull of VH are more challenging
-- the intersection of a VRep with an HRep is .. this is not right. We need to project the generators onto the set
--intersectVH :: (Metric f, Num a, Ord a) => VRep f a -> HRep f a -> VRep f a
-- intersectVH vs hs = filter (\v -> all (elemRH v) hs) vs


--findRay :: ConvCone f a -> Ray f a -> [Ray f a] -- this is some relative of iterate. Perhaps unfold
--findRay f r | Nothing <- f r = [r] 
--            | Just h <- f r  =  r : (unfoldr (\r' -> fmap (\h' -> (r', projectoff (dpolar h') r'))) (projectoff (dpolar h) r))
-- The simplest possible findRay. It might be a better idea to reflect over the returned support plane? If you project a polar, it will become zero. not good.
-- reflection is a orthonomral transfromation which is nice. 
-- This is porbably a first order method. Maybe keep a running sum of dual planes with decreasing coefficients 1/k? Line search?
-- The returned dual is a kind of subgradient.
-- don't we really want this to returen [(Ray, Dual Ray)] pairs?
-- what if the convex set is empty? only zero. Then all duals are in it's dual space. But the rays will never converge to zero?
findRay :: (Metric f, Fractional a) => ConvCone (Ray f) a -> Ray f a -> [Ray f a] -- this is some relative of iterate. Perhaps unfold
findRay f r | Nothing <- (f r) = [r] 
            | Just h  <- (f r) = r : (findRay f (projectoff (polar h) r))

orthogonalize :: forall v a. (Metric v, Fractional a) => [v a] -> [v a]
orthogonalize vs = foldl (\bs v -> (projectbasis bs v) : bs) [] vs where 
    projectbasis :: [v a] -> v a -> v a
    projectbasis basis v = foldr projectoff v basis 
prune :: Epsilon a => [a] -> [a]
prune = filter (not . nearZero) 





{-

If we take linear inhomogenous contraints to homgenous linear constraint in d+1, gram Schmidt becomes a solution method for
it also corresponds to working with the augmented matrix.
This is a more geometrical method of finding nullspaces. Gaussian alemination is always performed relative to an extrnal basis. And QR
Sparsity is also not a bisis independent thing

The convexSet representation doesn't lose us that much with resepnt to jus ta plane. We get the plane directly back
The functional representation of a linear map requires reconstitution of the linear map with a solve... what is my point?

linear maps can be thought of as linear relations. I dunno.

a linear subspace is also a convex set.
type LinearSubSpace :: Ray f a -> Maybe (Plane f a)

linearsub :: Plane f a -> 
Solvoing a system of linear equaltites is ocnverting from an Hrep to a VRep


We can also work external solvers into this framework. For example, we might use OSQP.

We could get further clues from the set rather than just a dual plane. We could get our approximate distance for example. The local approximation of the function
d(y) = min { |x - y|  | x \in Set} -> derivative and hessian.
This seems like it would roughly corresopnd to a newton / 2nd order method.
-}