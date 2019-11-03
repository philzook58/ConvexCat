{-# LANGUAGE DerivingVia, GeneralizedNewtypeDeriving, DeriveTraversable, ScopedTypeVariables, TypeApplications  #-}

module Cone where

import Linear.Metric
import Linear.Vector
import Linear.V4
import Linear.Epsilon
import Data.Coerce
import Linear.V1
import Linear.V2
import Linear.V3
import Data.List (unfoldr)
import Data.Functor.Product
-- import Data.Maybe (fromMaybe)



-- look at Linear.Affine for Point, which is what we're modelling off of. I took some typeclasses out that I don't know what they do.
-- There is a slot for the shape of the space f (_^7 7 dimensional or whatever)  and a slot for the underlying number a.  
newtype Ray f a = Ray (f a) deriving (Eq, Ord, Show, Read, Monad, Functor, Applicative, Foldable, Traversable, Additive, Metric, Fractional , Num, Epsilon)
newtype Dual f a = Dual (f a) deriving (Eq, Ord, Show, Read, Monad, Functor, Applicative, Foldable, Traversable, Additive, Metric, Fractional , Num, Epsilon)
type HalfSpace f a = Dual (Ray f) a -- should I be using Compose?

-- I'm cool with these.
absorbdual :: Dual (Dual f) a -> f a
absorbdual = coerce

injectdual :: f a -> Dual (Dual f) a
injectdual = coerce


-- The folllowing feel very ad hoc in their usage. The make me queasy
polar :: HalfSpace f a -> Ray f a
polar = coerce

dpolar :: Ray f a -> HalfSpace f a 
dpolar = coerce

dual :: f a -> (Dual f) a
dual = coerce

dual' :: (Dual f) a -> f a
dual' d = absorbdual (dual d)

-- ConvSet returns Nothing if (f a) is in the set and (Just Dual) giving a hyperplane in which the entire set is contained for which the original argument is not in the set.
-- In other words, it tells you if in set, or gives simple "proof" that you are not in the set. And a clue for possible future inquiries.
type ConvCone f a = f a -> Maybe (Dual f a)
-- forcing it to give you  something extremal in some sense is nice.
type ConvCone' f a = f a -> Maybe (Dual f a, f a) --

{-
Gives nearest cone to given ray and dual supporting plane. Hmm. Maybe I don't have to give f a back because you can just project it onto the given plane if you like.
I think standard projection works fine
We're implcitly using max  (u dot v) / |u| |v| s.t. u in K
This is indeed a reasonable kind of the conification of the support function.
The support function form of a convex set is a black box that can give the maximal point for any linear obective.
I've struggled with

-}

{-
Convex Relation need an intput and output

trans :: ConvRel (f :*: g) h -> ConvRel f (g :*: h)
untrans = inverse

converse :: ConvRel f g -> ConRel g f
compose :: ConvRel f g -> ConvRel g h -> ConvRel f h -- hmm.
meet :: ConvRel f g -> ConvRel f g -> ConvRel f g
join :: ConvRel f g -> ConvRel f g -> ConvRel f g
leftdiv :: --hmm
rightdiv :: --hmm

lift :: ConvCone f -> ConvRel f f 
dup :: ConvRel f (Product f f)
par :: ConRel f g -> ConvRel h k -> ConvRel (Product f h) (Product g k)
projEq  = converse dup 

-}
newtype ConvRel f g a = ConveRel (ConvCone (Product f g) a)

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

projectonto :: (Metric v, Fractional a) => Dual v a -> v a -> v a -- subtracts the first argument off of the second 
projectonto b u =  u ^-^ (project (dual' b) u)

reflectover :: (Metric v, Fractional a) => v a -> v a -> v a -- reflect the second argument about the first
reflectover b u =  u ^-^ (2 *^ (project b u))


score :: (Metric f, Num a) => Ray f a -> HalfSpace f a -> a
score v h = (polar h) `dot` v

elemRH :: (Metric f, Ord a, Num a) => f a -> (Dual f) a -> Bool
elemRH v h = (dual' h) `dot` v >= 0





halfcone :: (Metric f, Ord a, Num a) => HalfSpace f a -> ConvCone (Ray f) a
halfcone h r | r `elemRH` h = Nothing
             | otherwise  = Just h
{-
halfcone' :: (Metric f, Ord a, Num a) => HalfSpace f a -> ConvCone (Ray f) a
halfcone' h r | r `elemRH` h = Nothing
              | otherwise  = Just h
-}
-- From a collection of halfplanes, get the worst scoring one. Seems like a good start for a greedy method. proto simplex.
-- It might be wise to return a collection of the worst socring
-- Or to keep a heap rather than a list? Ehh.
-- I don't really need Ord f a, I need to only sort by the score
hrep' :: (Metric f, Ord a, Num a, Ord (f a)) => HRep f a -> ConvCone (Ray f) a
hrep' hs r = let (hurtiness, h) = minimum (map (\h -> (score r h, h)) hs) in if hurtiness >= 0 then Nothing else Just h



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
-- This is some kind of alternating projection method. With sort of chaotic ordering of the planes.
-- we can also overshoot or undershoot. As I was saying about reflection rather than projection. Supposedly a 1.5 can be pretty good rather than 2.
findRay :: (Metric f, Fractional a) => ConvCone (Ray f) a -> Ray f a -> [Ray f a] -- this is some relative of iterate. Perhaps unfold
findRay f r | Nothing <- (f r) = [r] 
            | Just h  <- (f r) = r : (findRay f (projectoff (polar h) r))

-- orthogonal set of vectors. If we don't ever peek under the cover, we cna guarantee that it is orhtogonalized and pruned.
newtype Orthogonal v a = Orthogonal [v a] -- This is the kind of thing Ghosts of departed proofs might be cool for. We oculd also type level tag orthogonalized according to different metrics
-- really I may want other containers. I could paramtrize over a general traversable perhaps?
-- I'm goingto very likely want a queue. Okasaki queue? Because I'll be popping off the end.
-- And really, I probably want a mutable version of f itself


orthogonalize :: forall v a. (Metric v, Fractional a) => [v a] -> [v a] -- Orthogonal v a
orthogonalize vs = foldl (\bs v -> (projectbasis bs v) : bs) [] vs where 
    projectbasis :: [v a] -> v a -> v a
    projectbasis basis v = foldr project v basis 
prune :: Epsilon a => [a] -> [a]
prune = filter (not . nearZero) 

orthogonalize' :: forall v a. (Epsilon (v a), Metric v, Fractional a) => [v a] -> Orthogonal v a
orthogonalize' hs = foldr appendOrthogonal nilOrthogonal hs

nilOrthogonal :: Orthogonal f a
nilOrthogonal = Orthogonal []
-- identity :: (Num a, Traversable t, Applicative t) => t (t a) 
orthobasis :: (Additive f, Traversable f, Num a) => Orthogonal f a
orthobasis = Orthogonal basis

-- This is where the money happens
appendOrthogonal :: (Epsilon (v a), Metric v, Fractional a) => v a -> Orthogonal v a -> Orthogonal v a
appendOrthogonal h (Orthogonal hs) = let h' = foldr projectoff h hs in if (nearZero h') then (Orthogonal hs) else Orthogonal (h' : hs)

-- This is safe. Might be needed if I don't export the Orthogonal constructor.
forgetOrthogonal :: Orthogonal f a -> [f a]
forgetOrthogonal = coerce
-- headOrthogonal
-- tail Orthogonal
-- dropLastOrhtoognal

-- Does this make sense? I'm not 100% sure.
-- projectOntoPlanes :: (HRep f a, Ray f a) -> (HRep f a, Ray f a) -- return the new orthogonal basis, pruned. We can completely avoid re-orthogonalizing also.
projectOntoPlanes :: (Metric f, Fractional a) => HRep f a -> Ray f a -> Ray f a -- could merge this into a single orthognlaize pass by placing the ray into the 
projectOntoPlanes hs r = let hs' = (orthogonalize hs) in foldr (projectonto) r hs'

projectOntoPlanes' :: (Metric f, Fractional a) => Orthogonal (Dual (Ray f)) a -> Ray f a -> Ray f a -- could merge this into a single orthognlaize pass by placing the ray into the 
projectOntoPlanes' (Orthogonal hs) r = foldr (projectonto) r hs

-- if we greedy ask for hrep', and then projectOntoPlanes 

admmstep :: (Metric f, Fractional a) => ConvCone (Ray f) a -> ConvCone (Ray f) a -> (Ray f a, Ray f a, Ray f a) -> (Ray f a, Ray f a, Ray f a)
admmstep f g (u1, u2, l) = let u2' = proj1 (u1 ^+^ l) in 
                           let u1'  = proj2 (u2' ^-^ l) in
                           let l'  = l ^+^ u1' ^-^ u2' in -- did I get these signs right? I think so.
                           (u1', u2', l') where
    proj1 upl = maybe upl (flip projectonto upl) (f upl) -- it does feel very likely to be some duplication of work here.
    proj2 uml = maybe uml (flip projectonto uml) (g uml)

admm f g = iterate (admmstep f g)

-- no. I'm not doing something right. I;m losing the planes so I don't get closure
-- intersect' :: ConvCone (Ray f) a -> ConvCone (Ray f) a -> ConvCone (Ray f) a
-- intersect' f g r = (take 20 (admm f g)

ex2 = admm f g (Ray $ V2 1 0, Ray $ V2 1 0, Ray $ V2 0 0) where
        f = halfcone (Dual $ Ray (V2 1 0))
        g = halfcone (Dual $ Ray (V2 (-1) 1))



data DD f a = DD {primalDD :: [f a], dualDD :: [Dual f a] }

plane ::  (Functor f, Num a) => HalfSpace f a -> HRep f a
plane h = [h,  (-1) *^ h]
{-
ddray :: (Traversable f, Additive f, Metric f, Fractional a) => Ray f a -> DD (Ray f) a
ddray r = DD [r] hs where hs = concatMap (plane . polar . (projectonto (dpolar r))) basis  -- slightly over complete. We could then orthogonalize to remove also not compiling.

-}
ddhalfspace :: (Traversable f, Additive f, Metric f, Fractional a) => HalfSpace f a -> DD (Ray f) a 
ddhalfspace h = DD rs [h] where rs = polar h : (map (projectonto h) (basisFor (polar h))) -- will be overcomplete. Could orthogonalize. 

ddpolar :: DD f a -> DD (Dual f) a
ddpolar (DD rs hs) = DD (coerce hs) (coerce rs) 



ddintersect :: (Metric f, Ord a, Num a) => DD f a -> DD f a -> DD f a
ddintersect (DD rs hs) (DD rs' hs') = DD  ([r |  r <- rs, all (elemRH r) hs'] ++ [r' |  r' <- rs', all (elemRH r') hs])  (hs ++ hs')

-- a simple pruning. Prune any point that isn't tight on at least one plane. Maybe prune the bigger one first? There can still be tight planes that nevertheless redundant.
-- There is presumably a large literaturee about these problems
-- We may also want to prune any vertices that are literally mutiples of each other. Prune any that are linear sums of others.
-- So any ray should be on at least (D-1?) planes to actually be a corner of the polyhedral cone
-- likewise any plane should touch at least D-1 ray generators.
ddprune (DD rs hs) = let rs' = [r | r <- rs, not (nearZero r), any (\h -> nearZero (score r h)) hs] in
                     let hs' = [h | h <- hs, not (nearZero h), any (\r -> nearZero (score r h)) rs] in
                     DD rs' hs' 


-- cone containement of Hrep. Ax >= 0 ==> Bx >= 0 if exists D >=0 s.t. B >= DA. D is a linear derivation of the second equation from the first  

-- ddhull (DD rs hs) (DD rs' hs') = DD (rs ++ rs') ?
    {-
There is an analong of admm that would work on lists of contraints, if that's what you're into.
if I find a fixed point of admmstep, 
    I kind of get the sense this is not ok. pushing the fixes out later.
    Is there some constraint that I need on the individual lambda?
    admmstep (halfcone  fix admmstep) (halfcone fix admmstep) <=> ???   fix (addmmstep  (admmsteo) (admmstep)) 

Maybe pack together u1 u2 into a larger product space, and then have the full constraint that u11 = u12, u21 = u22. 
Or we could arbitrarily select u1 or u2 from the left or right side 
    Would I want to have some kind of queue system for my convex constraints rather than just sweeping? It seems crazy to keep checking that some irrelevant constraint keeps getting satisfied.
Like count how many times it has not been obtained.


f upl is returning a halfplane, then we use that hyperplane to project. (since that is fast and easy)
could keep a running thing of the post recent hyperplanes we've seen. Dump them if the dimension is greater than f or if they are 0.



I feel like we should be keeping a record of our previous support planes. We're implicitly in the thrall of our previous planes.
And yes, we should keep a queue of the planes, such that the most recently disobeyed ones are at the top? Or the queue could be of fixed size?
Changing the ordering of the queue is a curious manipulation on the lambda.
The lambdas array is about connections between planes, so it should be of size one less that Hrep list.
Maybe as soon as something returns happy, we dump that plane, and combine it's lambdas somehow.
(HalfSpace, Ray, [(HalfSpace, Ray, Lambda)]) Gives inherently correct data structure.
selfadmmstep :: ConvCone -> (HRep, Ray, [Lambdas] -> (HRep, Ray, [Lambdas])
   proj1 (u1 ^+^ lprev ^_^ lnext) -- all except the very first one.
   lprev' = yaday

This makes sense. I might even be convinced that such a method is correct. We're implicitly secviring convex cones as the interesction of there half spaces. We are allow to 
"ask" sort of for particle useful planes, the ones that disprove particular rays are in the set.

This procedure probably does not finitely completely converge? Just gives increasingly good approximations?
Orthogonalization of our halfspaces might help. Then we can gurantee we aren't bouncing into each other's faces. Projection of one does not affect projection of the others.
Maybe a generalization of orthogonlaization. The two directions of a hyper plane are different things.
A complete set of rays under non negative multiples is of size 2*d. Any point can be written as a sum of these.
I dunno. Maybe this makes no sense. Orthogonalization might help constrain you to orthants.

I feel that if your cone has an interior, we might be good. We'll find it? But if it doesn't...
The paper mentions an alpha parameter alpha u  + (1-alpha) u' = u', alpha = 1 is pure. But upper or lower is over relxation and under relaxation
https://web.stanford.edu/~boyd/papers/pdf/scs_long.pdf
Seems very evocative of the symmettric over relaxation. Since ADMM is a guass jacobi style algorithm
 https://en.wikipedia.org/wiki/Successive_over-relaxation
tuning parameters make me ill. I guess guarantees have been out the window anyhow.

Perhaps ConvexSet should be able to return   Feasible | HalfSpace | Plane. Approximating a plane as the intersection of two halfspaces seems like trouble. How do we know we'll ever get both? Well, I guess it both are relevant we'll eventually probe there. No I think this is just an optimization. A good one probably.

 Q: is there a reasonable analog of GMRES for convex problems? GMRES is using a found Krylov space and doing least squares in it.
 We are also doing a found space. I have suppose you could send out our explored planes to an external LP solver.
 One can also prune redundant halfspace constraints using LP solves. Wheels within wheels.

Is what I'm trying to do obscene? Purely functional, no mutation all sorts of wasted doubled ops for the sake of composition.


did i mention yet that we should probably also include series acceleration? a la Hughes perhaps. Anderson, Richardson, who knows
hmm succesevive over relaxation says this is related to rhicardson extarpolation. Interesting.


A different interpetation of the ConvexCone function
What if it is required to return a dual plane such that.
    a. all of the convex cone is in the halfspace
    b. the projection of the given ray onto the halfspace is in the set.
b is a much stronger condition than just the support condition a.
The primitive halfspace function does satisfy this of course.
Fixed admm converges to satisfying this... uh. Maybe
Djiktra projection method. Is dykstra's projection exactly equilvanet to ADMM or not? Has an extra variable, but doesn't seem like it needs one? It does all told look more symmetrical


I guess my concern is that the thing could return an exact 

If in the interior of the cone, we could have it return nothing or perhaps a similar hyperplane?

https://web.stanford.edu/class/ee392o/alt_proj.pdf

-}
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

we need to start building a tower of embeddings.
self dual embeddings,
affine embeddings,
function to set problems
minimization to feasibility problems

-}

{-


Does the optimization package make all this ridiculous?

I like the penalty interior point method. Hmm. Actually, wait, they are using penalty, not interior point. Why is that?
Anyhow, that is good enough probably for a demo mixed integer program using logicT / list search.



If it's easy, maybe it makes sense to compile to glpk or cbc
https://kul-forbes.github.io/scs/
https://www.chrisstucchio.com/blog/2012/linear_programming.html


A couple algorithms:
Just iterated projection
ADMM styel. sort of damped projection
take SVD of all active planes. Kind of an interior point method
Do sinlge round of QR algorithm for svd. Low rank update of 
    Guass newton aplicaable?
    min_{|x|=1} (a^T x) ^2 / (a^T a)

SVD

-}

{-
# Convex Programming in Haskell Stuff

I feel pulled in multiple directions. I want to wait to write a blog post until the ideas in it have some level of completion. But sometimes that level never comes. My attention may wander off. And yet, I feel I often have some interesting half built things to show.
I have many loves in my life, but two of them are Haskell and Convex Programming.

In a previous post, I talked a bit about how reducing convex programming to questions about cones leads to elegance and simplification.

There are a couple of approaches to well-typed vectors in Haskell. One of them is to consder an f :: * -> * parameter as characterizing a vector space. We intend to fill in this with things like V4.
data V4 a = V4 a a a a 
A vector space can be consdiered $R^n$. The f shape is sort of the $-^n$ part. The part that tells us the shape/ size of the vector space. We can then choose to fill it in, making is a vector space over Doubles or Complex Doubles or something more exotic.
The fundamental construct we want to talk about for cones is the Ray, rather than the vector or point. A Ray is a direction. We can represent it with a vector, if we ignore the magnitude of the vector or if we choose to always work with normalized vectors (which is less elegant really).
A cone is a set that is closed under addition and non negative scalar multiplication. It is a convex set of rays.
The dual space to rays are the halfspace cones. They can be naturally parametrized also by vectors (basically the vectors normal to the planes pointing in the direction of halfspace you want). Any ray that has positive dot product with this vector is in the halfspace, which gives a simple test.

Polytopes have (at least) two natural representations. They can be expressed as a sum of generator rays (the corners of the polytope) or as the set obeying a set of halfplane constraints. The first is called the V representation and the latter the H representation. The two are highly interconnected by duality.

The difficulty and method for solving these qeustions can depend strongly on the representation you have

What are geoemtrically interesting questions to ask. 
1. Whether a ray is in a cone.
   + Easy given HRep. Just do all the dot products
   + 
2. Finding a ray that is in a cone
  - Any ray or an extremal ray? Extremal according to what?
  - Projecting to a cone. There are metrics available. The cosine metric. The Chord metric. Others? They feel a bit goofy. They feel far less elegant than the quadratic metric of a vector space. The dot product is a somewhat natural thing to talk about, but by design. Abs metric? You can use the unit plane. and put a metric on it. The cone space is a relative of projective space. which has similar problems. The cross ratio? Logarithm of dot product? Maybe you do need a plane to work on to have a reasonable deifition of extremal. Etremal with respect to a cone generalized inequality on the euclidean space. Extremal with respect to another cone? But a particular kind of cone. It contains (0,0,0,1)? A projection excluding (0,0,0,1)


4. Convert VRep to HRep and vice versa
   + take every possible choice out of set. Do linear solve. Check if linear solve matches everything else.
5. Intersections
  + Alternating Projecction. The most primitive and obvious.  
  + ADMM. 
6. Minkowski Sums = Convex Hull
7. Projection
8. Pruning Redundancy
9. Cone Containment Testing.

-}