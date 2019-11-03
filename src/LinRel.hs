{-# LANGUAGE ConstraintKinds, ScopedTypeVariables, TypeApplications, AllowAmbiguousTypes
#-}
module LinRel where

import Numeric.LinearAlgebra

type BEnum a = (Enum a, Bounded a) 
enumAll :: (BEnum a) => [a]
enumAll = [minBound .. maxBound]

card :: forall a. (BEnum a) => Int
card = (fromEnum (maxBound @a)) - (fromEnum (minBound @a))

-- LinRel holds A x = b constraint
-- HLinRel
-- VLinRel
data HLinRel a b = HLinRel (Matrix Double) (Vector Double)

-- x = A l + b
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
                                         let (VLinRel q p) = h2v (HLinRel a'' b'') in
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