module Control where
import Optimization.Constrained.Penalty
import Data.Functor.Compose
import Data.Functor.Product
import Linear.V

{-

Dynobud had a flavor like this right?
-}

-- control problem with state control vectors lasting for t time steps.
type ControlProb state control t a = (Compose (Product state control) (V t)) a

-- simpleDiscrete :: (forall s. Mode s => f (AD s a) -> AD s a) -> Opt f a -> Opt f a 
--simpleDiscrete :: (forall s. Mode s => f (AD s a) -> AD s a) -> Opt f a -> Opt f a 
--simpleDiscrete xdot = constrainEQ 

-- give lens
