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