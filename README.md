# LISA - HyperLogLog based approximation framework for very big Datasets

> **If it looks like a duck, swims like a duck, and quacks like a duck, then it probably is a duck.**

>  *This project is greatly inspired by the implementation of the HLL algorithm in the Julia library:
> https://github.com/jakobnissen/Probably.jl/blob/master/src/hyperloglog/hyperloglog.jl* 

## Introduction

HyperLogLog definition from the Wiki [1]:
> The basis of the HyperLogLog algorithm is the observation that the cardinality of a multiset of uniformly distributed random numbers can be estimated by calculating the maximum number of leading zeros in the binary representation of each number in the set. If the maximum number of leading zeros observed is n, an estimate for the number of distinct elements in the set is 2**n.[2]

> In the HyperLogLog algorithm, a hash function is applied to each element in the original multiset to obtain a multiset of uniformly distributed random numbers with the same cardinality as the original multiset. The cardinality of this randomly distributed set can then be estimated using the algorithm above.

> The simple estimate of cardinality obtained using the algorithm above has the disadvantage of a large variance. In the HyperLogLog algorithm, the variance is minimized by splitting the multiset into numerous subsets, calculating the maximum number of leading zeros in the numbers in each of these subsets, and using a harmonic mean to combine these estimates for each subset into an estimate of the cardinality of the whole set.[3]

## HyperLogLog as an approximation of datasets

The important part of the algorithm described in Introduction is a collection of subsets that we are using to calculate leading zeros.

Lets look at the process of collecting numbers of leading zeros a little bit closer.

In most implementations of HLL ([4], [5]) this collection presented as an array (vector) and named as registers. The size of the registers vector is usually power of 2. For example, Redis implementation [4] is using 2**13 = 8192 elements.

Here is an algorithm to build registries (regvector) using Python as a pseudo-code:

















![Alt text](1_8OoE-wb-4zHDcvWWq2_7xQ.webp)
