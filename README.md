beliefProp
==========

(Loopy) Belief Propagation on Markov Random Fields with Univariate/Bivariate Potential Functions

Basic usage is described in the code comments.

`beliefProp.m` operates in probability space directly.

`belifPropLog.m` operates in log probability space, using the log-sum-exp trick to prevent potential underflow (http://lingpipe-blog.com/2009/06/25/log-sum-of-exponentials/)

Files can operate in one of 3 modes:

  1. **Marginal Mode**: return the marginal probability of each node in the field.
  2. **Max-marginal Mode**: return the max-marginal probability of each node in the field, for finding best global state.
  3. **Partition Function Mode**: return the normalizing constant of distribution described by the field.

Example code that syllabifies Berber strings is provided:

```
brbrSyll([2 3 7 8 7 1])

ans = 0 1 0 1 0 1
```

For a string with sonority levels [2 3 7 8 7 1], the correct syllable nuclei are in positions [0 1 0 1 0 1].

The syllabification operates using a Markov Random Field that implements the following harmonic grammar:

 * `exp(0)`: unary potential associated with non-nucleus segment
 * `exp(2^s-1)`: unary potential associated with nucleus segment with sonority s (preference for more sonorous nuclei)
 * `exp(0)`: binary potential associated with two adjacent non-nuclei
 * `exp(-2^8)`: binary potential associated with two adjacent nuclei (heavy penalty for violating syllable structure)

