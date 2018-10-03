# OBDDFACT
An implementation of a naive integer factorization algorithm using ordered binary decision diagrams arising from a graduate school project. Primary purpose was to verify that the brute-force approach implemented would result in an exponential (in the logarithm of the number to be factored) time factorization algorithm. For semiprimes data suggests time and space complexities of O(sqrt(a)) where a is the number to be factored.

## What it does
Factors input integer. Gathers performance data (e.g. peak memory use, total CPU time) on use and adjoins to a file labeled "output.dat" in the same directory as the executable.
<p align="center">
  <img width="" height="" src="https://raw.githubusercontent.com/UsernameInstance/obddfact/master/usage.png">
</p>
<p align="center">usage</p>
