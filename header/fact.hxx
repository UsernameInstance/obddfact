#ifndef FACT_HXX
#define FACT_HXX
#include "bdd.h"
#include <iostream>
#include <limits> //for std::numeric_limits<int>::max();
#include <gmpxx.h> 
#include <gmp.h>

//for bddStat object printing
std::ostream& operator<<(std::ostream &out, const bddStat& stats);

//for bdd package printing
void printhandler(std::ostream &out, int variable_index);

//returns largest integer m such that -1<= m <=n and
//
//sum_{k=0}^{m}(2^k*sum_{j=0}^{k}(x_j y_{k-j})) < 2^{n}
//
//where if 0<=j<=Lx then x_j=1 else x_j=0, and
//if 0<=j<=Ly then y_j=1 else y_j=0
int index_bound(int n, int Lx, int Ly);

//returns BuDDy variable index corresponding to mth largest binary variable with respect to the linear order for the nth digit
//y_N < y_{N-1} < ... < y_{1+u} < x_0 < y_u < x_1 < ... < y_1 < x_u < y_0 < x_{1+u} < ... < x_{N-1} < x_N
//where the x_j and y_j are the binary variables and u = index_bound(n, Lx, Ly),
//Lx and Ly represent upper bounds on floor(log2(x)) and floor(log2(y))
//n represents digit index
//N is largest int such that bdd var corresponding to x_N or y_N exists. 2*(1+N) variables. 
int order(int n, int m, int N, int Lx, int Ly);

//builds the index n obdd 
//Assumes bdd_init and bdd_setvarnum have already been called, as well as bdd_setvarorder
//builds bdd for ( (A + \sum_{m=0}^{n}[x_m * \sum_{l=\max(1+u(n),m)}^{n}[2^{l-1-u(n)}*y_{l-m}]]) mod 2^{n-u(n)} ) div 2^{n-1-u(n)}
//using the recurrence c_m =  ( c_{m-1} + x_m * Y_m ) mod 2^{n-u(n)}, c_{-1} = A, where 
//Y_m = \sum_{l=\max(1+u(n),m)}^{n}[2^{l-1-u(n)}*y_{l-m}]]
bdd build(int n, mpz_class a, int Lx = std::numeric_limits<int>::max(), int Ly = std::numeric_limits<int>::max());

//factor the integer a.
mpz_class fact(mpz_class a, bool show_order = false);
#endif
