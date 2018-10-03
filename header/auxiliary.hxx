#ifndef AUXILIARY_HXX
#define AUXILIARY_HXX
#include <gmpxx.h> 
#include <gmp.h>
#include <cmath>
#include <cassert>

//sign
template <typename T> inline int sgn(T val)
{
	return (T(0) < val) - (val < T(0));
}

//absolute value
template <typename T> inline T abs(T val)
{
	return T(sgn(val))*val;
}

//floored integer division
inline long int int_div(long int x, long int b) //floor(x/b)
{
	if(sgn(x)*sgn(b)<0) return -((abs(x)+abs(b)-1)/abs(b));
	return abs(x)/abs(b);
}

inline mpz_class int_div(mpz_class x, mpz_class b) 
{
	mpz_class res;
	mpz_fdiv_q(res.get_mpz_t(), x.get_mpz_t(), b.get_mpz_t());

	return res;
}

//integer modulo
inline mpz_class int_mod(mpz_class x, mpz_class b) //x-b*floor(x/b), where 0*floor(x/0) == 0
{
	if(b==0) return x;
	return x-b*int_div(x,b); 
}


inline long int int_mod(long int x, long int b) 
{
	if(!b) return x;
	return x-b*int_div(x,b);
}

//floor of base 2 logarithm
inline long int floor_lg(mpz_class arg)//for arg positive integer type floor of logarithm base 2 of arg
{
	assert(arg>0); 
	long int ret = 0;
	while((arg = arg/2)) ++ret;
	
	return ret;
}

inline int floor_lg(double arg)
{
	return static_cast<int>(std::log2(arg));
}

//ceil of base 2 logarithm
inline int ceil_lg(double arg)
{
	return static_cast<int>(std::ceil(std::log2(arg)));
}

//integer power function
inline long int int_pow(long int i, long int j)
{
	return static_cast<long int>(std::pow(i,j));
}
#endif
