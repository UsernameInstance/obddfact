#include "fact.hxx"

int main(int argc, char* argv[])
{
	if(argc==2) 
	{
		mpz_class M(argv[1]);
		fact(M);
	}
	else if(argc==3)
	{
		mpz_class M(argv[1]);
		bool showorder = std::stoi(argv[2]);	
		fact(M,showorder);
	}
	else std::cout << "Usage:"
		"\n(1)\tobddfact M"
		"\n(2)\tobddfact M N"
			"\n\t\tM == an integer to factor which may be prefixed with"
				"\n\t\t\t\u201c0x\u201d to interpret as hexadecimal"
				"\n\t\t\t\u201c0b\u201d to interpret as binary"
				"\n\t\t\t\u201c0\u201d  to interpret as octal" 
				"\n\t\t\t\u201c\u201d   to interpret as decimal"
			"\n\t\tN == any integer, where a nonzero value prints variable order to screen after each conjunction"
			<< std::endl;
	return 0;
}
