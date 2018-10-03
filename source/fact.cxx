#include "fact.hxx"
#include "auxiliary.hxx"
#include <fstream>
#include <iomanip>
#include <cassert>
#include <vector> 
#include <ctime> //for std::clock, CLOCKS_PER_SEC
#include <windows.h> //for memory stats
#include <psapi.h> //for memory stats
#include <string>

std::ostream& operator<<(std::ostream &out, const bddStat& stats) 
{
	out << "produced \t==\t" << stats.produced << '\n';
	out << "nodenum \t==\t" << stats.nodenum << '\n';
	out << "max_node_num \t==\t" << stats.maxnodenum << '\n';
	out << "freenodes \t==\t" << stats.freenodes << '\n';
	out << "min_free_nodes \t==\t" << stats.minfreenodes << '\n';
	out << "varnum \t\t==\t" << stats.varnum << '\n';
	out << "cachesize \t==\t" << stats.cachesize << '\n';
	out << "gbcnum \t\t==\t" << stats.gbcnum << '\n';

	return out;
}

void printhandler(std::ostream &out, int variable_index)
{
	char c;
	if(variable_index%2) c = 'y';
	else c = 'x';
	out << c << (variable_index/2);
}

int index_bound(int n, int Lx, int Ly)
{
	assert(Lx >= 1 && Ly >= 1); 
	if(n<=0) return -1;
	if(n<=2) return 0;
	if(n > 1 + Lx + Ly) return n;

	int Lmin = ((Lx<Ly)?Lx:Ly); 
	int Lmax = ((Lx>Ly)?Lx:Ly); 

	if(Lmin+floor_lg(Lmin) >= n-1)
	{
		int m0 = n-2-floor_lg(n-2); //underestimate
		int m1 = n-2-floor_lg(m0); //overestimate
		int m = (1+m0+m1)/2; //now use binary search
		while(m0!=m1)
		{
			if(m + floor_lg(m) <= n-2)
				m0 = m;
			else
				m1 = --m; 
			m = (1+m0+m1)/2;
		}

		return m;
	}
	if (Lmax >= n - ceil_lg(1+Lmin))
	{
		return n-1-ceil_lg(1+Lmin);
	}

	int m0 = Lmax;
	int m1 = Lx+Ly-1;
	int m = (1+m0+m1)/2;
	while(m0 != m1)
	{
		if(ceil_lg(2+Lx+Ly-m)+m <= n-1)
			m0 = m;
		else
			m1 = --m;
		m = (1+m0+m1)/2;
	}
	return m;
}

int order(int n, int m, int N, int Lx, int Ly)
{
	assert(n>=0 && m>=0 && "order n index/rank restriction");
	assert(Lx >= 1 && Ly >= 1);
	if(N<=-2) N=n;
	int u_n = index_bound(n,Lx,Ly);

	if(u_n < N)
	{
		if(m <= N-1-u_n) return 2*(N-m)+1; //return v_ind(y_{N-m}) 
		else if(m <= N+1+u_n)
		{
			if((m-(N-u_n))%2 == 0) return m-(N-u_n); //return v_ind(x_{(m-(N-u_n) div 2})
			return N+u_n+2-m; //return v_ind(y_{1+u_n-((m-(N-1-u_n)) div 2)})
		}
		return 2*(m-N-1); //return v_ind(x_{m-(N+1)})
	}
	
	if(N <= u_n && u_n < 2*N)
	{
		int K = u_n-N;
		if(m<=K) return 2*m;
		if(2*N+1-K <= m) return 2*(2*N+1-m)+1;
		if((m-K)%2==0) return 2*(K + (m-K)/2);
		return 2*(u_n-(K + (m-K)/2))+1;
	}
	
	if(m <= N) return 2*m;
	return 2*(2*N+1-m)+1;
}

bdd build(int n, mpz_class a, int Lx, int Ly)
{
	assert(n>=0 && "bdd index must be nonnegative.");	
	int u_n = index_bound(n,Lx,Ly);
	int Lmin = ((Lx<Ly)?Lx:Ly);
	int Lmax = ((Lx<Ly)?Ly:Lx);
	bool aodd;
	{
		mpz_class temp = (a%2);
		aodd = static_cast<bool>(temp.get_si());
	}
	long int A;
	{
		mpz_class temp = int_div(int_mod(-1-a,std::pow(2,1+n)),std::pow(2,1+u_n));
		A = temp.get_si();
	}
	
	std::vector<std::vector<bdd>> aux(2*(n+1),std::vector<bdd>(int_pow(2,2*(n-u_n))));	//use heap
	//first index is level, second index is integer encoding (Y_m,c_m) <-> c_m+2^{n-u(n)}*Y_m
	//build terminal nodes and nodes with label x_n

	if(n==0) return bdd_ite(bdd_ithvar(1),(A==0)?bdd_ithvar(0):bdd_nithvar(0),(A==0)?bddfalse:bddtrue); //handle n==0 case here.

	for(int l=0; l<int_pow(2,n-u_n+1);++l)
	{
		auto temp = [&](int i, int x)->int{ return (((i%int_pow(2,n-u_n)) + x*int_pow(2,n-1-u_n)*(i/int_pow(2,n-u_n)))%int_pow(2,n-u_n))/int_pow(2,n-1-u_n); };
		int t0 = temp(l,0);
		int t1 = temp(l,1);
		if(n<=Lmin)
		{
			if((t0 == 0) && (t1 == 0)) aux[2*n+1][l] = bddfalse;
			else if((t0 == 0) && (t1 == 1)) aux[2*n+1][l] = bdd_ithvar(2*n);
			else if((t0 == 1) && (t1 == 0)) aux[2*n+1][l] = bdd_nithvar(2*n);
			else aux[2*n+1][l] = bddtrue;
		}
		else
		{
			if(t0 == 0) aux[2*n+1][l] = bddfalse;
			else aux[2*n+1][l] = bddtrue;
		}
	}
	//build levels with labels x_{1+u_n}, ..., x_{n-1}
	for(int m=n-1; 1+u_n<=m; --m) 
	{
		auto c = [&](int i, int x)->int{ return ((i%int_pow(2,n-u_n)) + x*int_pow(2,m-1-u_n)*(i/int_pow(2,n-u_n)))%int_pow(2,n-u_n); };
		int Y;
		for(int l=0; l<int_pow(2,2*n-u_n+1-m); ++l)
		{
			Y = int_pow(2,n-u_n)*((l/int_pow(2,n-u_n))%int_pow(2,n-m));
			aux[n+1+m][l] = ( (m<=Lmin) ? bdd_ite(bdd_ithvar(2*m), aux[n+2+m][c(l,1)+Y], aux[n+2+m][c(l,0)+Y]) : aux[n+2+m][c(l,0)+Y] );
		}
	}	
	//build levels y_{u_n}, x_1, ..., x_{u_n}, y_0
	for(int m=n+1+u_n; n-u_n<m; --m)
	{
		if((n+u_n-m)%2)
		{
			int j = (n+1+u_n-m)/2; //y_j here followed by x_{1+u_n-j}
			for(int l=0; l<int_pow(2,2*(n-u_n)-1); ++l)
			{
				auto temp = [&](int y)->int{ return (l%int_pow(2,n-u_n)) + int_pow(2,n-u_n)*((l/int_pow(2,n-u_n))*2+y); };
				if(j || (!aodd)) aux[m][l] = ( (j<=Lmax) ? bdd_ite(bdd_ithvar(2*j+1), aux[m+1][temp(1)], aux[m+1][temp(0)]) : aux[m+1][temp(0)] );
				else aux[m][l] = aux[m+1][temp(1)];
			}
		}
		else
		{
			int j = u_n - ((n+u_n-m)/2); // x_j here followed by y_{1+u_n-(j+1)}
			for(int l=0; l<int_pow(2,2*(n-u_n)); ++l)
			{
				auto c = [&](int x)->int{ return ( (l%int_pow(2,n-u_n)) + x*(l/int_pow(2,n-u_n)) )%int_pow(2,n-u_n); };
				aux[m][l] = ( (j<=Lmin) ? bdd_ite(bdd_ithvar(2*j), aux[m+1][c(1)+int_pow(2,n-u_n)*( (l/int_pow(2,n-u_n))%int_pow(2,n-u_n-1) )], aux[m+1][c(0)+int_pow(2,n-u_n)*( (l/int_pow(2,n-u_n))%int_pow(2,n-u_n-1) )]) : aux[m+1][c(0)+int_pow(2,n-u_n)*( (l/int_pow(2,n-u_n))%int_pow(2,n-u_n-1) )] );
			}
		}
	}

	//build levels y_{1+u_n}, x_0
	for(int l=0; l<int_pow(2,n-u_n); ++l) //x_0
	{
		auto c = [&](int x)->int{ return ( A + x*l )%int_pow(2,n-u_n); };//c_1
		if(!aodd) aux[n-u_n][l] = bdd_ite(bdd_ithvar(0), aux[n-u_n+1][c(1)+int_pow(2,n-u_n)*( l%int_pow(2,n-u_n-1) )], aux[n-u_n+1][c(0)+int_pow(2,n-u_n)*( l%int_pow(2,n-u_n-1) )]);
		else aux[n-u_n][l] = aux[n-u_n+1][c(1)+int_pow(2,n-u_n)*( l%int_pow(2,n-u_n-1) )];
	}
	for(int l=0; l<int_pow(2,n-u_n-1); ++l) //y_{1+u_n} 
	{
		aux[n-u_n-1][l] = ( (1+u_n <= Lmax) ? bdd_ite(bdd_ithvar(2*(1+u_n)+1), aux[n-u_n][2*l+1], aux[n-u_n][2*l]) : aux[n-u_n][2*l] );
	}

	//build levels y_n, ... y_{2+u_n}
	for(int m=2+u_n; m<=n; ++m)
	{
		for(int l=0; l<int_pow(2,n-m); ++l) aux[n-m][l] = ( (m<=Lmax) ? bdd_ite(bdd_ithvar(2*m+1), aux[n-m+1][2*l+1], aux[n-m+1][2*l]) : aux[n-m+1][2*l] );
	}
	
	return aux[0][0]; 
}


mpz_class fact(mpz_class a, bool show_order)
{
	std::cout << "\na == " << a << "\n";
	std::cout << "a == 0b" << a.get_str(2) << "\n";
	{
		std::string temp = a.get_str(10);
		std::cout << "floor(log(a)) == " << temp.length()-1 << "\n";
	}
	std::cout << "floor(lg(a)) == " << (a.get_str(2)).length()-1 << "\n\n";

	std::clock_t total_time = std::clock();

	std::cout << "Triviality check\n\n";

	if(a < 0) a = -a;
	if(a%2==0) return 2;
	int La = mpz_sizeinbase(a.get_mpz_t(),2)-1;
	int Lx = La/2;
	int Ly = La-1;
	int Lxy = 1+Lx+Ly;

	std::cout << "Initializing BDD package\n";
	int init_num_nodes;
	int max_increase;
	{
		mpz_class temp;
		mpz_root(temp.get_mpz_t(), a.get_mpz_t(), 2); //sets temp to truncated integer part of square root of a.
		++temp;
		double temp_double = mpz_get_d(temp.get_mpz_t());
		temp_double *= 3.6;
		if(temp_double > static_cast<double>(std::numeric_limits<int>::max()))
			init_num_nodes = std::numeric_limits<int>::max();
		else
			init_num_nodes = static_cast<int>(temp_double);
		max_increase = init_num_nodes;
	}
	int init_size_cache = ((init_num_nodes<64)?2:(init_num_nodes/32)); 
	int cache_ratio = ((init_num_nodes<64)?(init_num_nodes/2):32); //apparently cachesize must be at least 2. Therefore cache_ratio <= init_num_nodes/2;
	int max_node_num = 0; //0 interpreted as unlimited/infinity
	int min_free_nodes = 0;  //if percentage free nodes left is less than or equal to min_free_nodes then node table resize is initiated (on next garbage collection).
	bdd_init(init_num_nodes, init_size_cache);
	if(cache_ratio>0) bdd_setcacheratio(cache_ratio);
	bdd_setmaxincrease(max_increase);
	bdd_setmaxnodenum(max_node_num);
	bdd_setminfreenodes(min_free_nodes);
	bdd_setvarnum(2*La); 
	
	std::cout << "bdd_varnum \t==\t" << bdd_varnum() << '\n'; 
	std::cout << "init num nodes \t==\t" << init_num_nodes << "\nmax increase \t==\t" << max_increase << "\ncache ratio \t==\t";
	if(cache_ratio>0) std::cout << cache_ratio;
	else std::cout << "\u221E";
	std::cout << "\nmin free nodes \t==\t" << min_free_nodes << "\n\n";

	bdd_strm_hook(printhandler);
	bdd_gbc_hook(nullptr);

	int* current_order = new int[2*La];
 	auto set_order = [&](int n)->std::clock_t
	{
		std::clock_t order_time = std::clock(); 
		for(int i=0; i<2*La; ++i) 
			current_order[i] = order(n,i,La-1,Lx,Ly);
		bdd_setvarorder(current_order);
		return std::clock()-order_time;
	};
	auto print_current_order = [&]()->void{ 
		int var; std::cout << '['; for(int i=0; i<2*La; ++i) { var = bdd_level2var(i); std::cout << ((var%2)?'y':'x') << (var/2) << ((i<2*La-1)?' ':'\0'); } std::cout << ']'; };

	std::clock_t build_time, order_time(0); 
	size_t size; 
	size_t max_build_time(0), max_order_time(0), max_size(0);
	bddStat stats;
	
	bdd cbdd = bddfalse;
	for(int i = La/2; i > 0; --i)
	{
		bdd temp = bdd_ithvar(2*i);
		for(int k=i+1; k<=La/2; ++k) 
		{
			temp = (temp&bdd_nithvar(2*k));
		}
		
		for(int j=La-i; j >= ((i<=(La-1)/2) ? La-i-1 : La/2); --j)
		{
			bdd temp2 = bdd_ithvar(2*j+1);
			for(int k=j+1; k<=La-1; ++k)
			{
				temp2 = (temp2&bdd_nithvar(2*k+1));
			}
			
			cbdd = (cbdd | (temp & temp2));
		}
	}
	cbdd = cbdd&bdd_ithvar(0)&bdd_ithvar(1);

	std::cout << "Building diagram. Max number loops == " << 2+La << "\n";
	std::cout << std::setw(5) << std::right << "conj" << std::setw(13) << std::right << "order time" << std::setw(13) << std::right << "build time" << std::setw(13) << std::right << "size";
	if(show_order) std::cout << std::setw(13) << std::right << "order";
	std::cout << std::endl;

	for(int i=1+La; i>=0; --i)
	{
		order_time = set_order(i);	
		if(order_time > max_order_time) max_order_time = order_time;
		build_time = std::clock();
		cbdd = cbdd&build(i, a, Lx, Ly); 
		build_time = std::clock()-build_time;
		if(build_time > max_build_time) max_build_time = build_time;
		size = bdd_nodecount(cbdd);
		if(size > max_size) max_size = size;

		std::cout << std::setw(5) << std::right << i << std::setw(13) << std::right << order_time << std::setw(13) << std::right << build_time << std::setw(13) << std::right << size;
	        if(show_order)
		{
			std::cout << std::setw(5) << std::right << ' ';
			print_current_order();
		}	
		std::cout << std::endl;
	}
	std::cout << std::endl;

	bdd solution_bdd = bdd_satone(cbdd);

	std::cout << "solution \t==\t" << solution_bdd << "\n\n";
	total_time = std::clock() - total_time;

	//get solution corresponding to x
	mpz_class solution = 0;
	if(solution_bdd == bddtrue)
	{
		std::cout << "error: solution_bdd == bddtrue" << std::endl;
	}
	else if(solution_bdd == bddfalse)
	{
		solution = 1;
	}
	else
	{
		int index = 0;
		while(true)
		{
			if(solution_bdd == bddtrue)
			{
				break;
			}

			index = bdd_var(solution_bdd);

			if(bdd_low(solution_bdd) != bddfalse)
			{
				solution_bdd = bdd_low(solution_bdd);
				continue;
			}
			else 
			{
				solution_bdd = bdd_high(solution_bdd);
				if(int_mod(index+1,2))
				{
					mpz_class temp = 1;
					int exponent = int_div(index,2);
					for(int i=0; i<exponent; ++i)
					{
						temp *= 2;
					}
					solution += temp;
				}
			}
		}
	}
	
	std::cout << "Closing BDD package\n";

	bdd_stats(stats);
	std::cout << stats << "\n";

	bdd_done();

	delete[] current_order;

	std::cout << "Performance summary\n";
	std::cout << "Max order time \t==\t" << max_order_time << " clicks (" << static_cast<double>(max_order_time)/CLOCKS_PER_SEC << " secs)\n";
	std::cout << "Max build time \t==\t" << max_build_time << " clicks (" << static_cast<double>(max_build_time)/CLOCKS_PER_SEC << " secs)\n";
	std::cout << "Max conj size \t==\t" << max_size << " nodes (" << static_cast<double>(max_size)/50000 << " MB)\n"; 
	std::cout << "Total CPU time \t==\t" << total_time << " clicks (" << static_cast<double>(total_time)/CLOCKS_PER_SEC << " secs)\n";

	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess(), &info, sizeof(info));
	std::cout << "Peak memory \t==\t" << info.PeakWorkingSetSize << " bytes (" << static_cast<double>(info.PeakWorkingSetSize)/std::pow(10,6) << " MB)\n";

	std::cout << "\nRecording performance statistics\n";
	std::ofstream ofs("output.dat", std::ios::app);
	if(!ofs) std::cout << "Error, \u201coutput.dat\u201d could not be opened for writing.\n" << std::endl;
	else
	{
		ofs << a << ' ' << max_order_time << ' ' << max_build_time << ' ' << max_size << ' ' << total_time << ' ' << info.PeakWorkingSetSize << std::endl;
		std::cout << "\n";
	}
	ofs.close();

	return solution;
}

