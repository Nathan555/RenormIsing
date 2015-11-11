/**
* RenormIsing - Calculates the sums of Boltzman Factors consistent with some
*				coarse graining condition. Six spins in a 1-D periodic 
*				nearest neighbor Ising model are mapped into a 2 spin
*				periodic Ising model with no magnetic field in either case. 
*
*				Exp[kH] => Exp[A(k) + k'H'] 
*
*				Since H' =2 * s'1 * s'2  with all s and s' taking on ±1
*
*				Exp[A(k) + 2k'] = sum of Exp[kH] consistent with s'1  = s2'  
*				Exp[A(k) - 2k'] = sum of Exp[kH] consistent with s'1 != s2'
*
*				The rule to go from {s1...s6} to {s1', s2'} is "majority 
*				rule" such that if the majority of {s1,s2,s3} are +1, then
*				s'1 = +1, and s1' = -1 otherwise. Similiarly {s4,s4,s6} 
*				is mapped to s2'
*
*
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <bitset>
#include <string>
#include <map>
#include <iomanip>

//class for a row of possible configurations of 
//{s1,...s6}
class renorm_row
{
public:
	std::vector<int> spin1;
	int H1;
	std::vector<int> spin2;
	int H2;
};

//applies majority rule for s -> s'
template<typename It>
int majority_rule(It first, It last)
{
	int s = std::accumulate(first, last, 0);
	if (s>0) return 1;
	return -1;
}

//determines the energy of the configuration
//nearest neighbor ising model interaction
//with no magnetic field and periodic boundary conditions
template<typename It>
int energy(It first, It last)
{
	int s = 0;
	for (It k = first; k<last - 1; ++k)
	{
		s += (*k)*(*(k + 1));
	}
	s += (*(last - 1))*(*first);
	return s;
}

int main()
{
	using namespace std;
	vector<renorm_row> v;
	ofstream config_data("cfg_out.txt");//write out raw data
	config_data << "s1 s2 s3\t s1\'\ts4 s5 s6\ts2\'\t\tH\t  H'" << endl;

	//iterate through bitmask configurations
	for (int mask = 0; mask != (1 << 6); ++mask)
	{
		bitset<6> bits(mask);
		v.push_back(renorm_row());
		for (int i = 0; i<6; ++i)
		{
			//store ±1 values for {s1 , ... , s6}
			v.back().spin1.push_back(bits[5-i] ? 1 : -1);
		}
		v.back().H1 = energy(v.back().spin1.begin(), v.back().spin1.end());
		v.back().spin2.push_back(
			majority_rule(v.back().spin1.begin(), v.back().spin1.begin() + 3));
		config_data << setw(3) << v.back().spin1[0] << setw(3) << v.back().spin1[1] << setw(3) << v.back().spin1[2] << "\t";
		config_data << setw(3) << v.back().spin2[0] << '\t';
		v.back().spin2.push_back(
			majority_rule(v.back().spin1.begin() + 3, v.back().spin1.end()));
		v.back().H2 = energy(v.back().spin2.begin(), v.back().spin2.end());
		config_data << setw(3) << v.back().spin1[3] << setw(3) << v.back().spin1[4] << setw(3) << v.back().spin1[5] << "\t";
		config_data << setw(3) << v.back().spin2[1] << "\t\t" << setw(3) << v.back().H1 << '\t'
			<< setw(3) << v.back().H2 <<  endl;

	}
	
	map<int, int> counts_primed_spins_equal,
		counts_primed_spins_unequal;

	//count number of terms where Exp[kH] has the same exponent
	for (auto i = v.begin(); i != v.end(); ++i)
	{
		//s1' = s2' case
		if (i->spin2[0] == i->spin2[1])
		{
			++counts_primed_spins_equal[i->H1];
		}
		//s1' != s2' case
		else
		{
			++counts_primed_spins_unequal[i->H1];
		}
	}
	//output file for equations
	//Exp[A(k) + 2 k'] = sum of Exp[kH] consistent with s'1  = s2'
	//Exp[A(k) - 2 k'] = sum of Exp[kH] consistent with s'1 != s2'
	ofstream ofs("renorm_out.txt");
	ofs << "Exp[A(k)+2k\'] = ";
	bool first_iter = true;
	for (auto i = counts_primed_spins_equal.begin();
	i != counts_primed_spins_equal.end(); ++i)
	{
		if (first_iter)
			first_iter = false; 
		else
			ofs << "+ ";
		ofs <<i->second <<" Exp[" << i->first << " k]";
	}
	ofs << endl;
	ofs << "Exp[A(k)-2k\'] = ";
	first_iter = true;
	for (auto i = counts_primed_spins_unequal.begin();
	i != counts_primed_spins_unequal.end(); ++i)
	{
		if (first_iter)
			first_iter = false;
		else
			ofs << "+ ";
		ofs << i->second << " Exp[" << i->first << " k]";
	}
	ofs << endl;
	return 0;
}