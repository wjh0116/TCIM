#pragma once
#pragma once
#include <chrono>
#include <ctime>
#include <ratio>
#include <queue>
//#include "infgraph.h"
#define e exp(1)
#define c 2*(exp(1)-2)

using namespace std::chrono;

class Math {
public:
	static double log2(int n) {
		return log(n) / log(2);
	}
	static double logcnk(int n, int k) {
		double ans = 0;
		for (int i = n - k + 1; i <= n; i++)
		{
			ans += log(i);
		}
		for (int i = 1; i <= k; i++)
		{
			ans -= log(i);
		}
		return ans;
	}
};

class Sand
{
private:


public:

	static void SandTCIM(InfGraph& g, Argument& arg)
	{

		sfmt_t sfmtSeed;
		sfmt_init_gen_rand(&sfmtSeed, rand());

		double total_spread = 0;
		double total_time = 0;
		string rumor_file;
		rumor_file = arg.dataset + "rumorSet_" + to_string(arg.Rumor_num) + ".txt";
		ifstream inFile(rumor_file);
		if (inFile.is_open())
		{
			int number;
			while (inFile >> number)
			{
				g.rumorSet.push_back(number);
				g.isRumor[number] = true;
			}
			g.isRumor[g.n] = true;
			inFile.close();
		}
		else
		{
			std::cerr << "Unable to open input file!" << std::endl;
		}
		
		ofstream of(arg.res, ios::app);
		g.precomputePoisson(2,10000);
		g.init_hyper_graph_seed();
		double initial_inf;
		double after_blocking_inf;
		initial_inf = g.MC_based_estimate(10000, arg.T);
		cout << "The TD influence spread of S without blocking is: " << initial_inf  << endl;

		

		high_resolution_clock::time_point startTime = high_resolution_clock::now();	
		
		if(arg.algo=="GB")
			g.GB(arg);
		else
			g.SandTCIM(arg);
		
		high_resolution_clock::time_point endTime = high_resolution_clock::now();
		duration<double> interval = duration_cast<duration<double>>(endTime - startTime);
		total_time += (double)interval.count();
		cout << "time:" << interval.count() << endl;
		double ratio;
		//compute the objective value
		if(arg.algo=="GB")
		{
			after_blocking_inf = g.MC_based_estimate_with_removal(10000, arg.T, g.seedSet);
		}
		else
		{
			after_blocking_inf = g.determine_better_solution(arg);
			ratio=g.rep_approximation(arg);
			cout<<"The approximation ratio is : "<<ratio<<endl;
		}

		
		cout << "The TD influence spread of S after blocking is: " << after_blocking_inf  << endl;

		of <<initial_inf<<"\t"<<after_blocking_inf<< "\t"<<(double)interval.count()<< "\t"<< ratio << endl;
		of.close();
	}
};
