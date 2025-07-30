#pragma once
#pragma once
#include <chrono>
#include <ctime>
#include <ratio>
#include <queue>
#include <ctime>
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

class CP
{
private:


public:
	

	static void CP_based(InfGraph& g, Argument& arg)
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
		vector<int>CB;//out_neighbor set
		vector<int>is_ON(g.n, 0);
		for (int node : g.rumorSet)
		{
			for (int i = 0; i < g.gT[node].size(); i++)
			{
				int v = g.gT[node][i];
				if (!is_ON[v])
				{
					CB.push_back(v);
					is_ON[v] = 1;
				}
			}
		}
		ofstream of(arg.res, ios::app);
		g.init_hyper_graph();
		g.deadline_1=arg.T;
		double initial_inf, after_blocking_inf;
		// g.precomputePoisson(3,10000);
		g.precomputeGeometricForAllNodes(10000); // node-specific setting
		g.loadVotingTimesFromFile(arg.dataset + "voting_times.txt");
		// g.loadRealVotingTimesFromFile(arg.dataset + "real_voting_times.txt");
		cout<<"bingo!"<<endl;


		initial_inf = g.MC_based_estimate_inf(g.rumorSet, g.seedSet,  10000, g.deadline_1);
	
		g.get_reachable_node(g.Rnode);
		// cout << "reachable node size: " << g.Rnode.size() << endl;
		
		double initial_ADinf;
		initial_ADinf = g.MC_based_estimate_AD_inf(10000, arg.T) ;


		high_resolution_clock::time_point startTime = high_resolution_clock::now();
		int CPnum = 10000;
		
		
		g.seedSet.clear();

		// double inf = g.estimate_inf_byStop(g.rumorSet, g.seedSet, arg.beta, 1.0 / g.n / 6);

		double inf = initial_inf;

		// cout << "estimate influence by StopAlgorithm is: " << inf << endl;
		double OPT_LB = g.calculate_OPT_lower(arg.k, CB);

		g.opimc_sandwich(arg.k, arg.epsilon, 1.0 / g.n, arg, inf, OPT_LB, 0, 0);
		g.deg_based_heuristic(arg.k, CB);
		g.reset_pro();


		high_resolution_clock::time_point endTime = high_resolution_clock::now();

		double inf_upper = g.estimate_inf_byStop(g.rumorSet, g.UB_seedSet, arg.gamma, 1.0 / g.n);
		cout << "upper's influence:" << inf_upper << endl;

		double inf_lower = g.estimate_inf_byStop(g.rumorSet, g.LB_seedSet, arg.gamma, 1.0 / g.n);
		cout << "lower's influence:" << inf_lower << endl;

		double inf_or = g.estimate_inf_byStop(g.rumorSet, g.Or_seedSet, arg.gamma, 1.0 / g.n);
		cout << "orginal's influence:" << inf_or << endl;
		if (inf_upper <= inf_lower && inf_upper <= inf_or)
			g.seedSet = g.UB_seedSet;
		else if (inf_lower <= inf_upper && inf_lower <= inf_or)
			g.seedSet = g.LB_seedSet;
		else if (inf_or <= inf_lower && inf_or <= inf_upper)
			g.seedSet = g.Or_seedSet;

		g.reset_pro();
		double after_blocking_ADinf;



		after_blocking_ADinf = g.output_AD_influence(arg);
		
		
		duration<double> interval = duration_cast<duration<double>>(endTime - startTime);
		total_time += (double)interval.count();
		cout << "time:" << interval.count() << endl;

		
		// double num_activated_before_vote= g.MC_based_estimate_with_removal_2(10000, arg.T, g.LB_seedSet);
		// cout<<"num_activated_before_vote: "<< num_activated_before_vote<<endl;




		cout << "The AD influence spread of S without blocking is: " << initial_ADinf  << endl;
		cout << "The AD influence spread of S after blocking is: " << after_blocking_ADinf   << endl;

		cout<<"-------------------------"<<endl;
		cout << "The influence spread of S without blocking is: " << initial_inf  << endl;
		after_blocking_inf = g.MC_based_estimate_inf(g.rumorSet, g.seedSet,  10000, g.deadline_1);
		cout << "The influence spread of S after blocking is: " << after_blocking_inf  << endl;




		of <<initial_ADinf<<"\t"<<after_blocking_ADinf<< "\t"<<(double)interval.count()<<  endl;


		of.close();
	}
};
