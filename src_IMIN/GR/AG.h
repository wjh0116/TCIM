#pragma once
#pragma once
#include <chrono>
#include <ctime>
#include <ratio>
#include <algorithm>
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
	

	static void CP_based_AMR(InfGraph& g, const Argument& arg)
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
			inFile.close();
		}
		else
		{
			std::cerr << "Unable to open input file!" << std::endl;
		}
		int rumor_num = g.rumorSet.size();//rumor�ڵ�ĸ���
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
		high_resolution_clock::time_point startTime = high_resolution_clock::now();
		int seed_num = 0;
		double st, et;

		vector<bool>remove_flag;
		remove_flag.resize(g.n, 0);
		for (int node : g.rumorSet)
			remove_flag[node] = 1;
		double memory;
		if (arg.algo == "GR")
		{
			vector<int>tmp_space;
			//out-neighbor
			int k = min((int)CB.size(), arg.k);
			for (int i = 0; i < k; i++)
			{
				tmp_space.clear();
				int seed = g.DecreaseES_ON(10000, g.rumorSet, remove_flag, CB);
				remove_flag[seed] = 1;
				g.seedSet.push_back(seed);
				tmp_space.push_back(seed);
				g.Delete_Node(tmp_space);
			}
			//global
			vector<int>tmp_seedset;
			while (1)
			{
				if (g.seedSet.size() == 0)break;
				int u = g.seedSet.back();
				remove_flag[u] = 0;
				g.seedSet.pop_back();
				g.reset_pro();
				g.Delete_Node(g.seedSet);
				g.Delete_Node(tmp_seedset);
				int seed = g.DecreaseES(10000, g.rumorSet, remove_flag);
				remove_flag[seed] = 1;
				tmp_seedset.push_back(seed);
				if (seed == u)
					break;
			}
			for (int node : tmp_seedset)
			{
				g.seedSet.push_back(node);
			}
			disp_mem_usage(memory);
		}
		else if (arg.algo == "AG")
		{
			vector<int>tmp_space;
			for (int i = 0; i < arg.k; i++)
			{
				tmp_space.clear();
				int seed = g.DecreaseES(10000, g.rumorSet, remove_flag);
				remove_flag[seed] = 1;
				g.seedSet.push_back(seed);
				tmp_space.push_back(seed);
				g.Delete_Node(tmp_space);
			}
		}
		high_resolution_clock::time_point endTime = high_resolution_clock::now();
		duration<double> interval = duration_cast<duration<double>>(endTime - startTime);
		total_time += (double)interval.count();
		g.reset_pro();
		double inf1 = g.MC_based_estimate(g.rumorSet, 50000, arg);
		g.Delete_Node(g.seedSet);
		double inf2 = g.MC_based_estimate(g.rumorSet, 50000, arg);
		double inf3 = inf1 - inf2;
		of << inf1 << "\t" <<  inf2 << "\t" << inf3 << "\t" << (double)interval.count() << "\t"  <<memory<<  endl;
		
		of.close();
	}
};
