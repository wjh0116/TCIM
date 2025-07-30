#include "iheap.h"
#include <queue>	//priority_queue
#include <utility>  // pair
#include <numeric>
#include <algorithm>
#include <random>
class InfGraph: public Graph
{
private:
	//vector<bool> activated;
    vector<bool> visit;
    vector<int> visit_mark;
	
public:
    vector<vector<int>> hyperG;
    vector<vector<int>> hyperGT;
	vector<vector<int>> SampleG;
	vector<vector<int>> SampleGT;
	vector<vector<int>> hyperG_2;
	vector<vector<int>> hyperGT_2;
	vector<vector<int>> SampleG_2;
	vector<vector<int>> SampleGT_2;
	vector<int>Rnode;
	vector<vector<int>> LhyperG;//LRR1
	vector<vector<int>> LhyperGT;
	vector<vector<int>> LhyperG_2;//LRR2
	vector<vector<int>> LhyperGT_2;//LRR2
	vector<double>influenced;//for estimate opt
	vector<double>influenced_deg;//for heuristic
	vector<int> ActNode;
	int tmp_count;
	int S;
	double sumc;
	vector<int> Dfn; 
	vector<int> Ord; 
	vector<int> Parent; 
	vector<double>C;
	// vector<int> precomputedPoisson;
	int Stamp;
	vector<int> uni; 
	vector<int> mn;  
	vector<int> Sdom; 
	vector<int> Idom; 
	vector<vector<int>> SdomTree;

	sfmt_t sfmtSeed;
	vector<int> UB_seedSet;
	vector<int> LB_seedSet;
	vector<int> Or_seedSet;
	vector<int> seedSet;
	vector<int> rumorSet;
	vector<int> ReachableNode;
	vector<bool> isReach;
	vector<bool> isRumor;
	vector<bool> isSelect;
	vector<int> Dec;
	int deadline_1;
	vector<vector<int>> precomputedGeometric; // 为每个节点预计算的几何分布采样
	vector<double> stratifiedVotingTimes; // 分层采样的投票时间
	vector<double> realVotingTimes; // 存储真实投票时间

    InfGraph(string folder, string graph_file): Graph(folder, graph_file)
    {
        srand(time(NULL));
        sfmt_init_gen_rand(&sfmtSeed , rand());		
        visit = vector<bool> (n+1);
        visit_mark = vector<int> (n+1);
		Dec = vector<int>(n+1);
		isRumor = vector<bool>(n+1);
		C = vector<double>(n + 1, 0);
		influenced = vector<double>(n + 1, 0);
		influenced_deg = vector<double>(n + 1, 0);
		isSelect = vector<bool>(n+1);
		isReach = vector<bool>(n + 1);
		//activated = vector<bool>(n, false);
		hyperG.resize(n+1, vector<int>());
		LhyperG.resize(n + 1, vector<int>());
		hyperG_2.resize(n + 1, vector<int>());
		LhyperG_2.resize(n + 1, vector<int>());
		SampleG.resize(n+1, vector<int>());
		SampleG_2.resize(n + 1, vector<int>());
		precomputedGeometric = vector<vector<int>>(n);
		stratifiedVotingTimes = vector<double>(n);
		realVotingTimes = vector<double>(n); // 存储真实投票时间
    }

    void init_hyper_graph(){
		for (auto& hyper : hyperG)hyper.clear();
		for (auto& hyperT : hyperGT)vector<int>().swap(hyperT);
		hyperGT.clear();
		for (auto& hyper : LhyperG)hyper.clear();
		for (auto& hyperT : LhyperGT)vector<int>().swap(hyperT);
		LhyperGT.clear();
		for (auto& hyper : SampleG)hyper.clear();
		for (auto& hyperT : SampleGT)vector<int>().swap(hyperT);
		SampleGT.clear();

		for (auto& hyper : hyperG_2)hyper.clear();
		for (auto& hyperT : hyperGT_2)vector<int>().swap(hyperT);
		hyperGT_2.clear();
		for (auto& hyper : LhyperG_2)hyper.clear();
		for (auto& hyperT : LhyperGT_2)vector<int>().swap(hyperT);
		LhyperGT_2.clear();
		for (auto& hyper : SampleG_2)hyper.clear();
		for (auto& hyperT : SampleGT_2)vector<int>().swap(hyperT);
		SampleGT_2.clear();
		//seedSet.clear();
		for (int i = 0; i < n; i++)
			isSelect[i] = false;
    }

	char* map_file(const char* fname, size_t& length)
	{
		int fd = open(fname, O_RDONLY);
		if (fd == -1)
			handle_error("open");

		// obtain file size
		struct stat sb;
		if (fstat(fd, &sb) == -1)
			handle_error("fstat");

		length = sb.st_size;

		char* addr = static_cast<char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
		if (addr == MAP_FAILED)
			handle_error("mmap");

		// TODO close fd at some point in time, call munmap(...)
		close(fd);
		return addr;
	}

	/// Mente-Carlo simulations to estimate the influence, exclusing seeds themselves
	double MC_based_estimate_inf(vector<int> rumorNode, vector<int> Blockers,  int count, int deadline)
	{
		vector<int> isBlocker(n, 0);
		for(int node : Blockers) {
			isBlocker[node] = 1;
		}
		vector<int> vis;
		vis.resize(n);
		double sum = 0;
		for (int tag = 1; tag <= count; tag++)
		{
			priority_queue<Node2> Q;
			vector<int> AT(n, 10000); 
			vector<int> I(n, 0); 
			for (int u : rumorSet) {
				AT[u] = 0;
				Q.push({u, 0});
			}
			while (!Q.empty()) 
			{
				Node2 current = Q.top();
				Q.pop();
				int q = current.id;
				int q_AT = current.activate_time;
				if (I[q] == 1) continue;

				for(int i=0;i<gT[q].size();i++)
				{
					if(isBlocker[gT[q][i]]) continue; 
					int tmp= sfmt_genrand_uint32(&sfmtSeed) % precomputedGeometric[q].size();
					int diff_delay=precomputedGeometric[q][tmp]; 
					int u=gT[q][i];
					if(sfmt_genrand_real1(&sfmtSeed) <= probT[q][i] && AT[q] + diff_delay <= deadline)
					{
						if (AT[q] + diff_delay < AT[u]) {
							AT[u] = AT[q] + diff_delay;
							Q.push({u, AT[u]});
						} 
					}
				}	
				I[q] = 1;
				sum++;
			}
		}
		return sum / ((double)count)- rumorNode.size();
	}

	/// Mente-Carlo simulations to estimate the misinformation's TD influence
	double MC_based_estimate_AD_inf(int count, int deadline)
	{
		int cnt=0;
		double sum=0;
		while (cnt<count)
		{
			priority_queue<Node2> Q;
			vector<int> AT(n, 10000); 
			vector<int> I(n, 0); 
			for (int u : rumorSet) {
				AT[u] = 0;
				Q.push({u, 0});
			}
			while (!Q.empty()) 
			{
				Node2 current = Q.top();
				Q.pop();
				int q = current.id;
				int q_AT = current.activate_time;
				if (I[q] == 1) continue;

				for(int i=0;i<gT[q].size();i++)
				{
					int tmp= sfmt_genrand_uint32(&sfmtSeed) % precomputedGeometric[q].size();
					int diff_delay=precomputedGeometric[q][tmp]; 
					int u=gT[q][i];
					if(sfmt_genrand_real1(&sfmtSeed) <= probT[q][i] && AT[q] + diff_delay <= deadline)
					{
						if (AT[q] + diff_delay < AT[u]) {
							AT[u] = AT[q] + diff_delay;
							Q.push({u, AT[u]});
						} 
					}
				}	
				I[q] = 1;
				// sum+=U(deadline-AT[q]);
				sum+=computeU(q, AT[q]);
			}
			cnt++;
		}
		return sum/(double)count;
	}

	// void precomputePoisson(double lambda, int sampleSize) {
	// 	static thread_local std::mt19937 gen(std::random_device{}());
	// 	std::poisson_distribution<> poisson(lambda);

	// 	precomputedPoisson.reserve(sampleSize);
	// 	for (int i = 0; i < sampleSize; ++i) {
	// 		precomputedPoisson.push_back(poisson(gen) + 1);
	// 	}
	// }


	void Delete_Node(vector<int> Nodeset)
	{
		for (int node : Nodeset)//forward
		{
			for (int i = 0; i < (int)gT_reverse[node].size(); i++)//node��ȫ����ߵĸ�������Ϊ0
			{
				int v = gT_reverse[node][i];
				for (int j = 0; j < (int)gT[v].size(); j++)
				{
					if (gT[v][j] == node)
					{
						probT[v][j] = 0;
						break;
					}
				}
			}
		}
	}
	
	void reset_pro()
	{
		probT.assign(probT2.begin(),probT2.end());
	}

	void gen_sample_2(vector<vector<int>>& sample, vector<int>rumorset, vector<int>& ReachableNode, vector<int>& vis, int tag) {

		queue<int> q;

		for (auto x : rumorset)
			q.push(x), vis[x] = -1;
		while (!q.empty())
		{
			int x = q.front();
			q.pop();
			for (int i = 0; i < gT[x].size(); i++)
			{
				if (probT[x][i] == 0)
					continue;
				if (probT[x][i] >= sfmt_genrand_real1(&sfmtSeed))
				{
					sample[x].push_back(gT[x][i]);
					if (vis[gT[x][i]] != -1)
						q.push(gT[x][i]), vis[gT[x][i]] = -1, ReachableNode.push_back(gT[x][i]);
				}
			}
		}
	}

	double select_seed_low(const Argument& arg)
	{
		vector<int> coverage(n, 0); 
		size_t maxDeg = 0; 
		for (auto i = n; i--;) { 
			const auto deg = hyperG[i].size();
			coverage[i] = deg; 
			if (deg > maxDeg) maxDeg = deg; 
		}
		vector<vector<int>> degMap(maxDeg + 1); 
		for (auto i = n; i--;) { 
			degMap[coverage[i]].push_back(i); 
		}
		vector<int> sortedNode(n); 
		vector<int> nodePosition(n); 
		vector<int> degreePosition(maxDeg + 2); 
		uint32_t idxSort = 0;
		size_t idxDegree = 0;
		for (auto& nodes : degMap) { 
			degreePosition[idxDegree + 1] = degreePosition[idxDegree] + (uint32_t)nodes.size();
			idxDegree++;
			for (auto& node : nodes) {
				nodePosition[node] = idxSort;
				sortedNode[idxSort++] = node;
			}
		}
		std::vector<bool> edgeMark(hyperGT.size(), false);
		// record the total of top-k marginal gains
		size_t sumTopk = 0;
		for (auto deg = maxDeg + 1; deg--;)
		{
			if (degreePosition[deg] <= n - arg.k)
			{
				sumTopk += deg * (degreePosition[deg + 1] - (n - arg.k));
				break;
			}
			sumTopk += deg * (degreePosition[deg + 1] - degreePosition[deg]);
		}
		double __boundMin = 1.0 * sumTopk;
		size_t sumInf = 0;
		LB_seedSet.clear();
		for (auto k = arg.k; k--;) {
			const auto seed = sortedNode.back(); 
			sortedNode.pop_back(); 
			int newNumV = sortedNode.size();
			sumTopk += coverage[sortedNode[newNumV - arg.k]] - coverage[seed];
			sumInf += coverage[seed];
			LB_seedSet.push_back(seed); 
			if (k == 0)
			{
				return __boundMin;
			}
			coverage[seed] = 0; 
			for (auto edgeIdx : hyperG[seed]) { 
				if (hyperGT[edgeIdx].size() == 0 || edgeMark[edgeIdx]) continue; 
				edgeMark[edgeIdx] = true; 
				for (auto nodeIdx : hyperGT[edgeIdx]) { 
					if (coverage[nodeIdx] == 0) continue; 
					const auto currPos = nodePosition[nodeIdx]; 
					const auto currDeg = coverage[nodeIdx]; 
					const auto startPos = degreePosition[currDeg]; 
					const auto startNode = sortedNode[startPos]; 
					// Swap this node to the start position with the same degree, and update their positions in nodePosition
					std::swap(sortedNode[currPos], sortedNode[startPos]); 
					nodePosition[nodeIdx] = startPos; 
					nodePosition[startNode] = currPos;
					// Increase the start position of this degree by 1, and decrease the degree of this node by 1
					degreePosition[currDeg]++; 
					coverage[nodeIdx]--; 
					// If the start position of this degree is in top-k, reduce topk by 1
					if (startPos >= newNumV - arg.k) sumTopk--;
				}
			}
			double __boundLast = 1.0 * (sumInf + sumTopk);
			if (__boundMin > __boundLast) __boundMin = __boundLast;
			
		}
		return __boundMin;
	}

	double select_seed_upper(const Argument& arg)
	{
		vector<int> coverage(n, 0); 
		size_t maxDeg = 0; 
		for (auto i = n; i--;) { 
			const auto deg = LhyperG[i].size();
			coverage[i] = deg; 
			if (deg > maxDeg) maxDeg = deg; 
		}
		vector<vector<int>> degMap(maxDeg + 1); 
		for (auto i = n; i--;) { 
			//if (coverage[i] == 0) continue;
			degMap[coverage[i]].push_back(i); 
		}
		vector<int> sortedNode(n); 
		vector<int> nodePosition(n); 
		vector<int> degreePosition(maxDeg + 2); 
		uint32_t idxSort = 0;
		size_t idxDegree = 0;
		for (auto& nodes : degMap) { 
			degreePosition[idxDegree + 1] = degreePosition[idxDegree] + (uint32_t)nodes.size();
			idxDegree++;
			for (auto& node : nodes) {
				nodePosition[node] = idxSort;
				sortedNode[idxSort++] = node;
			}
		}
		std::vector<bool> edgeMark(LhyperGT.size(), false);
		// record the total of top-k marginal gains
		size_t sumTopk = 0;
		for (auto deg = maxDeg + 1; deg--;)
		{
			if (degreePosition[deg] <= n - arg.k)
			{
				sumTopk += deg * (degreePosition[deg + 1] - (n - arg.k));
				break;
			}
			sumTopk += deg * (degreePosition[deg + 1] - degreePosition[deg]);
		}
		double __boundMin = 1.0 * sumTopk;
		size_t sumInf = 0;
		UB_seedSet.clear();
		for (auto k = arg.k; k--;) {
			const auto seed = sortedNode.back(); 
			sortedNode.pop_back(); 
			int newNumV = sortedNode.size();
			sumTopk += coverage[sortedNode[newNumV - arg.k]] - coverage[seed];
			sumInf += coverage[seed];
			UB_seedSet.push_back(seed); 
			//cout << "seed: " << seed << "cover: " << coverage[seed] << endl;
			if (k == 0)
			{
				//cout << "k have equalled to 0!" << endl;
				return __boundMin;
			}
			coverage[seed] = 0; 
			for (auto edgeIdx : LhyperG[seed]) { 
				if (LhyperGT[edgeIdx].size() == 0 || edgeMark[edgeIdx]) continue; 
				edgeMark[edgeIdx] = true; 
				for (auto nodeIdx : LhyperGT[edgeIdx]) { 
					if (coverage[nodeIdx] == 0) continue; 
					const auto currPos = nodePosition[nodeIdx]; 
					const auto currDeg = coverage[nodeIdx]; 
					const auto startPos = degreePosition[currDeg]; 
					const auto startNode = sortedNode[startPos]; 
					// Swap this node to the start position with the same degree, and update their positions in nodePosition
					std::swap(sortedNode[currPos], sortedNode[startPos]); 
					nodePosition[nodeIdx] = startPos; 
					nodePosition[startNode] = currPos;
					// Increase the start position of this degree by 1, and decrease the degree of this node by 1
					degreePosition[currDeg]++; 
					coverage[nodeIdx]--;
					// If the start position of this degree is in top-k, reduce topk by 1
					if (startPos >= newNumV - arg.k) sumTopk--;
				}
			}
			double __boundLast = 1.0 * (sumInf + sumTopk);
			if (__boundMin > __boundLast) __boundMin = __boundLast;

		}
		return __boundMin;
		//return sumInf;
	}


	struct tl
	{
		int clk;
		vector<int> fa, idx, ridx, c, best, semi, idom, sum;
		vector<vector<int>> dt, G, rG;
		void init(int n, vector<vector<int>> e)
		{
			G.clear();
			G = e;
			fa.clear();
			fa.resize(n + 1);
			idx.clear();
			idx.resize(n + 1);
			ridx.clear();
			ridx.resize(n + 1);
			c.clear();
			c.resize(n + 1);
			best.clear();
			best.resize(n + 1);
			dt.clear();
			dt.resize(n + 1);
			rG.clear();
			rG.resize(n + 1);
			semi.clear();
			semi.resize(n + 1);
			idom.clear();
			idom.resize(n + 1);
			sum.clear();
			sum.resize(n + 1);
			clk = 0;
			for (int i = 0; i <= n; i++)
			{
				sum[i] = 0;
				c[i] = -1;
				semi[i] = i;
				best[i] = i;
				idx[i] = 0;
				for (int v : G[i])
					rG[v].push_back(i);
			}
		}
		void dfs(int u)
		{
			idx[u] = ++clk;
			ridx[clk] = u;
			for (int& v : G[u])
				if (!idx[v])
				{
					fa[v] = u;
					dfs(v);
				}
		}
		int fix(int x)
		{
			if (c[x] == -1)
				return x;
			int& f = c[x], rt = fix(f);
			if (idx[semi[best[x]]] > idx[semi[best[f]]])
				best[x] = best[f];
			return f = rt;
		}

		void compute(int x)
		{
			sum[x] = 1;
			for (auto y : dt[x])
			{
				compute(y);
				sum[x] += sum[y];
			}
		}
		void go(int rt)
		{
			dfs(rt);
			for (int i = clk; i > 1; i--)
			{
				int x = ridx[i], mn = clk + 1;
				for (int& u : rG[x])
				{
					if (!idx[u])
						continue;
					fix(u);
					mn = min(mn, idx[semi[best[u]]]);
				}
				c[x] = fa[x];
				dt[semi[x] = ridx[mn]].push_back(x);
				x = ridx[i - 1];
				for (int& u : dt[x])
				{
					fix(u);
					if (semi[best[u]] != x)
						idom[u] = best[u];
					else
						idom[u] = x;
				}
				dt[x].clear();
			}
			for (int i = 2; i <= clk; i++)
			{
				int u = ridx[i];
				if (idom[u] != semi[u])
					idom[u] = idom[idom[u]];
				dt[idom[u]].push_back(u);
				// cout<<idom[u]<<" "<<u<<endl;
			}
			//compute(rt);
		}
		void construct_CPset(int rt, vector<vector<int>>& hyperG, vector<vector<int>>& hyperGT, vector<bool>isRumor, int& j)
		{
			vector<int>tg;
			tg.resize(rt);//tg[u]=j,u denotes node, j denotes the tag of CPset
			for (int i = 2; i <= clk; i++)
			{
				int v = ridx[i];
				if (isRumor[v] == 0)
				{
					if (isRumor[idom[v]] == 1)
					{
						hyperG[v].push_back(j);
						hyperGT[j].push_back(v);
						tg[v] = j;
						j++;
					}
					else if (isRumor[idom[v]] == 0)
					{
						hyperG[v].push_back(j);
						hyperGT[j].push_back(v);
						int tag = tg[idom[v]];
						for (int node : hyperGT[tag])
						{
							hyperG[node].push_back(j);
							hyperGT[j].push_back(node);
						}
						tg[v] = j;
						j++;
					}
				}

			}
		}


		void construct_LRRset_v2(int n, int node, vector<vector<int>>& LhyperG, vector<vector<int>>& LhyperGT, int tag, vector<bool>isRumor, vector<vector<int>> sample_reverse)
		{
			auto uStart = node;
			unsigned int n_visit_mark = 0, curIdx = 0;
			vector<int>tmp_visit_mark;
			tmp_visit_mark.resize(n, 0);
			vector<bool>tmp_visit;
			tmp_visit.resize(n, 0);
			tmp_visit_mark[n_visit_mark++] = uStart;
			tmp_visit[uStart] = true;
			LhyperG[uStart].push_back(tag);
			while (curIdx < n_visit_mark) {
				int i = tmp_visit_mark[curIdx++];
				for (int j = 0; j < (int)sample_reverse[i].size(); j++) {
					int v = sample_reverse[i][j];
					if (tmp_visit[v] || isRumor[v] )continue;
					tmp_visit[v] = true;
					tmp_visit_mark[n_visit_mark++] = v;
					LhyperG[v].push_back(tag);
				}
			}
			LhyperGT.push_back(vector<int>(tmp_visit_mark.begin(), tmp_visit_mark.begin() + n_visit_mark));
		}
	} tree;

	int generateCP_local_global(int count, const Argument& arg, int presize, int pre_subsize, int Uflag, int Lflag)
	{
		int j = pre_subsize;//number of CP sets
		int num = presize;//number of CP sequences
		int LRRnum = LhyperGT.size();
		while (num < count)
		{
			ReachableNode.clear();
			vector<vector<int>> sample;
			sample.resize(n);
			vector<vector<int>> sample_reverse;
			sample_reverse.resize(n);
			vector<int> vis;
			vis.resize(n, 0);
			gen_sample_2(sample, rumorSet, ReachableNode, vis, num);
			if (Uflag == 0)
			{
				int ran = sfmt_genrand_uint32(&sfmtSeed) % Rnode.size();
				int Idx = Rnode[ran];
				if (vis[Idx] == -1)
				{

					for (int i = 0; i < n; i++)
					{
						for (int v : sample[i])
							sample_reverse[v].push_back(i);
					}
					tree.construct_LRRset_v2(n, Idx, LhyperG, LhyperGT, LRRnum, isRumor, sample_reverse);
					LRRnum++;
				}
			}
			if (Lflag == 0)
			{
				sample.push_back(rumorSet);
				tree.init(n, sample);
				tree.go(n);


				for (int i = 0; i < ReachableNode.size(); i++)
				{
					hyperGT.push_back(vector<int>());
				}
				tree.construct_CPset(n, hyperG, hyperGT, isRumor, j);
			}
			num++;
		}
		return j;
	}
	int generateCP_local_global_2(int count, const Argument& arg, int presize, int pre_subsize, int Uflag, int Lflag)
	{
		int j = pre_subsize;
		int num = presize;
		int LRRnum = LhyperGT_2.size();
		while (num < count)
		{
			ReachableNode.clear();
			vector<int>temp(n);
			int tmp = 0;
			vector<vector<int>> sample;
			sample.resize(n);
			vector<vector<int>> sample_reverse;
			sample_reverse.resize(n);
			vector<int> vis;
			vis.resize(n, 0);
			gen_sample_2(sample, rumorSet, ReachableNode, vis, num);
			if (Uflag == 0)
			{
				int ran = sfmt_genrand_uint32(&sfmtSeed) % Rnode.size();
				int Idx = Rnode[ran];
				if (vis[Idx] == -1)
				{
					for (int i = 0; i < n; i++)
					{
						for (int v : sample[i])
							sample_reverse[v].push_back(i);
					}
					tree.construct_LRRset_v2(n, Idx, LhyperG_2, LhyperGT_2, LRRnum, isRumor, sample_reverse);
					LRRnum++;
				}
			}
			if (Lflag == 0)//generate CP2
			{
				sample.push_back(rumorSet);
				tree.init(n, sample);
				tree.go(n);

				for (int i = 0; i < ReachableNode.size(); i++)
				{
					hyperGT_2.push_back(vector<int>());
				}
				tree.construct_CPset(n, hyperG_2, hyperGT_2, isRumor, j);
			}
			num++;
		}

		return j;

	}

	static bool cmp(const pair<int, double>& a, const pair<int, double>& b) {
		return a.second > b.second;
	}

	void get_reachable_node(vector<int>& ReachableNode)
	{
		vector<int> vis;
		vis.resize(n ,0);
		
		
		queue<int> q;
		for (int x : rumorSet)
			q.push(x), vis[x] = 1;
		while (!q.empty())
		{
			int x = q.front();
			q.pop();
			for (int i = 0; i < gT[x].size(); i++)
			{
				if (vis[gT[x][i]] != 1)
					q.push(gT[x][i]), vis[gT[x][i]] = 1, ReachableNode.push_back(gT[x][i]);
				
			}
		}

	}
	static inline double logcnk(const size_t n, size_t k)
	{
		k = k < n - k ? k : n - k;
		double res = 0;
		for (auto i = 1; i <= k; i++) res += log(double(n - k + i) / i);
		return res;
	}
	static inline double pow2(const double t)
	{
		return t * t;
	}
	double self_inf_cal_CP(vector<int>seedset, int RRnum)
	{
		std::vector<bool> vecBoolVst = std::vector<bool>(RRnum);
		std::vector<bool> vecBoolSeed(n);
		for (auto seed : seedset) vecBoolSeed[seed] = true;
		for (auto seed : seedset)
		{
			for (auto node : hyperG_2[seed])
			{
				vecBoolVst[node] = true;
			}
		}
		return 1.0 * std::count(vecBoolVst.begin(), vecBoolVst.end(), true);
	}

	double self_inf_cal_LRR(vector<int>seedset, int RRnum)
	{
		std::vector<bool> vecBoolVst = std::vector<bool>(RRnum);
		std::vector<bool> vecBoolSeed(n);
		for (auto seed : seedset) vecBoolSeed[seed] = true;
		for (auto seed : seedset)
		{
			for (auto node : LhyperG_2[seed])
			{
				vecBoolVst[node] = true;
			}
		}
		return 1.0 * std::count(vecBoolVst.begin(), vecBoolVst.end(), true);
	}
	
	int generate_one_RRset_WC(vector<int>isSeed, vector<int>isBlocker)
	{
		bool flag = false;//RR set�Ƿ񸲸�seed
		vector<int>tmp_visit_mark;//visit_mark
		tmp_visit_mark.resize(n + 1);
		vector<bool>tmp_visit;
		tmp_visit.resize(n + 1, 0);

		const auto uStart = sfmt_genrand_uint32(&sfmtSeed) % n;
		if (isBlocker[uStart] == 1)
			return 0;
		if (isSeed[uStart] == 1)
			return 1;
		unsigned int n_visit_mark = 0, curIdx = 0; //�ѷ��ʽ��ĸ�������ǰ���ʵ��±�
		tmp_visit_mark[n_visit_mark++] = uStart; //�����ѷ��ʵĽ��
		tmp_visit[uStart] = true; //����Ƿ񱻷��ʹ��ı��

		while (curIdx < n_visit_mark) {//SUBSIM
			int i = tmp_visit_mark[curIdx++];
			if ((int)gT_reverse[i].size() == 0) continue;
			double p = probT_reverse[i][0];
			if (p < 1)
			{
				double log2Prob = log2f(1 - p);
				double prob = (sfmt_genrand_real1(&sfmtSeed) * 0.999999) + 0.000001;
				int startPos = log2f(prob) / log2Prob;
				int endPos = (int)gT_reverse[i].size();
				while (startPos < endPos) {
					int v = gT_reverse[i][startPos];
					if (tmp_visit[v] || isBlocker[v] == 1) {
						int increment = log2f((sfmt_genrand_real1(&sfmtSeed) * 0.999999) + 0.000001) / log2Prob;
						startPos += (increment + 1);
						continue;
					}
					if (isSeed[v] == 1)
						return 1;
					tmp_visit[v] = true;
					tmp_visit_mark[n_visit_mark++] = v;
					int increment = log2f((sfmt_genrand_real1(&sfmtSeed) * 0.999999) + 0.000001) / log2Prob;
					startPos += increment + 1;
				}
			}
			else
			{
				int v = gT_reverse[i][0];
				if (tmp_visit[v] || isBlocker[v] == 1) continue;
				if (isSeed[v] == 1)
					return 1;
				tmp_visit[v] = true;
				tmp_visit_mark[n_visit_mark++] = v;
			}
		}
		return 0;
	}

	double estimate_inf_byStop(vector<int>nodeset, vector<int>blockerset, double epsilon, double delta)
	{
		vector<int>isSeed;
		isSeed.resize(n, 0);
		vector<int>isBlocker;
		isBlocker.resize(n, 0);
		for (int node : blockerset)
		{
			isBlocker[node] = 1;
		}
		for (int node : nodeset)
		{
			isSeed[node] = 1;
		}
		//double upperBound = 1 + (1 + epsilon) * 4 * (exp(1) - 2) * log(2 / delta) / (epsilon * epsilon);//naive threshold
		double T = 2 * (1 + epsilon) * (1 + 1 / 3.0 * epsilon) * log(2 / delta) / (epsilon * epsilon);//2-hop estimation
		cout << "Threshold: " << T << endl;
		long long sum = 0, num = 0;
		while (sum < T)
		{
			int flag = generate_one_RRset_WC(isSeed, isBlocker);//SUBSIM
			sum += flag;
			num++;//the number of RR set
		}
		cout << "RR num: " << num << endl;
		num = num - sum + T;
		return n * T / num * 1.0;
	}

	struct Node {
		int a;
		double b;

		Node(int a, double b) : a(a), b(b) {}
	};
	struct Compare {
		bool operator()(const Node& n1, const Node& n2) {
			return n1.b > n2.b;
		}
	};
	struct Compare1 {
		bool operator()(const Node& n1, const Node& n2) {
			return n1.b < n2.b;
		}
	};
	double calculate_OPT_lower(int budget, vector<int>CB)
	{
		vector<int>nodeset;
		for (int i = 0; i < CB.size(); i++)
		{
			int v = CB[i];
			influenced[v] = 1;
			for (int j = 0; j < gT_reverse[v].size(); j++)
			{
				int node = gT_reverse[v][j];
				influenced[v] *= (1 - probT_reverse[v][j]);
				//cout << influenced[v] << endl;
			}
			influenced[v] = 1 - influenced[v];
			influenced_deg[v] = influenced[v] * deg[v] * 1.0;
			//cout << influenced_deg[v] <<"  deg:" << deg[v] << endl;
		}
		priority_queue<Node, std::vector<Node>, Compare1> pq;
		for (int i = 0; i < CB.size(); i++) {
			pq.push(Node(CB[i], influenced[CB[i]]));
		}
		int seed;
		double OPT = 0;
		while (budget > 0)
		{
			seed = pq.top().a;
			pq.pop();
			if (isRumor[seed] == false)
			{
				budget--;
				OPT += influenced[seed];
				//cout << influenced[seed] << endl;
			}
		}
		return OPT;
	}


	int opimc_sandwich(const int targetSize, const double epsilon, const double delta, const Argument& arg, double inf, double OPT_L, int Uflag, int Lflag)
	{
		const double e = exp(1);
		const double approx = 1 - 1.0 / e;
		const double alpha = sqrt(log(12.0 / delta));
		const double beta = sqrt((1 - 1.0 / e) * (logcnk(n - rumorSet.size(), targetSize) + log(12.0 / delta)));
		size_t numRbase, maxNumR, numIter;

		numRbase = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
		maxNumR= size_t(2.0 * Rnode.size() * pow2((1 - 1 / e) * sqrt(log(6.0 / delta)) + sqrt((1 - 1.0 / e) * (logcnk(n - rumorSet.size(), targetSize) + log(6.0 / delta)))) / OPT_L / pow2(epsilon)) + 1;
		numIter = (size_t)log2(maxNumR / numRbase) + 1;
		const double a1 = log(numIter * 3.0 / delta);
		const double a2 = log(numIter * 3.0 / delta);

	
		long long CP_num = 0;
		long long subCP_num1 = 0;
		long long subCP_num2 = 0;
		cout << numRbase << " __ " << maxNumR << endl;
		//int Uflag = 0, Lflag = 0;
		int opim_num;
		for (auto idx = 0; idx < numIter; idx++)
		{
			//generate sample phase
			auto numR = numRbase << idx;
			int pre_subsize1 = generateCP_local_global(numR, arg, CP_num, subCP_num1, Uflag, Lflag);//CP1 and LRR1
			subCP_num1 = pre_subsize1;
			int pre_subsize2 = generateCP_local_global_2(numR, arg, CP_num, subCP_num2, Uflag, Lflag);//CP2 and LRR2
			subCP_num2 = pre_subsize2;
			CP_num = numR;//presize
			cout << "Sample num: " << numR << endl;
			//test the quality of seeds phase, (i) LB (ii) UB
			if (Lflag == 0)
			{
				auto boundMin = select_seed_low(arg);
				auto upperBound = boundMin * (1 + arg.beta) / inf;

				int CPset2_num = hyperGT_2.size();
				auto infVldt = self_inf_cal_CP(LB_seedSet, CPset2_num);

				if (infVldt * (1 - arg.beta) / inf >= 5.0 * a1 / 18)
					infVldt = infVldt * (1 - arg.beta) / inf;
				else if (infVldt * (1 + arg.beta) / inf <= 5.0 * a1 / 18)
					infVldt = infVldt * (1 + arg.beta) / inf;
				else
					infVldt = 0;
				auto lowerSelect = pow2(sqrt(infVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
				auto upperOPT = pow2(sqrt(upperBound + a2 / 2.0) + sqrt(a2 / 2.0));
				auto approxOPIMC = lowerSelect / upperOPT;
				cout << "lower ratio:" << approxOPIMC << endl;
				if (approxOPIMC >= approx - epsilon)
				{
					cout << "RR_num of lower bound:" << numR << " max_num:" << maxNumR << endl;
					Lflag = 1;
					opim_num = numR;
					disp_mem_usage();
					relax_memory(hyperG);
					relax_memory(hyperGT);
					relax_memory(hyperG_2);
					relax_memory(hyperGT_2);
				}
			}
			if (Uflag == 0)
			{
				auto boundMin = select_seed_upper(arg);
				double upperBound;
				int infVldt;
				double lowerSelect, upperOPT;

				int LRR2_num = LhyperGT_2.size();
				infVldt = self_inf_cal_LRR(UB_seedSet, LRR2_num);
				upperBound = boundMin;

				lowerSelect = pow2(sqrt(infVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
				upperOPT = pow2(sqrt(upperBound + a2 / 2.0) + sqrt(a2 / 2.0));
				

				auto approxOPIMC = lowerSelect / upperOPT;
				cout << "upper ratio:" << approxOPIMC << endl;
				if (approxOPIMC >= approx - epsilon)
				{
					cout << "RR_num of upper bound: " << numR << endl;
					Uflag = 1;
				}
			}
			if (Uflag == 1 && Lflag == 1)
				break;
		}
		return opim_num;
	}

	void relax_memory(vector<vector<int>> myVector)
	{
		for (auto& innerVector : myVector) {
			innerVector.clear();
			vector<int>().swap(innerVector);
		}

		myVector.clear();
		vector<vector<int>>().swap(myVector);
	}
	void deg_based_heuristic(int budget, vector<int>CB)
	{
		vector<int>nodeset;

		priority_queue<Node, std::vector<Node>, Compare1> pq;
		for (int i = 0; i < CB.size(); i++) {
			pq.push(Node(CB[i], influenced_deg[CB[i]]));
		}
		int seed;
		while(budget>0)
		{
			seed = pq.top().a;
			pq.pop();
			if (isRumor[seed] == false)
			{
				Or_seedSet.push_back(seed);
				budget--;
				//cout << influenced_deg[seed] << endl;
			}
		}
	}


	struct Node2 {// min-heap
    int id;
    int activate_time;
    bool operator<(const Node2 &other) const {
        return activate_time > other.activate_time; // Min-heap
    }
	};

	double generate_TRU_set_with_blocker(vector<int>&blocker,int deadline)
	{
		vector<int> isBlocker(n, 0);
		for (int node : blocker) {
			isBlocker[node] = 1;
		}
		int source_node = sfmt_genrand_uint32(&sfmtSeed) % n; 
		if(isRumor[source_node]) return computeU(source_node, 0);
		if(isBlocker[source_node]) return 0;

		priority_queue<Node2> Q;
		vector<int> AT(n, 10000);
		vector<int> I(n, 0);

		AT[source_node] = 0;
		Q.push({source_node, 0});
		
		while (!Q.empty()) {
			Node2 current = Q.top();
			Q.pop();

			int q = current.id;
			int q_AT = current.activate_time;
			if (I[q] == 1 || isBlocker[q]) continue; 
			if(isRumor[q]) return computeU(source_node,AT[q]);
			I[q] = 1;

			for (int i = 0; i < gT_reverse[q].size(); i++) {
				int u = gT_reverse[q][i];
				if (isBlocker[u]) continue;

				int tmp= sfmt_genrand_uint32(&sfmtSeed) % precomputedGeometric[u].size();
				int diff_delay=precomputedGeometric[u][tmp]; 
				

				if (sfmt_genrand_real1(&sfmtSeed) <= probT_reverse[q][i] && AT[q] + diff_delay <= deadline) {
					if (AT[q] + diff_delay < AT[u]) {
						AT[u] = AT[q] + diff_delay;
						Q.push({u, AT[u]});
					}
				}
			}	
		}
			
		return 0;
	}

	double stopping_rule(vector<int>& blocker, int deadline, double epsilon, double delta)
	{
		double x_max=1;
		double T = 2 * x_max * (1 + epsilon) * (1 + 1 / 3.0 * epsilon) * log(2 / delta) / (epsilon * epsilon);

		double sum = 0, theta = 0;
		double flag;
		while (sum < T)
		{
			flag = generate_TRU_set_with_blocker(blocker, deadline);
			sum += flag;
			theta++;
		}
	
		theta = theta - (sum - T)/(double)flag;
		return n * T / (double) theta * 1.0;
	}

	// /// Utility funtion
	// int U(int x) {
	// 	if(x<0) return 0;
	// 	return deadline_1/(double)pow((deadline_1-x+1),1);
	// 	// return std::exp(x/32.0)+1;
	// }
	double computeU(int v, double t_v) {
		double x_v = stratifiedVotingTimes[v];
		double B= 0.1 * deadline_1;
		// 计算分段函数 U_v(T - t_v)
		if (0 <= t_v && t_v <= std::max(x_v - B, 0.0)) {
			return 1.0;
		} else if (std::max(x_v - B, 0.0) < t_v && t_v <= std::min(x_v + B, static_cast<double>(deadline_1))) {
			return 1.0 - (t_v - (x_v - B)) / (2.0 * B);
		} else if (t_v > std::min(x_v + B, static_cast<double>(deadline_1))) {
			return 0.0;
		}
		
		return 0.0; // 默认情况
	}

	double output_AD_influence(Argument& arg)
	{
		int T=arg.T;
		double gamma=arg.gamma;
		double delta=1.0/(double)n;
		double inf_after_block = stopping_rule(LB_seedSet, T, gamma, delta);
		// cout << "Estimated AD influence after blocking: " << inf_after_block - rumorSet.size() * U(deadline_1) << endl;
		return inf_after_block;
	}


	double MC_based_estimate_with_removal_2(int count, int deadline, const vector<int>& removedNodes) {
		int cnt = 0;
		double sum = 0;
		vector<int> isBlocker(n, 0);
		for (int node : removedNodes) {
			isBlocker[node] = 1;
		}

		while (cnt < count) {
			priority_queue<Node2> Q;
			vector<int> AT(n, 10000);
			vector<int> I(n, 0);

			for (int u : rumorSet) {
				AT[u] = 0;
				Q.push({u, 0});
			}

			while (!Q.empty()) {
				Node2 current = Q.top();
				Q.pop();

				int q = current.id;
				int q_AT = current.activate_time;

				if (I[q] == 1 || isBlocker[q]) continue; 

				for (int i = 0; i < gT[q].size(); i++) {
					int tmp= sfmt_genrand_uint32(&sfmtSeed) % precomputedGeometric[q].size();
					int diff_delay=precomputedGeometric[q][tmp]; 
					int u = gT[q][i];
					if (isBlocker[u]) continue;

					if (sfmt_genrand_real1(&sfmtSeed) <= probT[q][i] && AT[q] + diff_delay <= deadline) {
						if (AT[q] + diff_delay < AT[u]) {
							AT[u] = AT[q] + diff_delay;
							Q.push({u, AT[u]});
						}
					}
				}

				I[q] = 1;
				if(AT[q]<=stratifiedVotingTimes[q]) sum++;
				// sum+=U(deadline-AT[q]);
			}
			cnt++;
		}

		return sum / (double)count;
	}

	void precomputeGeometricForAllNodes(int sampleSize) {
		static thread_local std::mt19937 gen(std::random_device{}());
		precomputedGeometric.clear();
		precomputedGeometric.resize(n);
		for (int u = 0; u < n; ++u) {
			int out_deg = deg[u];
			double p = 5.0 / (out_deg + 5.0);
			std::geometric_distribution<> geom(p);
			precomputedGeometric[u].reserve(sampleSize);
			for (int i = 0; i < sampleSize; ++i) {
				precomputedGeometric[u].push_back(geom(gen) + 1); // +1 保证最小为1
			}
		}
	}

	void loadVotingTimesFromFile(const std::string& filename) {
		std::ifstream inFile(filename);
		if (!inFile.is_open()) {
			std::cerr << "无法打开文件: " << filename << std::endl;
			return;
		}

		stratifiedVotingTimes.clear();
		stratifiedVotingTimes.resize(n);
		
		double time;
		int v = 0;
		while (inFile >> time && v < n) {
			stratifiedVotingTimes[v] = time;
			v++;
		}

		inFile.close();
		std::cout << "已从文件加载 " << v << " 个节点的投票时间" << std::endl;
	}

	void loadRealVotingTimesFromFile(const std::string& filename) {
		std::ifstream inFile(filename);
		if (!inFile.is_open()) {
			std::cerr << "无法打开文件: " << filename << std::endl;
			return;
		}
		realVotingTimes.clear();
		double t;
		while (inFile >> t) {
			realVotingTimes.push_back(t);
		}
		inFile.close();
		std::cout << "已从文件加载 " << realVotingTimes.size() << " 个节点的真实投票时间" << std::endl;
	}

};