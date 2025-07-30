#include "iheap.h"
#include <queue>	//priority_queue
#include <utility>  // pair
#include <numeric>
#include <algorithm>
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
	double MC_based_estimate(vector<int> rumorNode, int count)
	{
		vector<int> vis;
		vis.resize(n);
		double sum = 0;
		for (int tag = 1; tag <= count; tag++)
		{
			//vector<vector<int>> eSample;
			//eSample.resize(n);
			queue<int> q;
			for (int x : rumorNode)
				q.push(x), vis[x] = -tag, sum += 1.0;
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
						//eSample[x].push_back(gT[x][i]);
						if (vis[gT[x][i]] != -tag)
							q.push(gT[x][i]), vis[gT[x][i]] = -tag, sum += 1.0, C[gT[x][i]]++;
					}
				}
			}
		}
		return sum / ((double)count);
	}


	void Delete_Node(vector<int> Nodeset)
	{
		for (int node : Nodeset)//forward
		{
			for (int i = 0; i < (int)gT_reverse[node].size(); i++)
			{
				// reverse ..
				probT_reverse[node][i] = 0;

				//正的概率 set to 0
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
		probT_reverse.assign(probT_reverse2.begin(),probT_reverse2.end());
	}

	void gen_sample_2(vector<vector<int>>& sample,vector<vector<int>>& sample_UB, vector<int>rumorset, vector<int>& ReachableNode, vector<int>& vis, int tag) {

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
						q.push(gT[x][i]), vis[gT[x][i]] = -1, ReachableNode.push_back(gT[x][i]), sample_UB[x].push_back(gT[x][i]);
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

	double select_seed_up(const Argument& arg)
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
			if (k == 0)
			{
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
			// compute(rt);
		}
		void construct_CPset(int uStart, vector<vector<int>>& hyperG, vector<vector<int>>& hyperGT, vector<bool>isRumor, int & tag)
		{
			int current_v=uStart;
			hyperG[current_v].push_back(tag);
			hyperGT[tag].push_back(current_v);
			while(isRumor[idom[current_v]]==0){
				current_v = idom[current_v];
				hyperG[current_v].push_back(tag);
				hyperGT[tag].push_back(current_v);
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

	void generateCP_local_global(int count, const Argument& arg, int presize, int Uflag, int Lflag)
	{
		int num = presize;
		while (num < count)
		{
			ReachableNode.clear();
			vector<vector<int>> sample;
			sample.resize(n);

			vector<vector<int>> sample_UB;
			sample_UB.resize(n);
			
			vector<vector<int>> sample_UB_reverse;
			sample_UB_reverse.resize(n);

			vector<int> vis;
			vis.resize(n, 0);
			gen_sample_2(sample, sample_UB, rumorSet, ReachableNode, vis, num);
			if (Uflag == 0)
			{
				for (int i = 0; i < n; i++)
				{
					for (int v : sample_UB[i])
						sample_UB_reverse[v].push_back(i);
				}
				LhyperGT.push_back(vector<int>());
				const auto random_num = sfmt_genrand_uint32(&sfmtSeed) % ReachableNode.size(); 
				int uStart = ReachableNode[random_num];
				int current_v = uStart;
				LhyperG[current_v].push_back(num);
				LhyperGT[num].push_back(current_v);

				while(isRumor[sample_UB_reverse[current_v][0]]==0){
					current_v = sample_UB_reverse[current_v][0];
					LhyperG[current_v].push_back(num);
					LhyperGT[num].push_back(current_v);
				}
			}
			if (Lflag == 0)
			{
				sample.push_back(rumorSet);
				tree.init(n, sample);
				tree.go(n);
				hyperGT.push_back(vector<int>());
				const auto random_num = sfmt_genrand_uint32(&sfmtSeed) % ReachableNode.size(); 
				int uStart = ReachableNode[random_num];
				tree.construct_CPset(uStart, hyperG, hyperGT, isRumor, num);
			}
			num++;
		}
		//return j;
	}
	void generateCP_local_global_2(int count, const Argument& arg, int presize, int Uflag, int Lflag)
	{
		int num = presize;
		while (num < count)
		{
			ReachableNode.clear();
			vector<vector<int>> sample;
			sample.resize(n);

			vector<vector<int>> sample_UB;
			sample_UB.resize(n);
			
			vector<vector<int>> sample_UB_reverse;
			sample_UB_reverse.resize(n);

			vector<int> vis;
			vis.resize(n, 0);
			gen_sample_2(sample, sample_UB, rumorSet, ReachableNode, vis, num);
			if (Uflag == 0)
			{
				for (int i = 0; i < n; i++)
				{
					for (int v : sample_UB[i])
						sample_UB_reverse[v].push_back(i);
				}
				LhyperGT_2.push_back(vector<int>());
				const auto random_num = sfmt_genrand_uint32(&sfmtSeed) % ReachableNode.size(); 
				int uStart = ReachableNode[random_num];
				int current_v = uStart;
				LhyperG_2[current_v].push_back(num);
				LhyperGT_2[num].push_back(current_v);

				while(isRumor[sample_UB_reverse[current_v][0]]==0){
					current_v = sample_UB_reverse[current_v][0];
					LhyperG_2[current_v].push_back(num);
					LhyperGT_2[num].push_back(current_v);
				}
			}
			if (Lflag == 0)
			{
				sample.push_back(rumorSet);
				tree.init(n, sample);
				tree.go(n);
				hyperGT_2.push_back(vector<int>());
				const auto random_num = sfmt_genrand_uint32(&sfmtSeed) % ReachableNode.size(); 
				int uStart = ReachableNode[random_num];
				tree.construct_CPset(uStart, hyperG_2, hyperGT_2, isRumor, num);
			}
			num++;
		}

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
	
	int generate_one_RRset_stop(vector<int>& isRumor)
	{
		vector<int>tmp_visit_mark;//visit_mark
		tmp_visit_mark.resize(n + 1);
		vector<bool>tmp_visit;
		tmp_visit.resize(n + 1, 0);

		const auto uStart = sfmt_genrand_uint32(&sfmtSeed) % n; 
		unsigned int n_visit_mark = 0, curIdx = 0;
		tmp_visit_mark[n_visit_mark++] = uStart; 
		tmp_visit[uStart] = true; 
		if(isRumor[uStart]==1){
			//cout<<"222"<<endl;
			return 1;
		}
		//cout<<curIdx<<"  "<<n_visit_mark<<endl;
		while (curIdx < n_visit_mark) {
			int i = tmp_visit_mark[curIdx++]; 
			for (int j = 0; j < (int)gT_reverse[i].size(); j++) { 
				//cout<<"node is:"<<i<<endl;
				int v = gT_reverse[i][j]; 
				if (tmp_visit[v])continue; 
				double randDouble = sfmt_genrand_real1(&sfmtSeed); 
				if (randDouble > probT_reverse[i][j])continue; 
				tmp_visit[v] = true; 
				tmp_visit_mark[n_visit_mark++] = v; 
				if(isRumor[v]==1){
					//cout<<"000"<<endl;
					return 1;
				}
			}
		}
		// for (unsigned int i = 0; i < n_visit_mark; i++)visit[visit_mark[i]] = false;
		//cout<<"11"<<endl;
		return 0;
	}

	// double estimate_inf_byStop(vector<int>nodeset, vector<int>blockerset, double epsilon, double delta)
	// {
	// 	vector<int>isSeed;
	// 	isSeed.resize(n, 0);
	// 	vector<int>isBlocker;
	// 	isBlocker.resize(n, 0);
	// 	for (int node : blockerset)
	// 	{
	// 		isBlocker[node] = 1;
	// 	}
	// 	for (int node : nodeset)
	// 	{
	// 		isSeed[node] = 1;
	// 	}
	// 	//double upperBound = 1 + (1 + epsilon) * 4 * (exp(1) - 2) * log(2 / delta) / (epsilon * epsilon);//naive threshold
	// 	double T = 2 * (1 + epsilon) * (1 + 1 / 3.0 * epsilon) * log(2 / delta) / (epsilon * epsilon);//2-hop estimation
	// 	cout << "Threshold: " << T << endl;
	// 	long long sum = 0, num = 0;
	// 	while (sum < T)
	// 	{
	// 		int flag = generate_one_RRset_WC(isSeed, isBlocker);//SUBSIM
	// 		sum += flag;
	// 		num++;//the number of RR set
	// 	}
	// 	cout << "RR num: " << num << endl;
	// 	num = num - sum + T;
	// 	return n * T / num * 1.0;
	// }
	double estimate_inf_byStop(vector<int> & rumorset, double epsilon, double delta)
	{
		vector<int>isRumor;
		isRumor.resize(n, 0);
		for (int node : rumorset)
		{
			isRumor[node] = 1;
		}
		//double upperBound = 1 + (1 + epsilon) * 4 * (exp(1) - 2) * log(2 / delta) / (epsilon * epsilon);//naive threshold
		double T = 2 * (1 + epsilon) * (1 + 1 / 3.0 * epsilon) * log(2 / delta) / (epsilon * epsilon);//2-hop estimation
		//cout << "Threshold: " << T << endl;
		long long sum = 0, num = 0;
		while (sum < T)
		{
			int flag = generate_one_RRset_stop(isRumor);
			sum += flag;
			num++;//the number of RR set
		}
		//cout<<"sum: "<<sum<<endl;
		//cout << "RR num: " << num << endl;
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


	int opimc_sandwich(const int targetSize, const double epsilon, const double delta, const Argument& arg, double inf, double OPT_L, int Uflag, int Lflag, double & cov_ub)
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
		cout << numRbase << " __ " << maxNumR << endl;
		//int Uflag = 0, Lflag = 0;
		int opim_num;
		for (auto idx = 0; idx < numIter; idx++)
		{
			//generate sample phase
			auto numR = numRbase << idx;
			generateCP_local_global(numR, arg, CP_num, Uflag, Lflag);//CP1 and LRR1
			
			generateCP_local_global_2(numR, arg, CP_num, Uflag, Lflag);//CP2 and LRR2
			
			CP_num = numR;//presize
			cout << "Sample num: " << numR << endl;
			//test the quality of seeds phase, (i) LB (ii) UB
			if (Lflag == 0)
			{
				auto boundMin = select_seed_low(arg);
				auto upperBound = boundMin;

				int CPset2_num = hyperGT_2.size();
				auto infVldt = self_inf_cal_CP(LB_seedSet, CPset2_num);
				//cout<<"lb_cov: "<<infVldt<<endl;

				auto lowerSelect = pow2(sqrt(infVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
				auto upperOPT = pow2(sqrt(upperBound + a2 / 2.0) + sqrt(a2 / 2.0));
				auto approxOPIMC = lowerSelect / upperOPT;
				//cout << "lower ratio:" << approxOPIMC << endl;
				if (approxOPIMC >= approx - epsilon)
				{
					cout << "RR_num of lower bound:" << numR << " max_num:" << maxNumR << endl;
					Lflag = 1;
					opim_num = numR;
					// disp_mem_usage();
					// relax_memory(hyperG);
					// relax_memory(hyperGT);
					// relax_memory(hyperG_2);
					// relax_memory(hyperGT_2);
				}
			}
			if (Uflag == 0)
			{
				auto boundMin = select_seed_up(arg);
				auto upperBound = boundMin;

				int LRRset2_num = LhyperGT_2.size();
				auto infVldt = self_inf_cal_LRR(UB_seedSet, LRRset2_num);
				//cout<<"ub_cov: "<<infVldt<<endl;
			
				auto lowerSelect = pow2(sqrt(infVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
				auto upperOPT = pow2(sqrt(upperBound + a2 / 2.0) + sqrt(a2 / 2.0));
				auto approxOPIMC = lowerSelect / upperOPT;
				//cout << "UPPER ratio:" << approxOPIMC << endl;
				if (approxOPIMC >= approx - epsilon)
				{
					cout << "sample_num of upper bound:" << numR << " max_num:" << maxNumR << endl;
					Uflag = 1;
					opim_num = numR;
					cov_ub= infVldt/(double)LRRset2_num;
					// cout<<"cov_ub is: "<<cov_ub<<endl;
					// disp_mem_usage();
					// relax_memory(hyperG);
					// relax_memory(hyperGT);
					// relax_memory(hyperG_2);
					// relax_memory(hyperGT_2);
				}
				// auto boundMin = select_seed_upper(arg);
				// double upperBound;
				// int infVldt;
				// double lowerSelect, upperOPT;

				// int LRR2_num = LhyperGT_2.size();
				// infVldt = self_inf_cal_LRR(UB_seedSet, LRR2_num);
				// upperBound = boundMin;

				// lowerSelect = pow2(sqrt(infVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
				// upperOPT = pow2(sqrt(upperBound + a2 / 2.0) + sqrt(a2 / 2.0));
				

				// auto approxOPIMC = lowerSelect / upperOPT;
				// cout << "upper ratio:" << approxOPIMC << endl;
				// if (approxOPIMC >= approx - epsilon)
				// {
				// 	cout << "RR_num of upper bound: " << numR << endl;
				// 	Uflag = 1;
				// }
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

	// report approximation
	int generate_lrr_stop(vector<int> & isUB){

		ReachableNode.clear();
		vector<vector<int>> sample;
		sample.resize(n);

		vector<vector<int>> sample_UB;
		sample_UB.resize(n);
		
		vector<vector<int>> sample_UB_reverse;
		sample_UB_reverse.resize(n);

		vector<int> vis;
		vis.resize(n, 0);
		gen_sample_2(sample, sample_UB, rumorSet, ReachableNode, vis, 0);
		
		for (int i = 0; i < n; i++)
		{
			for (int v : sample_UB[i])
				sample_UB_reverse[v].push_back(i);
		}
		
		const auto random_num = sfmt_genrand_uint32(&sfmtSeed) % ReachableNode.size(); 
		int uStart = ReachableNode[random_num];
		int current_v = uStart;
		if(isUB[current_v]==1) return 1;


		// LhyperG[current_v].push_back(num);
		// LhyperGT[num].push_back(current_v);

		while(isRumor[sample_UB_reverse[current_v][0]]==0){
			current_v = sample_UB_reverse[current_v][0];
			if(isUB[current_v]==1) return 1;
		}
			
		return 0;	
	}

	double estimate_upper_byStop(Argument& arg, double rumor_inf){
		vector<int>isUB;
		isUB.resize(n, 0);
		for (int node : UB_seedSet)
		{
			isUB[node] = 1;
		}
		double delta = 1.0 / n;
		double T = 2 * (1 + arg.gamma) * (1 + 1 / 3.0 * arg.gamma) * log(2 / delta) / (arg.gamma * arg.gamma);
		//cout << "Threshold: " << T << endl;
		long long sum = 0, num = 0;
		while (sum < T)
		{
			int flag = generate_lrr_stop(isUB);
			sum += flag;
			num++;//the number of RR set
		}
		//cout<<"sum: "<<sum<<endl;
		//cout << "RR num: " << num << endl;
		num = num - sum + T;
		return rumor_inf * T / num * 1.0;
	}
	int one_lrr_sequence_stop(vector<int>& isUB){

		
		ReachableNode.clear();
		vector<vector<int>> sample;
		sample.resize(n);

		vector<vector<int>> sample_UB;
		sample_UB.resize(n);
		
		vector<vector<int>> sample_UB_reverse;
		sample_UB_reverse.resize(n);

		vector<int> vis;
		vis.resize(n, 0);
		gen_sample_2(sample, sample_UB, rumorSet, ReachableNode, vis, 0);
		int remain;
		queue<int> q;
		int cnt=0;
		vector<int> vis2;
		vis2.resize(n, 0);
		for (auto x : UB_seedSet)
			q.push(x), vis2[x]=1, cnt++;
		while (!q.empty())
		{
			int x = q.front();
			q.pop();
			for (int i = 0; i < sample_UB[x].size(); i++)
			{
				int v=sample_UB[x][i];
				if(vis2[x]==0){
					q.push(v);
					cnt++;
					vis2[x]=1;
				}		
			}
		}
		remain = ReachableNode.size()+rumorSet.size()-cnt;
		//cout<<"remain:  "<<remain<<endl;
		return remain;
	}

	double estimate_upper(){
		double num=30000;
		vector<int>isUB;
		isUB.resize(n, 0);
		for (int node : UB_seedSet)
		{
			isUB[node] = 1;
		}
		int cnt=0;
		for(int i=0;i<num;i++){
			int tmp = one_lrr_sequence_stop(isUB);
			cnt+=tmp;
		}
		//cout<<"cnt: "<<cnt<<endl;
		double out=cnt/(double)num;
		return out;
	}

	double report_approximation_ratio(Argument& arg, double rumor_inf){
		double tmp=pow((1-arg.gamma),2)/(double)(pow((1+arg.gamma),2));
		const double e = exp(1);
		tmp=tmp*(1-1/e-arg.epsilon);
		cout<<"tmp: "<<tmp<<endl;
		// UB_seedSet
		
		// double Dsu_bu = estimate_upper_byStop(arg, rumor_inf);
		
		// double Dsu_bu = coverage * rumor_inf;
		double remain_ub = estimate_upper();
		
		double Dsu_bu = rumor_inf-remain_ub;
		cout<<"Dsu_bu: "<<Dsu_bu<<endl;

		reset_pro();
		Delete_Node(UB_seedSet);
		double Ds_bu = estimate_inf_byStop(rumorSet, arg.gamma, 1.0 / n);
		Ds_bu = rumor_inf-Ds_bu;
		cout<<"Ds_bu"<<Ds_bu<<endl;
		
		tmp = tmp * Ds_bu/Dsu_bu;
		
		cout<<"approximation ratio: "<<tmp<<endl;
		return tmp;
	}
	// to test the upper bound 一定是 upper bound
	void test(){
		ReachableNode.clear();

			vector<vector<int>> sample;
			sample.resize(n);

			vector<vector<int>> sample_UB;
			sample_UB.resize(n);
			
			vector<vector<int>> sample_UB_reverse;
			sample_UB_reverse.resize(n);

			vector<int> vis;
			vis.resize(n, 0);
			gen_sample_2(sample, sample_UB, rumorSet, ReachableNode, vis, 0);

			vector<int>sum_ub(n,0);
			vector<int>sum(n,0);
			// upper bound
			for(int i=0;i<ReachableNode.size();i++){
				int v=ReachableNode[i];
				if(isRumor[v]==1) continue;
				sum_ub[v]++;
				queue<int> q;
				q.push(v);
				while (!q.empty())
				{
					int x = q.front();
					q.pop();
					for (int i = 0; i < sample_UB[x].size(); i++)
					{
						int u =sample_UB[x][i];
						q.push(u),sum_ub[v]++;
					}
				}
			}
			

			
			sample.push_back(rumorSet);
			tree.init(n, sample);
			tree.go(n);
			for (int i = 0; i < n; i++)
				sum[i] += (double)tree.sum[i];
			int count=0;

			for(int i=0;i<ReachableNode.size();i++){
				int v=ReachableNode[i];
				if(isRumor[v]==1) continue;
				cout<<"upper bound: "<<sum_ub[v]<<"---- original: "<<sum[v]<<endl;
				if(sum_ub[v]>sum[v])count++;
				if(sum_ub[v]<sum[v])cout<<"fuck!"<<endl;
			}
			cout<<"the number of larger is : "<<count<<endl;
	}
	double rep_approximation(Argument& arg){
		vector<int>isUB;
		isUB.resize(n, 0);
		for (int node : UB_seedSet)
		{
			isUB[node] = 1;
		}
		double sum_ori=0;
		double sum_up=0;
		for(int i=0;i<50000;i++){
			ReachableNode.clear();
			vector<vector<int>> sample;
			sample.resize(n);

			vector<vector<int>> sample_UB;
			sample_UB.resize(n);
			
			vector<vector<int>> sample_UB_reverse;
			sample_UB_reverse.resize(n);

			vector<int> vis;
			vis.resize(n, 0);
			gen_sample_2(sample, sample_UB, rumorSet, ReachableNode, vis, 0);

			int decrease_upper=0, decrease_ori=0;
			queue<int> q;
			int activated_num_ori=0;//表示没有ub以后ori function的值
			vector<int> vis3;
			vis3.resize(n, 0);
			for (auto x : rumorSet)
				q.push(x), vis3[x]=1;
			while (!q.empty())
			{
				int x = q.front();
				q.pop();
				for(int v:sample[x]){
					if(isUB[v]==1 || vis3[v]==1)continue;
					
					q.push(v);
					activated_num_ori++;
					vis3[v]=1;
					
				}		
			}
			decrease_ori=ReachableNode.size()-activated_num_ori;
			//cout<<"decrease_ori: "<<decrease_ori<<endl;
			sum_ori+=decrease_ori;
			
		
			vector<int> vis2;
			vis2.resize(n, 0);
			for (auto x : UB_seedSet)
				q.push(x), vis2[x]=1, decrease_upper++;
			while (!q.empty())
			{
				int x = q.front();
				q.pop();
				for(int v:sample_UB[x]){
					if(vis2[v]==0){
						q.push(v);
						decrease_upper++;
						vis2[v]=1;
					}
				}
	
			}
			//cout<<"decrease_upper is: "<<decrease_upper<<endl;
			sum_up+=decrease_upper;
			
			if(decrease_ori>decrease_upper){
				cout<<"fk!"<<endl;
				break;
			}
		}

		double t1=sum_ori/(double)50000;

		double t2=sum_up/(double)50000;

		double tmp=pow((1-arg.gamma),2)/(double)(pow((1+arg.gamma),2));
		const double e = exp(1);
		tmp=tmp*(1-1/e-arg.epsilon)*t1/t2;
		cout<<"approximation: "<<tmp<<endl;


		return tmp;
	}
};