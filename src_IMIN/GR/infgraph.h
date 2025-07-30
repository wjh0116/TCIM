#include "iheap.h"
#include <queue>	//priority_queue
#include <utility>  // pair
#include <numeric>

class InfGraph: public Graph
{
private:
	//vector<bool> activated;
    vector<bool> visit;
    vector<int> visit_mark;
	
public:
    vector<vector<int>> hyperG;
    vector<vector<int>> hyperGT;
	vector<vector<int>> hyperG_2;
	vector<vector<int>> hyperGT_2;
	vector<int> ActNode;

	int S;
	vector<int> Dfn; // Dfn[i]表示节点i的dfs序
	vector<int> Ord; // Ord[i]表示排名第i位的节点
	vector<int> Parent; // Parent[i]表示深搜树上i节点的父节点
	int Stamp;
	vector<int> uni; // 带权并查集的father数组
	vector<int> mn;  // mn[i]表示i的所有祖先中，半支配点的dfs序最小的那个祖先
	vector<int>rank;
	vector<int> Sdom; // Sdom[i]表示i的半支配点
	vector<int> Idom; // Idom[i]表示i的直接支配点
	vector<vector<int>> SdomTree;
	vector<int> ans;
	vector<double>D;
	sfmt_t sfmtSeed;
	vector<int> seedSet;
	vector<int> rumorSet;
	vector<bool> isRumor;
	vector<bool> isSelect;
	vector<int> Dec;//选择每个节点以后的activity的减少值

    InfGraph(string folder, string graph_file): Graph(folder, graph_file)
    {
        srand(time(NULL));
        sfmt_init_gen_rand(&sfmtSeed , rand());		
        visit = vector<bool> (n+1);
        visit_mark = vector<int> (n+1);
		Dec = vector<int>(n+1);
		isRumor = vector<bool>(n+1);
		isSelect = vector<bool>(n+1);
		//activated = vector<bool>(n, false);
		hyperG.resize(n+1, vector<int>());		
		hyperG_2.resize(n+1, vector<int>());
		ans = vector<int>(n + 2);
		D = vector<double>(n + 1);
    }

    void init_hyper_graph(){
		for (auto& hyper : hyperG)hyper.clear();
		for (auto& hyperT : hyperGT)vector<int>().swap(hyperT);
		hyperGT.clear();
		for (auto& hyper : hyperG_2)hyper.clear();
		for (auto& hyperT : hyperGT_2)vector<int>().swap(hyperT);
		hyperGT_2.clear();
		seedSet.clear();
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
	double MC_based_estimate(vector<int> rumorNode, int count, const Argument& arg)
	{
		vector<int> vis;
		vis.resize(n);
		double sum = 0;
		for (int tag = 1; tag <= count; tag++)
		{
			vector<vector<int>> eSample;
			eSample.resize(n);
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
						eSample[x].push_back(gT[x][i]);
						if (vis[gT[x][i]] != -tag)
							q.push(gT[x][i]), vis[gT[x][i]] = -tag, sum += 1.0;
					}
				}
			}
		}
		return sum / ((double)count);
	}

	void Delete_Node(vector<int> Nodeset)
	{
		for (int node : Nodeset)//5
		{
			for (int i = 0; i < (int)gT_reverse[node].size(); i++)//node的全部入边的概率设置为0
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
			compute(rt);
		}
	} tree;
	/*void Lengauer_Tarjan(int s) {
		
		/// 确定dfs序
		Dfn.assign(n + 3, Stamp = 0);
		Ord.assign(n + 3, 0);
		Parent.assign(n + 3, 0);
		dfs(S = s);

		/// 求解半支配点与支配点的中间答案
		Sdom.assign(n + 3, 0);
		uni.assign(n + 3, 0);
		rank.assign(n + 3, 0);
		mn.assign(n + 3, 0);
		SdomTree.assign(n + 3, vector<int>());
		Idom.assign(n + 3, -1);
		calcSdom();

		/// 最后确定支配点
		for (int i = 2; i <= Stamp; ++i) {
			int u = Ord[i];
			if (Idom[u] != Sdom[u])
				Idom[u] = Idom[Idom[u]];
			SdomTree[Idom[u]].push_back(u);
		}
		return;
	}
	void dfs(int u) {
		Ord[Dfn[u] = ++Stamp] = u;//ord起码n+2，dfn起码n+1
		for (auto v : sample[u]) {
			if (0 == Dfn[v]) {
				Parent[v] = u;
				dfs(v);
			}
		}
		return;
	}
	void calcSdom() {
		for (int i = 0; i <= n; ++i) Sdom[i] = mn[i] = uni[i] = i;//5.21
		for (int i = Stamp; i >= 2; --i) {//5.21
			int u = Ord[i]; // 排名第i位的节点是u
			int tmp = Stamp + 1;
			for (int v : sample_reverse[u]) {//u的入节点
				if (0 == Dfn[v]) continue;//没被DFS遍历，但是一般不会有这种情况
				uni_query(v);
				tmp = min(tmp, Dfn[Sdom[mn[v]]]);
			}
			uni[u] = Parent[u]; // 将u向上合并   1.25->24

			SdomTree[Sdom[u] = Ord[tmp]].push_back(u);
			u = Ord[i - 1];
			for (int& v : SdomTree[u])
			{
				uni_query(v);
				if (Sdom[mn[v]] != u)
					Idom[v] = mn[v];
				else
					Idom[v] = u;
			}
			SdomTree[u].clear();
		}
		return;
	}



	int uni_query(int u) {
		if (u == uni[u]) return u;
		int ans = uni_query(uni[u]);
		if (Dfn[Sdom[mn[uni[u]]]] < Dfn[Sdom[mn[u]]]) mn[u] = mn[uni[u]];
		if (rank[u] < rank[uni[u]]) swap(u, uni[u]); // 按秩合并
		rank[u] += rank[uni[u]];
		return uni[u] = ans;
	}*/


	void gen_sample(vector<vector<int>>& sample, vector<int>rumorset) {

		queue<int> q;
		vector<int> vis;
		vis.resize(n, 0);
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
						q.push(gT[x][i]), vis[gT[x][i]] = -1;
				}
			}
		}
	}

	int DecreaseES(int count, vector<int>rumorset, vector<bool> & remove_flag)
	{
		D.assign(n + 1, 0);
		vector<double> sum;
		sum.resize(n, 0);
		for (int j = 0; j < count; j++)//j代表了CP set的编号
		{
			//disp_mem_usage();
			vector<vector<int>> sample;
			sample.resize(n);
			vector<vector<int>> sample_reverse;
			sample_reverse.resize(n);
			gen_sample(sample, rumorset);
			sample.push_back(rumorset);
			tree.init(n, sample);
			tree.go(n);
			for (int i = 0; i < n; i++)
				sum[i] += (double)tree.sum[i];
		}
		double Max = 0;
		int remove_node = 0;
		for (int i = 0; i < n; i++)
			if (!remove_flag[i])
				if (sum[i] / ((double)count) > Max)
					Max = sum[i] / ((double)count), remove_node = i;
		return remove_node;
	}
	int DecreaseES_ON(int count, vector<int>rumorset, vector<bool>& remove_flag, vector<int> CB)
	{
		D.assign(n + 1, 0);
		vector<double> sum;
		sum.resize(n, 0);
		for (int j = 0; j < count; j++)//j代表了CP set的编号
		{
			//disp_mem_usage();
			vector<vector<int>> sample;
			sample.resize(n);
			vector<vector<int>> sample_reverse;
			sample_reverse.resize(n);
			gen_sample(sample, rumorset);
			sample.push_back(rumorset);
			tree.init(n, sample);
			tree.go(n);
			for (int i = 0; i < n; i++)
				sum[i] += (double)tree.sum[i];
		}
		double Max = 0;
		int remove_node = 0;
		for (int i : CB)
		{
			if (!remove_flag[i])
				if (sum[i] / ((double)count) > Max)
					Max = sum[i] / ((double)count), remove_node = i;
		}
		return remove_node;
	}
};