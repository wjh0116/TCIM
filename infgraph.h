#include "iheap.h"
#include <queue>	//priority_queue
#include <utility>  // pair
#include <numeric>
#include <algorithm>
#include <random>
#include <chrono>
#include <ctime>
using namespace std::chrono;
class InfGraph: public Graph
{
private:
	//vector<bool> activated;
    vector<bool> visit;
    vector<int> visit_mark;
	// int cound_fun;
	// double total_time_ms;
	// double total_time_WS;
public:
	// TRW set_1
	vector<vector<std::pair<int,double> > > node_Lhyperedges_1;
	vector<vector<std::pair<int,double> > > Lhyperedges_1;
	vector<double> node_Lhyperedge_weight_1;
	// TRW set_2
	vector<vector<std::pair<int,double> > > node_Lhyperedges_2;
	vector<vector<std::pair<int,double> > > Lhyperedges_2;
	vector<double> node_Lhyperedge_weight_2;
	// TRW path_1
	vector<vector<std::pair<int,double> > > node_Uhyperedges_1;
	vector<vector<std::pair<int,double> > > Uhyperedges_1;
	vector<double> node_Uhyperedge_weight_1;
	// TRW path_2
	vector<vector<std::pair<int,double> > > node_Uhyperedges_2;
	vector<vector<std::pair<int,double> > > Uhyperedges_2;
	vector<double> node_Uhyperedge_weight_2;
	vector<int> precomputedPoisson;


	

	sfmt_t sfmtSeed;

	

    InfGraph(string folder, string graph_file): Graph(folder, graph_file)
    {
        srand(time(NULL));
        sfmt_init_gen_rand(&sfmtSeed , rand());		
        visit = vector<bool> (n+1);
        visit_mark = vector<int> (n+1);
		isRumor = vector<bool>(n+1);
		// cound_fun=0;
		node_Lhyperedges_1=vector<vector<std::pair<int,double>>>(n);
		node_Lhyperedge_weight_1=vector<double>(n);
		node_Lhyperedges_2=vector<vector<std::pair<int,double>>>(n);
		node_Lhyperedge_weight_2=vector<double>(n);

		node_Uhyperedges_1=vector<vector<std::pair<int,double>>>(n);
		node_Uhyperedge_weight_1=vector<double>(n);
		node_Uhyperedges_2=vector<vector<std::pair<int,double>>>(n);
		node_Uhyperedge_weight_2=vector<double>(n);
    }

    void init_hyper_graph_seed(){
		seedSet.clear();
		UB_seedSet.clear();
		LB_seedSet.clear();
    }

	void precomputePoisson(double lambda, int sampleSize) {
		static thread_local std::mt19937 gen(std::random_device{}());
		std::poisson_distribution<> poisson(lambda);

		precomputedPoisson.reserve(sampleSize);
		for (int i = 0; i < sampleSize; ++i) {
			precomputedPoisson.push_back(poisson(gen) + 1);
		}
	}

	// int got_cnt()
	// {
	// 	return cound_fun;
	// }
	// double got_time()
	// {
	// 	return total_time_ms;
	// }
	// double got_ws_time()
	// {
	// 	return total_time_WS;
	// }

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

	/// @brief DT-based algorithms
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
			
		}
		void construct_CPset(int n, int source_node, vector<int>& CPSET, vector<bool>isRumor)
		{
			int current_v=source_node;
			while(isRumor[idom[current_v]]==0){
				current_v = idom[current_v];
				CPSET.push_back(current_v);
			}
		}
		
	} tree;

	/// Utility funtion
	int U(int x) {
		if(x<0) return 0;
		return x + 1;
	}

	/// Mente-Carlo simulations to estimate the misinformation's TD influence
	double MC_based_estimate(int count, int deadline)
	{
		int cnt=0;
		double sum=0;
		while (cnt<count)
		{
			priority_queue<Node> Q;
			vector<int> AT(n, 10000); 
			vector<int> I(n, 0); 
			for (int u : rumorSet) {
				AT[u] = 0;
				Q.push({u, 0});
			}
			while (!Q.empty()) 
			{
				Node current = Q.top();
				Q.pop();
				int q = current.id;
				int q_AT = current.activate_time;
				if (I[q] == 1) continue;

				for(int i=0;i<gT[q].size();i++)
				{
					int tmp= sfmt_genrand_uint32(&sfmtSeed) % precomputedPoisson.size();
					int diff_delay=precomputedPoisson[tmp]; 
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
				sum+=U(deadline-AT[q]);
			}
			cnt++;
		}
		return sum/(double)count;
	}

	double MC_based_estimate_with_removal(int count, int deadline, const vector<int>& removedNodes) {
		int cnt = 0;
		double sum = 0;
		vector<int> isBlocker(n, 0);
		for (int node : removedNodes) {
			isBlocker[node] = 1;
		}

		while (cnt < count) {
			priority_queue<Node> Q;
			vector<int> AT(n, 10000);
			vector<int> I(n, 0);

			for (int u : rumorSet) {
				AT[u] = 0;
				Q.push({u, 0});
			}

			while (!Q.empty()) {
				Node current = Q.top();
				Q.pop();

				int q = current.id;
				int q_AT = current.activate_time;

				if (I[q] == 1 || isBlocker[q]) continue; 

				for (int i = 0; i < gT[q].size(); i++) {
					int tmp= sfmt_genrand_uint32(&sfmtSeed) % precomputedPoisson.size();
					int diff_delay=precomputedPoisson[tmp]; 
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
				sum+=U(deadline-AT[q]);
			}
			cnt++;
		}

		return sum / (double)count;
	}



	/// sample a diffusion delay for each edge
	int sampleShiftedPoisson(double lambda) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::poisson_distribution<> poisson(lambda);

		return poisson(gen) + 1; 
	}

	/// @brief compute weight for TRW set, via the Dijkstra's algorithm
	double compute_AT(vector<vector<int>>& sample_graph, vector<vector<int>>& edge_length, int blockedNode, int target) {

    	vector<double> dist(n, 1000000); 
		priority_queue<pair<double, int>, vector<pair<double, int>>, greater<>> pq; 

		for (int seed : rumorSet) {
			dist[seed] = 0.0;
			pq.emplace(0.0, seed);
		}

		while (!pq.empty()) {
			auto [d, current] = pq.top();
			pq.pop();

			if (current == target) return d; // 如果到达目标节点，返回最短距离
			if (d > dist[current]) continue; // 如果当前距离大于已知最短距离，跳过

			// 遍历当前节点的所有出边
			for (int i = 0; i < sample_graph[current].size(); i++) {
				int next = sample_graph[current][i];
				double length = edge_length[current][i];

				// 跳过阻塞节点
				if (next == blockedNode) continue;

				// 更新距离
				if (dist[current] + length < dist[next]) {
					dist[next] = dist[current] + length;
					pq.emplace(dist[next], next);
				}
			}
		}
		return -1;// cannot activate due to the blocking
	}

	
	struct Node {// min-heap
    int id;
    int activate_time;
    bool operator<(const Node &other) const {
        return activate_time > other.activate_time; // Min-heap
    }
	};

	/// @brief generating a random TRW set and a random TRW path
	void WS_TRW_R1(Argument& arg)
	{
		priority_queue<Node> Q;
		vector<int> V; 
		vector<vector<int>>sample_graph(n);
		vector<vector<int>>edge_length(n);
		vector<int> AT(n, 10000); 
		vector<vector<int>> pred(n+1); 
		vector<int> Flag(n, 0); 
		vector<int> I(n, 0); 
		int deadline = arg.T;
		
		// auto start_ws = high_resolution_clock::now();
		for (int u : rumorSet) {
			AT[u] = 0;
			Q.push({u, 0});
    	}

		while (!Q.empty()) 
		{
			Node current = Q.top();
			Q.pop();
			int q = current.id;
			int q_AT = current.activate_time;
			if (I[q] == 1) continue;

			for(int i=0;i<gT[q].size();i++)
			{
				int tmp= sfmt_genrand_uint32(&sfmtSeed) % precomputedPoisson.size();
				int diff_delay=precomputedPoisson[tmp];  
				// int diff_delay=sampleShiftedPoisson(2); // 
				// cout<<"diff_delay: "<<diff_delay<<endl;
				int u=gT[q][i];
				if(sfmt_genrand_real1(&sfmtSeed) <= probT[q][i] && AT[q] + diff_delay <= deadline)
				{
					sample_graph[q].push_back(u);
					edge_length[q].push_back(diff_delay);
					if (AT[q] + diff_delay < AT[u]) {
						AT[u] = AT[q] + diff_delay;
						pred[u].clear();
						pred[u].push_back(q);
						Q.push({u, AT[u]});
					} else if (AT[q] + diff_delay == AT[u]) {
						pred[u].push_back(q);
					}
				}
			}	

			I[q] = 1;
			if(!isRumor[q]) V.push_back(q);

			if (!pred[q].empty()) {
				Flag[q] = Flag[pred[q][0]] + pred[q].size() - 1;
        	}
    	}
		// auto end_ws = high_resolution_clock::now();
		// duration<double> elapsed_ws= end_ws - start_ws;
		// total_time_WS += elapsed_ws.count();

		const auto random_num = sfmt_genrand_uint32(&sfmtSeed) % V.size(); 
		int source_node = V[random_num];

		// TRW set construction
		vector<std::pair<int,double>> Lhyperedge;
		int L_node=source_node;
		int T=arg.T;
		Lhyperedge.push_back(make_pair(source_node, U(T-AT[source_node]))); // source node

		if(Flag[L_node]==0)
		{
			// cout<<"1"<<endl;
			while(!isRumor[pred[L_node][0]])// need attention
			{
				L_node=pred[L_node][0];

				// auto start = high_resolution_clock::now();
				double aft_AT = compute_AT(sample_graph, edge_length, L_node, source_node);
				// cout<<"aft_AT: "<<aft_AT<<endl;
				// auto end = high_resolution_clock::now();
				// duration<double> elapsed= end - start;
				// cout<<"one time: "<<elapsed.count()<<endl;
    			// total_time_ms += elapsed.count();
				// cound_fun++;
				double weight;
				if(aft_AT==-1)
					weight=U(T-AT[source_node]);
				else
					weight=U(T-AT[source_node])-U(T-aft_AT);
				// cout<<"weight: "<<weight<<endl;
				Lhyperedge.push_back(make_pair(L_node, weight));
			}
		}
		else
		{
			for(int node:rumorSet)pred[node].push_back(n);
			vector<int>CP_set;
			//reverse pred
			vector<vector<int>> pred_reverse(n+1);
			for (int i = 0; i <= n; i++)
			{
				for (int v : pred[i])
					pred_reverse[v].push_back(i);
			}

			tree.init(n, pred_reverse);
			tree.go(n);
			tree.construct_CPset(n, source_node, CP_set,isRumor);
			for(int node:CP_set)
			{
				double aft_AT = compute_AT(sample_graph, edge_length, node, source_node);
				// cound_fun++;
				double weight;
				if(aft_AT==-1)
					weight=U(T-AT[source_node]);
				else
					weight=U(T-AT[source_node])-U(T-aft_AT);
				// cout<<"weight: "<<weight<<endl;
				Lhyperedge.push_back(make_pair(node, weight));
			}
		}
		
		// TRW path construction
		vector<std::pair<int,double>> Uhyperedge;
		int U_node=source_node;
		while(!isRumor[U_node])
		{
			Uhyperedge.push_back(make_pair(U_node, U(T - AT[source_node])));
			int nextV = numeric_limits<int>::max();
			for (int u : pred[U_node]) {
				nextV = min(nextV, u); 
			}
			U_node = nextV;
		}
	
		// add edges to R1
		Lhyperedges_1.push_back(Lhyperedge);
		Uhyperedges_1.push_back(Uhyperedge);

    	unsigned int index = Lhyperedges_1.size() - 1;

		for (unsigned int i = 0; i < Lhyperedge.size(); i++) {
			node_Lhyperedges_1[Lhyperedge[i].first].push_back(make_pair(index, Lhyperedge[i].second));
			node_Lhyperedge_weight_1[Lhyperedge[i].first] += Lhyperedge[i].second;	
	    }
		for (unsigned int i = 0; i < Uhyperedge.size(); i++) {
			node_Uhyperedges_1[Uhyperedge[i].first].push_back(make_pair(index, Uhyperedge[i].second));
			node_Uhyperedge_weight_1[Uhyperedge[i].first] += Uhyperedge[i].second;
		}

	}


	void WS_TRW_R2(Argument& arg)
	{
		priority_queue<Node> Q;
		vector<int> V; 
		vector<vector<int>>sample_graph(n);
		vector<vector<int>>edge_length(n);
		vector<int> AT(n, 10000); 
		vector<vector<int>> pred(n+1); 
		vector<int> Flag(n, 0); 
		vector<int> I(n, 0); 
		int deadline = arg.T;

		for (int u : rumorSet) {
			AT[u] = 0;
			Q.push({u, 0});
    	}

		while (!Q.empty()) 
		{
			Node current = Q.top();
			Q.pop();
			int q = current.id;
			int q_AT = current.activate_time;
			if (I[q] == 1) continue;

			for(int i=0;i<gT[q].size();i++)
			{
				int tmp= sfmt_genrand_uint32(&sfmtSeed) % precomputedPoisson.size();
				int diff_delay=precomputedPoisson[tmp];
				int u=gT[q][i];
				if(sfmt_genrand_real1(&sfmtSeed) <= probT[q][i] && AT[q] + diff_delay <= deadline)
				{
					sample_graph[q].push_back(u);
					edge_length[q].push_back(diff_delay);
					if (AT[q] + diff_delay < AT[u]) {
						AT[u] = AT[q] + diff_delay;
						pred[u].clear();
						pred[u].push_back(q);
						Q.push({u, AT[u]});
					} else if (AT[q] + diff_delay == AT[u]) {
						pred[u].push_back(q);
					}
				}
			}	

			I[q] = 1;
			if(!isRumor[q]) V.push_back(q);

			if (!pred[q].empty()) {
				Flag[q] = Flag[pred[q][0]] + pred[q].size() - 1;
        	}
    	}
		const auto random_num = sfmt_genrand_uint32(&sfmtSeed) % V.size(); 
		int source_node = V[random_num];

		// TRW set construction
		vector<std::pair<int,double>> Lhyperedge;
		int L_node=source_node;
		int T=arg.T;
		Lhyperedge.push_back(make_pair(source_node, U(T-AT[source_node]))); // source node

		if(Flag[L_node]==0)
		{
			// cout<<"1"<<endl;
			while(!isRumor[pred[L_node][0]])// need attention
			{
				L_node=pred[L_node][0];
				double aft_AT = compute_AT(sample_graph, edge_length, L_node, source_node);
				// cound_fun++;
				double weight;
				if(aft_AT==-1)
					weight=U(T-AT[source_node]);
				else
					weight=U(T-AT[source_node])-U(T-aft_AT);
				// cout<<"weight: "<<weight<<endl;
				Lhyperedge.push_back(make_pair(L_node, weight));
			}
		}
		else
		{
			for(int node:rumorSet)pred[node].push_back(n);
			vector<int>CP_set;
			//reverse pred
			vector<vector<int>> pred_reverse(n+1);
			for (int i = 0; i <= n; i++)
			{
				for (int v : pred[i])
					pred_reverse[v].push_back(i);
			}

			tree.init(n, pred_reverse);
			tree.go(n);
			tree.construct_CPset(n, source_node, CP_set,isRumor);
			for(int node:CP_set)
			{
				double aft_AT = compute_AT(sample_graph, edge_length, node, source_node);
				// cound_fun++;
				double weight;
				if(aft_AT==-1)
					weight=U(T-AT[source_node]);
				else
					weight=U(T-AT[source_node])-U(T-aft_AT);
				// cout<<"weight: "<<weight<<endl;
				Lhyperedge.push_back(make_pair(node, weight));
			}
		}
		
		// TRW path construction
		vector<std::pair<int,double>> Uhyperedge;
		int U_node=source_node;
		while(!isRumor[U_node])
		{
			Uhyperedge.push_back(make_pair(U_node, U(T - AT[source_node])));
			int nextV = numeric_limits<int>::max();
			for (int u : pred[U_node]) {
				nextV = min(nextV, u); 
			}
			U_node = nextV;
		}
	
		// add edges to R2
		Lhyperedges_2.push_back(Lhyperedge);
		Uhyperedges_2.push_back(Uhyperedge);

    	unsigned int index = Lhyperedges_2.size() - 1;

		for (unsigned int i = 0; i < Lhyperedge.size(); i++) {
			node_Lhyperedges_2[Lhyperedge[i].first].push_back(make_pair(index, Lhyperedge[i].second));
			node_Lhyperedge_weight_2[Lhyperedge[i].first] += Lhyperedge[i].second;	
	    }
		for (unsigned int i = 0; i < Uhyperedge.size(); i++) {
			node_Uhyperedges_2[Uhyperedge[i].first].push_back(make_pair(index, Uhyperedge[i].second));
			node_Uhyperedge_weight_2[Uhyperedge[i].first] += Uhyperedge[i].second;
		}

	}

	void generate_rr_set(int source_node, vector<vector<int>>&sample_graph, vector<int>&reachable_nodes)
	{
		vector<int> vis;
		vis.resize(n,0);
	
			
		queue<int> q;
		q.push(source_node);
		vis[source_node]=1;
		reachable_nodes.push_back(source_node);

		while (!q.empty())
		{
			int x = q.front();
			q.pop();
			for (int i = 0; i < sample_graph[x].size(); i++)
			{
				if (vis[sample_graph[x][i]] != 1)
					q.push(sample_graph[x][i]), vis[sample_graph[x][i]] = 1, reachable_nodes.push_back(sample_graph[x][i]);
			}
		}
	}

	void WS_TRW_wo_PR_R1(Argument& arg)
	{
		priority_queue<Node> Q;
		vector<int> V; 
		vector<vector<int>>sample_graph(n);
		vector<vector<int>>edge_length(n);
		vector<int> AT(n, 10000); 
		vector<vector<int>> pred(n+1); 
		vector<int> Flag(n, 0); 
		vector<int> I(n, 0); 
		int deadline = arg.T;

		for (int u : rumorSet) {
			AT[u] = 0;
			Q.push({u, 0});
    	}

		while (!Q.empty()) 
		{
			Node current = Q.top();
			Q.pop();
			int q = current.id;
			int q_AT = current.activate_time;
			if (I[q] == 1) continue;

			for(int i=0;i<gT[q].size();i++)
			{
				int tmp= sfmt_genrand_uint32(&sfmtSeed) % precomputedPoisson.size();
				int diff_delay=precomputedPoisson[tmp]; 
				// int diff_delay=sampleShiftedPoisson(2); // 
				// cout<<"diff_delay: "<<diff_delay<<endl;
				int u=gT[q][i];
				if(sfmt_genrand_real1(&sfmtSeed) <= probT[q][i] && AT[q] + diff_delay <= deadline)
				{
					sample_graph[q].push_back(u);
					edge_length[q].push_back(diff_delay);
					if (AT[q] + diff_delay < AT[u]) {
						AT[u] = AT[q] + diff_delay;
						pred[u].clear();
						pred[u].push_back(q);
						Q.push({u, AT[u]});
					} else if (AT[q] + diff_delay == AT[u]) {
						pred[u].push_back(q);
					}
				}
			}	

			I[q] = 1;
			if(!isRumor[q]) V.push_back(q);

			if (!pred[q].empty()) {
				Flag[q] = Flag[pred[q][0]] + pred[q].size() - 1;
        	}
    	}
		const auto random_num = sfmt_genrand_uint32(&sfmtSeed) % V.size(); 
		int source_node = V[random_num];

		// TRW set construction
		vector<std::pair<int,double>> Lhyperedge;
		int L_node=source_node;
		int T=arg.T;
		// Lhyperedge.push_back(make_pair(source_node, U(T-AT[source_node]))); // source node

		//determine V(R)
		vector<int>reachable_nodes;
		vector<vector<int>> sample_graph_reverse(n);
		for (int i = 0; i < n; i++)
		{
			for (int v : sample_graph[i])
				sample_graph_reverse[v].push_back(i);
		}
		generate_rr_set(source_node, sample_graph_reverse, reachable_nodes);
		// cout<<"reachable_nodes.size(): "<<reachable_nodes.size()<<endl;
		for(int node:reachable_nodes)
		{
			// auto start = high_resolution_clock::now();
			double aft_AT = compute_AT(sample_graph, edge_length, node, source_node);
			// cout<<"aft_AT1: "<<aft_AT<<endl;
			// auto end = high_resolution_clock::now();
			// duration<double> elapsed= end - start;
			// total_time_ms += elapsed.count();
		
			// cound_fun++;
			double weight;
			if(aft_AT==-1)
				weight=U(T-AT[source_node]);
			else
				weight=U(T-AT[source_node])-U(T-aft_AT);
			// cout<<"weight: "<<weight<<endl;
			Lhyperedge.push_back(make_pair(node, weight));
		}

	
		// TRW path construction
		vector<std::pair<int,double>> Uhyperedge;
		int U_node=source_node;
		while(!isRumor[U_node])
		{
			Uhyperedge.push_back(make_pair(U_node, U(T - AT[source_node])));
			int nextV = numeric_limits<int>::max();
			for (int u : pred[U_node]) {
				nextV = min(nextV, u); 
			}
			U_node = nextV;
		}
	
		// add edges to R1
		Lhyperedges_1.push_back(Lhyperedge);
		Uhyperedges_1.push_back(Uhyperedge);

    	unsigned int index = Lhyperedges_1.size() - 1;

		for (unsigned int i = 0; i < Lhyperedge.size(); i++) {
			node_Lhyperedges_1[Lhyperedge[i].first].push_back(make_pair(index, Lhyperedge[i].second));
			node_Lhyperedge_weight_1[Lhyperedge[i].first] += Lhyperedge[i].second;	
	    }
		for (unsigned int i = 0; i < Uhyperedge.size(); i++) {
			node_Uhyperedges_1[Uhyperedge[i].first].push_back(make_pair(index, Uhyperedge[i].second));
			node_Uhyperedge_weight_1[Uhyperedge[i].first] += Uhyperedge[i].second;
		}

	}

	void WS_TRW_wo_PR_R2(Argument& arg)
	{
		priority_queue<Node> Q;
		vector<int> V; 
		vector<vector<int>>sample_graph(n);
		vector<vector<int>>edge_length(n);
		vector<int> AT(n, 10000); 
		vector<vector<int>> pred(n+1); 
		vector<int> Flag(n, 0); 
		vector<int> I(n, 0); 
		int deadline = arg.T;

		for (int u : rumorSet) {
			AT[u] = 0;
			Q.push({u, 0});
    	}

		while (!Q.empty()) 
		{
			Node current = Q.top();
			Q.pop();
			int q = current.id;
			int q_AT = current.activate_time;
			if (I[q] == 1) continue;

			for(int i=0;i<gT[q].size();i++)
			{
				int tmp= sfmt_genrand_uint32(&sfmtSeed) % precomputedPoisson.size();
				int diff_delay=precomputedPoisson[tmp];
				int u=gT[q][i];
				if(sfmt_genrand_real1(&sfmtSeed) <= probT[q][i] && AT[q] + diff_delay <= deadline)
				{
					sample_graph[q].push_back(u);
					edge_length[q].push_back(diff_delay);
					if (AT[q] + diff_delay < AT[u]) {
						AT[u] = AT[q] + diff_delay;
						pred[u].clear();
						pred[u].push_back(q);
						Q.push({u, AT[u]});
					} else if (AT[q] + diff_delay == AT[u]) {
						pred[u].push_back(q);
					}
				}
			}	

			I[q] = 1;
			if(!isRumor[q]) V.push_back(q);

			if (!pred[q].empty()) {
				Flag[q] = Flag[pred[q][0]] + pred[q].size() - 1;
        	}
    	}
		const auto random_num = sfmt_genrand_uint32(&sfmtSeed) % V.size(); 
		int source_node = V[random_num];

		// TRW set construction
		vector<std::pair<int,double>> Lhyperedge;
		int L_node=source_node;
		int T=arg.T;
		// Lhyperedge.push_back(make_pair(source_node, U(T-AT[source_node]))); // source node

		//determine V(R)
		vector<int>reachable_nodes;
		vector<vector<int>> sample_graph_reverse(n);
		for (int i = 0; i < n; i++)
		{
			for (int v : sample_graph[i])
				sample_graph_reverse[v].push_back(i);
		}
		generate_rr_set(source_node, sample_graph_reverse, reachable_nodes);
		// cout<<"reachable_nodes.size(): "<<reachable_nodes.size()<<endl;
		for(int node:reachable_nodes)
		{
			double aft_AT = compute_AT(sample_graph, edge_length, node, source_node);
			// cound_fun++;
			double weight;
			if(aft_AT==-1)
				weight=U(T-AT[source_node]);
			else
				weight=U(T-AT[source_node])-U(T-aft_AT);
			// cout<<"weight: "<<weight<<endl;
			Lhyperedge.push_back(make_pair(node, weight));
		}

	
		// TRW path construction
		vector<std::pair<int,double>> Uhyperedge;
		int U_node=source_node;
		while(!isRumor[U_node])
		{
			Uhyperedge.push_back(make_pair(U_node, U(T - AT[source_node])));
			int nextV = numeric_limits<int>::max();
			for (int u : pred[U_node]) {
				nextV = min(nextV, u); 
			}
			U_node = nextV;
		}
	
		// add edges to R2
		Lhyperedges_2.push_back(Lhyperedge);
		Uhyperedges_2.push_back(Uhyperedge);

    	unsigned int index = Lhyperedges_2.size() - 1;

		for (unsigned int i = 0; i < Lhyperedge.size(); i++) {
			node_Lhyperedges_2[Lhyperedge[i].first].push_back(make_pair(index, Lhyperedge[i].second));
			node_Lhyperedge_weight_2[Lhyperedge[i].first] += Lhyperedge[i].second;	
	    }
		for (unsigned int i = 0; i < Uhyperedge.size(); i++) {
			node_Uhyperedges_2[Uhyperedge[i].first].push_back(make_pair(index, Uhyperedge[i].second));
			node_Uhyperedge_weight_2[Uhyperedge[i].first] += Uhyperedge[i].second;
		}

	}



	void US_TRW_R1(Argument& arg)
	{
		priority_queue<Node> Q;
		vector<int> V; 
		vector<vector<int>>sample_graph(n);
		vector<vector<int>>edge_length(n);
		vector<int> AT(n, 10000); 
		vector<vector<int>> pred(n+1); 
		vector<int> Flag(n, 0); 
		vector<int> I(n, 0); 
		int deadline = arg.T;

		int source_node;
		int f1=1;

		while(f1)
		{
			
			//reset
			V.clear();
			for(int i=0;i<n;i++)
			{
				sample_graph[i].clear();
				edge_length[i].clear();
				I[i]=0;
				Flag[i]=0;
				pred[i].clear();
				AT[i]=10000;
			}
			pred[n].clear();

			do {
				source_node = sfmt_genrand_uint32(&sfmtSeed) % n; 
			} while (isRumor[source_node]==1);

			for (int u : rumorSet) {
				AT[u] = 0;
				Q.push({u, 0});
			}
			while (!Q.empty()) 
			{
				Node current = Q.top();
				Q.pop();
				int q = current.id;
				int q_AT = current.activate_time;
				if (I[q] == 1) continue;

				for(int i=0;i<gT[q].size();i++)
				{
					int tmp= sfmt_genrand_uint32(&sfmtSeed) % precomputedPoisson.size();
					int diff_delay=precomputedPoisson[tmp]; 
					int u=gT[q][i];
					if(sfmt_genrand_real1(&sfmtSeed) <= probT[q][i] && AT[q] + diff_delay <= deadline)
					{
						sample_graph[q].push_back(u);
						edge_length[q].push_back(diff_delay);
						if (AT[q] + diff_delay < AT[u]) {
							AT[u] = AT[q] + diff_delay;
							pred[u].clear();
							pred[u].push_back(q);
							Q.push({u, AT[u]});
						} else if (AT[q] + diff_delay == AT[u]) {
							pred[u].push_back(q);
						}
					}
				}	

				I[q] = 1;
				if(q==source_node) f1=0;
				if(!isRumor[q]) V.push_back(q);

				if (!pred[q].empty()) {
					Flag[q] = Flag[pred[q][0]] + pred[q].size() - 1;
				}
			}
		}

		// cout<<"source_node: "<<source_node<<endl;
		
		// const auto random_num = sfmt_genrand_uint32(&sfmtSeed) % V.size(); 
		// int source_node = V[random_num];

		// TRW set construction
		vector<std::pair<int,double>> Lhyperedge;
		int L_node=source_node;
		int T=arg.T;
		Lhyperedge.push_back(make_pair(source_node, U(T-AT[source_node]))); // source node

		if(Flag[L_node]==0)
		{
			// cout<<"1"<<endl;
			while(!isRumor[pred[L_node][0]])// need attention
			{
				L_node=pred[L_node][0];
				double aft_AT = compute_AT(sample_graph, edge_length, L_node, source_node);
				double weight;
				if(aft_AT==-1)
					weight=U(T-AT[source_node]);
				else
					weight=U(T-AT[source_node])-U(T-aft_AT);
				// cout<<"weight: "<<weight<<endl;
				Lhyperedge.push_back(make_pair(L_node, weight));
			}
		}
		else
		{
			for(int node:rumorSet)pred[node].push_back(n);
			vector<int>CP_set;
			//reverse pred
			vector<vector<int>> pred_reverse(n+1);
			for (int i = 0; i <= n; i++)
			{
				for (int v : pred[i])
					pred_reverse[v].push_back(i);
			}

			tree.init(n, pred_reverse);
			tree.go(n);
			tree.construct_CPset(n, source_node, CP_set,isRumor);
			for(int node:CP_set)
			{
				double aft_AT = compute_AT(sample_graph, edge_length, node, source_node);
				double weight;
				if(aft_AT==-1)
					weight=U(T-AT[source_node]);
				else
					weight=U(T-AT[source_node])-U(T-aft_AT);
				// cout<<"weight: "<<weight<<endl;
				Lhyperedge.push_back(make_pair(node, weight));
			}
		}
		
		// TRW path construction
		vector<std::pair<int,double>> Uhyperedge;
		int U_node=source_node;
		while(!isRumor[U_node])
		{
			Uhyperedge.push_back(make_pair(U_node, U(T - AT[source_node])));
			int nextV = numeric_limits<int>::max();
			for (int u : pred[U_node]) {
				nextV = min(nextV, u); 
			}
			U_node = nextV;
		}
	
		// add edges to R1
		Lhyperedges_1.push_back(Lhyperedge);
		Uhyperedges_1.push_back(Uhyperedge);

    	unsigned int index = Lhyperedges_1.size() - 1;

		for (unsigned int i = 0; i < Lhyperedge.size(); i++) {
			node_Lhyperedges_1[Lhyperedge[i].first].push_back(make_pair(index, Lhyperedge[i].second));
			node_Lhyperedge_weight_1[Lhyperedge[i].first] += Lhyperedge[i].second;	
	    }
		for (unsigned int i = 0; i < Uhyperedge.size(); i++) {
			node_Uhyperedges_1[Uhyperedge[i].first].push_back(make_pair(index, Uhyperedge[i].second));
			node_Uhyperedge_weight_1[Uhyperedge[i].first] += Uhyperedge[i].second;
		}

	}

	void US_TRW_R2(Argument& arg)
	{
		priority_queue<Node> Q;
		vector<int> V; 
		vector<vector<int>>sample_graph(n);
		vector<vector<int>>edge_length(n);
		vector<int> AT(n, 10000); 
		vector<vector<int>> pred(n+1); 
		vector<int> Flag(n, 0); 
		vector<int> I(n, 0); 
		int deadline = arg.T;

		int source_node;
		int f1=1;

		while(f1)
		{
			
			//reset
			V.clear();
			for(int i=0;i<n;i++)
			{
				sample_graph[i].clear();
				edge_length[i].clear();
				I[i]=0;
				Flag[i]=0;
				pred[i].clear();
				AT[i]=10000;
			}
			pred[n].clear();

			do {
				source_node = sfmt_genrand_uint32(&sfmtSeed) % n; 
			} while (isRumor[source_node]==1);

			for (int u : rumorSet) {
				AT[u] = 0;
				Q.push({u, 0});
			}
			while (!Q.empty()) 
			{
				Node current = Q.top();
				Q.pop();
				int q = current.id;
				int q_AT = current.activate_time;
				if (I[q] == 1) continue;

				for(int i=0;i<gT[q].size();i++)
				{
					int tmp= sfmt_genrand_uint32(&sfmtSeed) % precomputedPoisson.size();
					int diff_delay=precomputedPoisson[tmp];
					int u=gT[q][i];
					if(sfmt_genrand_real1(&sfmtSeed) <= probT[q][i] && AT[q] + diff_delay <= deadline)
					{
						sample_graph[q].push_back(u);
						edge_length[q].push_back(diff_delay);
						if (AT[q] + diff_delay < AT[u]) {
							AT[u] = AT[q] + diff_delay;
							pred[u].clear();
							pred[u].push_back(q);
							Q.push({u, AT[u]});
						} else if (AT[q] + diff_delay == AT[u]) {
							pred[u].push_back(q);
						}
					}
				}	

				I[q] = 1;
				if(q==source_node) f1=0;
				if(!isRumor[q]) V.push_back(q);

				if (!pred[q].empty()) {
					Flag[q] = Flag[pred[q][0]] + pred[q].size() - 1;
				}
			}
		}

		// cout<<"source_node: "<<source_node<<endl;
		
		// const auto random_num = sfmt_genrand_uint32(&sfmtSeed) % V.size(); 
		// int source_node = V[random_num];

		// TRW set construction
		vector<std::pair<int,double>> Lhyperedge;
		int L_node=source_node;
		int T=arg.T;
		Lhyperedge.push_back(make_pair(source_node, U(T-AT[source_node]))); // source node

		if(Flag[L_node]==0)
		{
			// cout<<"1"<<endl;
			while(!isRumor[pred[L_node][0]])// need attention
			{
				L_node=pred[L_node][0];
				double aft_AT = compute_AT(sample_graph, edge_length, L_node, source_node);
				double weight;
				if(aft_AT==-1)
					weight=U(T-AT[source_node]);
				else
					weight=U(T-AT[source_node])-U(T-aft_AT);
				// cout<<"weight: "<<weight<<endl;
				Lhyperedge.push_back(make_pair(L_node, weight));
			}
		}
		else
		{
			for(int node:rumorSet)pred[node].push_back(n);
			vector<int>CP_set;
			//reverse pred
			vector<vector<int>> pred_reverse(n+1);
			for (int i = 0; i <= n; i++)
			{
				for (int v : pred[i])
					pred_reverse[v].push_back(i);
			}

			tree.init(n, pred_reverse);
			tree.go(n);
			tree.construct_CPset(n, source_node, CP_set,isRumor);
			for(int node:CP_set)
			{
				double aft_AT = compute_AT(sample_graph, edge_length, node, source_node);
				double weight;
				if(aft_AT==-1)
					weight=U(T-AT[source_node]);
				else
					weight=U(T-AT[source_node])-U(T-aft_AT);
				// cout<<"weight: "<<weight<<endl;
				Lhyperedge.push_back(make_pair(node, weight));
			}
		}
		
		// TRW path construction
		vector<std::pair<int,double>> Uhyperedge;
		int U_node=source_node;
		while(!isRumor[U_node])
		{
			Uhyperedge.push_back(make_pair(U_node, U(T - AT[source_node])));
			int nextV = numeric_limits<int>::max();
			for (int u : pred[U_node]) {
				nextV = min(nextV, u); 
			}
			U_node = nextV;
		}
	
		// add edges to R2
		Lhyperedges_2.push_back(Lhyperedge);
		Uhyperedges_2.push_back(Uhyperedge);

    	unsigned int index = Lhyperedges_2.size() - 1;

		for (unsigned int i = 0; i < Lhyperedge.size(); i++) {
			node_Lhyperedges_2[Lhyperedge[i].first].push_back(make_pair(index, Lhyperedge[i].second));
			node_Lhyperedge_weight_2[Lhyperedge[i].first] += Lhyperedge[i].second;	
	    }
		for (unsigned int i = 0; i < Uhyperedge.size(); i++) {
			node_Uhyperedges_2[Uhyperedge[i].first].push_back(make_pair(index, Uhyperedge[i].second));
			node_Uhyperedge_weight_2[Uhyperedge[i].first] += Uhyperedge[i].second;
		}

	}

	void US_TRW_wo_PR_R1(Argument& arg)
	{
		priority_queue<Node> Q;
		vector<int> V; 
		vector<vector<int>>sample_graph(n);
		vector<vector<int>>edge_length(n);
		vector<int> AT(n, 10000); 
		vector<vector<int>> pred(n+1); 
		vector<int> Flag(n, 0); 
		vector<int> I(n, 0); 
		int deadline = arg.T;

		int source_node;
		int f1=1;

		while(f1)
		{
			
			//reset
			V.clear();
			for(int i=0;i<n;i++)
			{
				sample_graph[i].clear();
				edge_length[i].clear();
				I[i]=0;
				Flag[i]=0;
				pred[i].clear();
				AT[i]=10000;
			}
			pred[n].clear();

			do {
				source_node = sfmt_genrand_uint32(&sfmtSeed) % n; 
			} while (isRumor[source_node]==1);

			for (int u : rumorSet) {
				AT[u] = 0;
				Q.push({u, 0});
			}
			while (!Q.empty()) 
			{
				Node current = Q.top();
				Q.pop();
				int q = current.id;
				int q_AT = current.activate_time;
				if (I[q] == 1) continue;

				for(int i=0;i<gT[q].size();i++)
				{
					int tmp= sfmt_genrand_uint32(&sfmtSeed) % precomputedPoisson.size();
					int diff_delay=precomputedPoisson[tmp];
					int u=gT[q][i];
					if(sfmt_genrand_real1(&sfmtSeed) <= probT[q][i] && AT[q] + diff_delay <= deadline)
					{
						sample_graph[q].push_back(u);
						edge_length[q].push_back(diff_delay);
						if (AT[q] + diff_delay < AT[u]) {
							AT[u] = AT[q] + diff_delay;
							pred[u].clear();
							pred[u].push_back(q);
							Q.push({u, AT[u]});
						} else if (AT[q] + diff_delay == AT[u]) {
							pred[u].push_back(q);
						}
					}
				}	

				I[q] = 1;
				if(q==source_node) f1=0;
				if(!isRumor[q]) V.push_back(q);

				if (!pred[q].empty()) {
					Flag[q] = Flag[pred[q][0]] + pred[q].size() - 1;
				}
			}
		}

		// TRW set construction
		vector<std::pair<int,double>> Lhyperedge;
		int L_node=source_node;
		int T=arg.T;
	
		vector<int>reachable_nodes;
		vector<vector<int>> sample_graph_reverse(n);
		for (int i = 0; i < n; i++)
		{
			for (int v : sample_graph[i])
				sample_graph_reverse[v].push_back(i);
		}
		generate_rr_set(source_node, sample_graph_reverse, reachable_nodes);
		// cout<<"reachable_nodes.size(): "<<reachable_nodes.size()<<endl;
		for(int node:reachable_nodes)
		{
			double aft_AT = compute_AT(sample_graph, edge_length, node, source_node);
			double weight;
			if(aft_AT==-1)
				weight=U(T-AT[source_node]);
			else
				weight=U(T-AT[source_node])-U(T-aft_AT);
			// cout<<"weight: "<<weight<<endl;
			Lhyperedge.push_back(make_pair(node, weight));
		}
		
		// TRW path construction
		vector<std::pair<int,double>> Uhyperedge;
		int U_node=source_node;
		while(!isRumor[U_node])
		{
			Uhyperedge.push_back(make_pair(U_node, U(T - AT[source_node])));
			int nextV = numeric_limits<int>::max();
			for (int u : pred[U_node]) {
				nextV = min(nextV, u); 
			}
			U_node = nextV;
		}
	
		// add edges to R1
		Lhyperedges_1.push_back(Lhyperedge);
		Uhyperedges_1.push_back(Uhyperedge);

    	unsigned int index = Lhyperedges_1.size() - 1;

		for (unsigned int i = 0; i < Lhyperedge.size(); i++) {
			node_Lhyperedges_1[Lhyperedge[i].first].push_back(make_pair(index, Lhyperedge[i].second));
			node_Lhyperedge_weight_1[Lhyperedge[i].first] += Lhyperedge[i].second;	
	    }
		for (unsigned int i = 0; i < Uhyperedge.size(); i++) {
			node_Uhyperedges_1[Uhyperedge[i].first].push_back(make_pair(index, Uhyperedge[i].second));
			node_Uhyperedge_weight_1[Uhyperedge[i].first] += Uhyperedge[i].second;
		}

	}

	void US_TRW_wo_PR_R2(Argument& arg)
	{
		priority_queue<Node> Q;
		vector<int> V; 
		vector<vector<int>>sample_graph(n);
		vector<vector<int>>edge_length(n);
		vector<int> AT(n, 10000); 
		vector<vector<int>> pred(n+1); 
		vector<int> Flag(n, 0); 
		vector<int> I(n, 0); 
		int deadline = arg.T;

		int source_node;
		int f1=1;

		while(f1)
		{
			
			//reset
			V.clear();
			for(int i=0;i<n;i++)
			{
				sample_graph[i].clear();
				edge_length[i].clear();
				I[i]=0;
				Flag[i]=0;
				pred[i].clear();
				AT[i]=10000;
			}
			pred[n].clear();

			do {
				source_node = sfmt_genrand_uint32(&sfmtSeed) % n; 
			} while (isRumor[source_node]==1);

			for (int u : rumorSet) {
				AT[u] = 0;
				Q.push({u, 0});
			}
			while (!Q.empty()) 
			{
				Node current = Q.top();
				Q.pop();
				int q = current.id;
				int q_AT = current.activate_time;
				if (I[q] == 1) continue;

				for(int i=0;i<gT[q].size();i++)
				{
					int tmp= sfmt_genrand_uint32(&sfmtSeed) % precomputedPoisson.size();
					int diff_delay=precomputedPoisson[tmp];
					int u=gT[q][i];
					if(sfmt_genrand_real1(&sfmtSeed) <= probT[q][i] && AT[q] + diff_delay <= deadline)
					{
						sample_graph[q].push_back(u);
						edge_length[q].push_back(diff_delay);
						if (AT[q] + diff_delay < AT[u]) {
							AT[u] = AT[q] + diff_delay;
							pred[u].clear();
							pred[u].push_back(q);
							Q.push({u, AT[u]});
						} else if (AT[q] + diff_delay == AT[u]) {
							pred[u].push_back(q);
						}
					}
				}	

				I[q] = 1;
				if(q==source_node) f1=0;
				if(!isRumor[q]) V.push_back(q);

				if (!pred[q].empty()) {
					Flag[q] = Flag[pred[q][0]] + pred[q].size() - 1;
				}
			}
		}

		// TRW set construction
		vector<std::pair<int,double>> Lhyperedge;
		int L_node=source_node;
		int T=arg.T;
	
		vector<int>reachable_nodes;
		vector<vector<int>> sample_graph_reverse(n);
		for (int i = 0; i < n; i++)
		{
			for (int v : sample_graph[i])
				sample_graph_reverse[v].push_back(i);
		}
		generate_rr_set(source_node, sample_graph_reverse, reachable_nodes);
		// cout<<"reachable_nodes.size(): "<<reachable_nodes.size()<<endl;
		for(int node:reachable_nodes)
		{
			double aft_AT = compute_AT(sample_graph, edge_length, node, source_node);
			double weight;
			if(aft_AT==-1)
				weight=U(T-AT[source_node]);
			else
				weight=U(T-AT[source_node])-U(T-aft_AT);
			// cout<<"weight: "<<weight<<endl;
			Lhyperedge.push_back(make_pair(node, weight));
		}
		
		// TRW path construction
		vector<std::pair<int,double>> Uhyperedge;
		int U_node=source_node;
		while(!isRumor[U_node])
		{
			Uhyperedge.push_back(make_pair(U_node, U(T - AT[source_node])));
			int nextV = numeric_limits<int>::max();
			for (int u : pred[U_node]) {
				nextV = min(nextV, u); 
			}
			U_node = nextV;
		}
	
		// add edges to R2
		Lhyperedges_2.push_back(Lhyperedge);
		Uhyperedges_2.push_back(Uhyperedge);

    	unsigned int index = Lhyperedges_2.size() - 1;

		for (unsigned int i = 0; i < Lhyperedge.size(); i++) {
			node_Lhyperedges_2[Lhyperedge[i].first].push_back(make_pair(index, Lhyperedge[i].second));
			node_Lhyperedge_weight_2[Lhyperedge[i].first] += Lhyperedge[i].second;	
	    }
		for (unsigned int i = 0; i < Uhyperedge.size(); i++) {
			node_Uhyperedges_2[Uhyperedge[i].first].push_back(make_pair(index, Uhyperedge[i].second));
			node_Uhyperedge_weight_2[Uhyperedge[i].first] += Uhyperedge[i].second;
		}

	}

	/*
	* linear pass over coverages to find node with maximum marginal coverage
	* also maintains top k marginals for computation of improved upper bound
	*/
	int getMaxIndex(int n, vector<int> &node_weight, vector<int> &k_max_mc) {
		int max_ind = -1;
		int max_cov = 0;

		for (int i = 0; i < n; i++) {
			if (node_weight[i] > max_cov) {
				max_ind = i;
				max_cov = node_weight[i];
			}
			if (node_weight[i] > k_max_mc[0]) {
				k_max_mc[0] = node_weight[i];
				sort(k_max_mc.begin(), k_max_mc.end());
			}
		} 

		return max_ind;
	}

	int buildSeedSet(vector<vector<std::pair<int,double> > >& hyperedges,vector<vector<std::pair<int,double> > > &node_hyperedges,vector<double> &node_hyperedge_weight, 
	unsigned int k, vector<int>& seedSet)
	{	
		unsigned int i, j;
		int diff, coverage_weight, cur_cov_ub, max_index;
		bool alive;
		vector<pair<int,double> > edge_list, node_list;
		vector<int> k_max_mc(k,0);

		vector<int> node_weight(n,0);
		
		for (i = 0; i < n; i++) {
			node_weight[i] = node_hyperedge_weight[i];
		}
		
		int cur_coverage = 0;
		int improved_cov_ub = INT_MAX;
		long long numEdge = hyperedges.size();

		// check if an edge is removed
		vector<bool> edge_removed(numEdge, false);
		
		unsigned int cur_seed = 0;
		// building each seed at a time
		seedSet.clear();
		while(cur_seed < k) {
			
			max_index = getMaxIndex(n, node_weight, k_max_mc);
			if (max_index == -1) break; // all sets have been covered 

			cur_cov_ub = cur_coverage;
			for (i = 0; i < k; i++) {
				cur_cov_ub += k_max_mc[i];
				k_max_mc[i] = 0; // reset for next iteration
			}
			if (cur_cov_ub < improved_cov_ub) improved_cov_ub = cur_cov_ub;

			seedSet.push_back(max_index);
			cur_coverage += node_weight[max_index];

			edge_list = node_hyperedges[max_index];// the sample sets convered by 'max_index' 
			for (i = 0; i < edge_list.size(); i++) {
				if (edge_removed[edge_list[i].first]) continue;
				alive = false;
				coverage_weight = edge_list[i].second; //'max_index'  在某个sample set上 cover 的 weight
				node_list = hyperedges[edge_list[i].first];// 这个sample set中的<node, weight>
				for (j = 0; j < node_list.size(); j++) {
					if (node_list[j].second > 0) { // only process those nodes with a non-zero marginal for this hyperedge
						diff = node_list[j].second - coverage_weight;// 这个节点剩余对这个sample的marginal gain
						if (diff <= 0) {
							node_weight[node_list[j].first] -= node_list[j].second;
							node_list[j].second = 0;
						} else {
							node_weight[node_list[j].first] -= coverage_weight;
							node_list[j].second = diff;
						}
						if (!alive && node_list[j].second > 0) alive = true; // still marginal gain achievable from this hyperedge
					}
				}
				if (!alive) edge_removed[edge_list[i].first] = true;
			}
			cur_seed++;
		}

		getMaxIndex(n, node_weight, k_max_mc);
		cur_cov_ub = cur_coverage;
		for (i = 0; i < k; i++) {
			cur_cov_ub += k_max_mc[i];
		}
		if (cur_cov_ub < improved_cov_ub) improved_cov_ub = cur_cov_ub;

		return improved_cov_ub;
	}

	// 
	int getMaxIndexBaseline(int n, vector<int> &node_weight) {
		int max_ind = -1;
		int max_cov = 0;

		for (int i = 0; i < n; i++) {
			if (node_weight[i] > max_cov) {
				max_ind = i;
				max_cov = node_weight[i];
			}
		}

		return max_ind;
	}

	/*
	* greedy algorithm for weighted max cover over collection of RDR sets w/o improved UB computation
	*/
	int buildSeedSetBaseline(vector<vector<std::pair<int,double> > > hyperedges,vector<vector<std::pair<int,double> > > node_hyperedges,vector<double> node_hyperedge_weight, 
	unsigned int k, vector<int>& seedSet)
	{	
		unsigned int i, j;
		int diff, coverage_weight, max_index;
		bool alive;
		vector<pair<int,double> > edge_list, node_list;

		vector<int> node_weight(n,0);
		for (i = 0; i < n; i++) {
			node_weight[i] = node_hyperedge_weight[i];
		}

		int coverage = 0;
		long long numEdge = hyperedges.size();

		// check if an edge is removed
		vector<bool> edge_removed(numEdge, false);
		
		unsigned int cur_seed = 0;
		// building each seed at a time
		seedSet.clear();
		while(cur_seed < k) {
			max_index = getMaxIndexBaseline(n, node_weight);
			if (max_index == -1) break; // all sets have been covered 
			seedSet.push_back(max_index);
			coverage += node_weight[max_index];
			edge_list = node_hyperedges[max_index];
			for (i = 0; i < edge_list.size(); i++) {
				if (edge_removed[edge_list[i].first]) continue;
				alive = false;
				coverage_weight = edge_list[i].second; // 当前的seed 所覆盖的RR sets中的某一个RR set的associated weight
				node_list = hyperedges[edge_list[i].first];// 一个 RR set 中的节点
				for (j = 0; j < node_list.size(); j++) {
					if (node_list[j].second > 0) { // only process those nodes with a non-zero marginal for this hyperedge
						diff = node_list[j].second - coverage_weight;
						if (diff <= 0) {
							node_weight[node_list[j].first] -= node_list[j].second;
							node_list[j].second = 0;
						} else {
							node_weight[node_list[j].first] -= coverage_weight;
							node_list[j].second = diff;
						}
						if (!alive && node_list[j].second > 0) alive = true; // still marginal 	gain achievable from this hyperedge
					}
				}
				if (!alive) edge_removed[edge_list[i].first] = true;
			}
			cur_seed++;
		}

		return coverage;
	}

	double compute_cov_r2(vector<vector<std::pair<int,double> > > hyperedges,vector<vector<std::pair<int,double> > > node_hyperedges,vector<double> node_hyperedge_weight, 
	unsigned int k, vector<int>& seedSet) {
		unsigned int i, j;
		int diff, coverage_weight;
		bool alive;
		vector<pair<int,double> > edge_list, node_list;

		vector<int> node_weight(n,0);
		for (i = 0; i < n; i++){
			node_weight[i] = node_hyperedge_weight[i];
		}
		long long numEdge = hyperedges.size();

		vector<bool> edge_removed(numEdge, false);	
		unsigned int cur_seed = 0;
		int coverage = 0;
		
		while(cur_seed < seedSet.size()) {
			coverage += node_weight[seedSet[cur_seed]];
			edge_list = node_hyperedges[seedSet[cur_seed]]; 
			for (i = 0; i < edge_list.size(); i++) {
				if (edge_removed[edge_list[i].first]) continue;
				alive = false;
				coverage_weight = edge_list[i].second;
				node_list = hyperedges[edge_list[i].first];
				for (j = 0; j < node_list.size(); j++) {
					if (node_list[j].second > 0) { // only process those nodes with a non-zero marginal for this hyperedge
						diff = node_list[j].second - coverage_weight;
						if (diff <= 0) {
							node_weight[node_list[j].first] -= node_list[j].second;
							node_list[j].second = 0;
						} else {
							node_weight[node_list[j].first] -= coverage_weight;
							node_list[j].second = diff;
						}
						if (!alive && node_list[j].second > 0) alive = true; // still marginal gain achievable from this hyperedge
					}
				}
				if (!alive) edge_removed[edge_list[i].first] = true;
			}
			cur_seed++;
		}

		return coverage;
	}

	double generate_TRU_set_with_blocker(vector<int>&blocker,int deadline)
	{
		vector<int> isBlocker(n, 0);
		for (int node : blocker) {
			isBlocker[node] = 1;
		}
		int source_node = sfmt_genrand_uint32(&sfmtSeed) % n; 
		if(isRumor[source_node]) return U(deadline);
		if(isBlocker[source_node]) return 0;

		priority_queue<Node> Q;
		vector<int> AT(n, 10000);
		vector<int> I(n, 0);

		AT[source_node] = 0;
		Q.push({source_node, 0});
		
		while (!Q.empty()) {
			Node current = Q.top();
			Q.pop();

			int q = current.id;
			int q_AT = current.activate_time;
			if (I[q] == 1 || isBlocker[q]) continue; 
			if(isRumor[q]) return U(deadline-AT[q]);
			I[q] = 1;

			for (int i = 0; i < gT_reverse[q].size(); i++) {
				int tmp= sfmt_genrand_uint32(&sfmtSeed) % precomputedPoisson.size();
				int diff_delay=precomputedPoisson[tmp]; 
				int u = gT_reverse[q][i];
				if (isBlocker[u]) continue;

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
		double x_max=U(deadline-1);
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

	

	void generate_sample_TRtree(vector<vector<int>>& sample,vector<vector<int>>& sample_UB, double &sum, int deadline, vector<int>& AT)
	{
		priority_queue<Node> Q;
		
		vector<vector<int>> pred(n+1); 
		vector<int> I(n, 0); 
	
		for (int u : rumorSet) {
			AT[u] = 0;
			Q.push({u, 0});
    	}

		while (!Q.empty()) 
		{
			Node current = Q.top();
			Q.pop();
			int q = current.id;
			int q_AT = current.activate_time;
			if (I[q] == 1) continue;

			for(int i=0;i<gT[q].size();i++)
			{
				int tmp= sfmt_genrand_uint32(&sfmtSeed) % precomputedPoisson.size();
				int diff_delay=precomputedPoisson[tmp];  
				int u=gT[q][i];
				if(sfmt_genrand_real1(&sfmtSeed) <= probT[q][i] && AT[q] + diff_delay <= deadline)
				{
					sample[q].push_back(u);
					if (AT[q] + diff_delay < AT[u]) {
						AT[u] = AT[q] + diff_delay;
						pred[u].clear();
						pred[u].push_back(q);
						Q.push({u, AT[u]});
					} else if (AT[q] + diff_delay == AT[u]) {
						pred[u].push_back(q);
					}
				}
			}	

			I[q] = 1;
			if(!isRumor[q]) sum+=U(deadline-AT[q]);
			if(pred[q].size()!=0)
			{
				int pre=pred[q][0];
				sample_UB[pre].push_back(q);
			}

    	}
	}

	//under the sample ranlization to estimate the upper bound and original
	double rep_approximation(Argument& arg){
		vector<int>isUB;
		isUB.resize(n, 0);
		for (int node : UB_seedSet)
		{
			isUB[node] = 1;
		}
		double sum_ori=0;
		double sum_up=0;
		for(int i=0;i<10000;i++){
			
			double ori_ulitity=0;// excluding S
			vector<vector<int>> sample;
			sample.resize(n);

			vector<vector<int>> sample_UB;
			sample_UB.resize(n);
			vector<int> AT(n, 10000); 
			generate_sample_TRtree(sample, sample_UB, ori_ulitity, arg.T, AT);

			int decrease_upper=0, decrease_ori=0;
			queue<int> q;
			double aft_ulitity=0;//表示没有ub以后ori function的值
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
					aft_ulitity+=U(arg.T-AT[v]);
					// activated_num_ori++;
					vis3[v]=1;
					
				}		
			}
			decrease_ori=ori_ulitity-aft_ulitity;
			//cout<<"decrease_ori: "<<decrease_ori<<endl;
			sum_ori+=decrease_ori;
			
		
			vector<int> vis2;
			vis2.resize(n, 0);
			for (auto x : UB_seedSet)
				q.push(x), vis2[x]=1, decrease_upper+=U(arg.T-AT[x]);
			while (!q.empty())
			{
				int x = q.front();
				q.pop();
				for(int v:sample_UB[x]){
					if(vis2[v]==0){
						q.push(v);
						decrease_upper+=U(arg.T-AT[v]);
						vis2[v]=1;
					}
				}
	
			}
			//cout<<"decrease_upper is: "<<decrease_upper<<endl;
			sum_up+=decrease_upper;
			
			if(decrease_ori>decrease_upper){
				break;
			}
		}

		double t1=sum_ori/(double)10000;

		double t2=sum_up/(double)10000;

		double tmp=pow((1-arg.gamma),2)/(double)(pow((1+arg.gamma),2));
		const double e = exp(1);
		tmp=tmp*(1-1/e-arg.epsilon)*t1/t2;
		cout<<"the gap is: "<<t1/t2<<endl;


		return tmp;
	}

	double determine_better_solution(Argument& arg)
	{
		int T=arg.T;
		double gamma=arg.gamma;
		double delta=1.0/(double)n;
		double inf_after_block_LB = stopping_rule(LB_seedSet, T, gamma, delta);
		std::cout << "inf_after_block_LB: " << inf_after_block_LB << std::endl;

		double inf_after_block_UB = stopping_rule(UB_seedSet, T, gamma, delta);
		std::cout << "inf_after_block_UB: " << inf_after_block_UB << std::endl;
	
		if(inf_after_block_LB<inf_after_block_UB) 
		{
			seedSet=LB_seedSet;
			return inf_after_block_LB;
		}
		else
		{
			seedSet=UB_seedSet;
			return inf_after_block_UB;
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

	void SandTCIM(Argument& arg)
	{
		double epsilon=arg.epsilon;
		double delta=1.0/(double)n;
		int targetSize = arg.k;
		int T=arg.T;
		int X_max=U(T-1);
		const double e = exp(1);
		const double approx = 1 - 1.0 / e;
		double gamma=arg.gamma;

		const double alpha = sqrt(log(6.0 / delta));
		const double beta = sqrt((1 - 1.0 / e) * (logcnk(n - rumorSet.size(), targetSize) + log(6.0 / delta)));
		size_t numRbase, maxNumR, numIter;

		numRbase = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta)); // changable!
		double OPT_L=1;  // changable!
		
		maxNumR= size_t(2.0 * n * X_max * pow2((1 - 1 / e) * sqrt(log(6.0 / delta)) + sqrt((1 - 1.0 / e) * (logcnk(n - rumorSet.size(), targetSize) + log(6.0 / delta)))) / OPT_L / pow2(epsilon)) + 1;

		numIter = (size_t)log2(maxNumR / numRbase) + 1;
		const double a1 = log(numIter * 3.0 / delta);
		const double a2 = log(numIter * 3.0 / delta);

		double x_max=(double)U(arg.T-1);
		for (auto idx = 0; idx < numIter; idx++)
		{
			auto numR = numRbase << idx;
			cout<<"numR is: "<<numR<<endl;

			if(arg.ws_flag==0 && arg.pr_flag==1){//SandTCIM-WS
				while(Lhyperedges_1.size()<numR) US_TRW_R1(arg);
				while(Lhyperedges_2.size()<numR) US_TRW_R2(arg);
			}
			else if(arg.ws_flag==1 && arg.pr_flag==0){//SandTCIM-PR
				while(Lhyperedges_1.size()<numR) WS_TRW_wo_PR_R1(arg);
				while(Lhyperedges_2.size()<numR) WS_TRW_wo_PR_R2(arg);
			}
			else if(arg.ws_flag==0 && arg.pr_flag==0){//SandTCIM-
				while(Lhyperedges_1.size()<numR) US_TRW_wo_PR_R1(arg);
				while(Lhyperedges_2.size()<numR) US_TRW_wo_PR_R2(arg);
			}
			else if(arg.ws_flag==1 && arg.pr_flag==1){//SandTCIM
				while(Lhyperedges_1.size()<numR) WS_TRW_R1(arg);
				while(Lhyperedges_2.size()<numR) WS_TRW_R2(arg);
			}
			
			double cov_ub;
			double cov_r2;
			double lowerSelect;
			double upperOPT;

			//LOWER BOUND
			cov_ub = buildSeedSet(Lhyperedges_1, node_Lhyperedges_1, node_Lhyperedge_weight_1, arg.k, LB_seedSet);
			// cout<<"cov_ub: "<<cov_ub<<endl;
			cov_r2=compute_cov_r2(Lhyperedges_2,node_Lhyperedges_2,node_Lhyperedge_weight_2, arg.k, LB_seedSet);
			lowerSelect = pow2(sqrt(cov_r2/x_max + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
			upperOPT = pow2(sqrt(cov_ub/x_max + a2 / 2.0) + sqrt(a2 / 2.0));
			auto approx_LB = lowerSelect / upperOPT;


			//UPPER BOUND
			cov_ub = buildSeedSet(Uhyperedges_1, node_Uhyperedges_1, node_Uhyperedge_weight_1, arg.k, UB_seedSet);
			// cout<<"cov_ub: "<<cov_ub<<endl;
			cov_r2=compute_cov_r2(Uhyperedges_2,node_Uhyperedges_2,node_Uhyperedge_weight_2, arg.k, UB_seedSet);
			lowerSelect = pow2(sqrt(cov_r2/x_max + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
			upperOPT = pow2(sqrt(cov_ub/x_max + a2 / 2.0) + sqrt(a2 / 2.0));
			auto approx_UB = lowerSelect / upperOPT;
			cout<<"LB approximation ratio: "<<approx_LB<<endl;
			cout<<"UB approximation ratio: "<<approx_UB<<endl;
			if (approx_LB >= approx - epsilon && approx_UB >= approx - epsilon)
			{
				break;
			}

		}
	}

	double GB(Argument& arg)
	{
		int k=arg.k;
		int T=arg.T;
		vector<int>isBlocker(n,0);
		double min,min_id;
		for(int i=0;i<k;i++)
		{
			double ob_value;
			min=999999999,min_id=0;
			for(int j=0;j<n;j++)
			{
				if(isRumor[j] || isBlocker[j]) continue;
				seedSet.push_back(j);
				ob_value=MC_based_estimate_with_removal(10000, T, seedSet);
				seedSet.pop_back();
				if(ob_value<min){
					min=ob_value;
					min_id=j;
				}
			}
			seedSet.push_back(min_id);
			isBlocker[min_id]=1;
			cout<<"blocker is : "<<min_id<<endl;
		}
		return min;
		
	}

	//use US for selecting source node
	double SandTCIM_US(Argument& arg)
	{
		double epsilon=arg.epsilon;
		double delta=1.0/(double)n;
		int targetSize = arg.k;
		int T=arg.T;
		int X_max=U(T-1);
		const double e = exp(1);
		const double approx = 1 - 1.0 / e;
		double gamma=arg.gamma;

		const double alpha = sqrt(log(6.0 / delta));
		const double beta = sqrt((1 - 1.0 / e) * (logcnk(n - rumorSet.size(), targetSize) + log(6.0 / delta)));
		size_t numRbase, maxNumR, numIter;

		numRbase = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta)); // changable!
		double OPT_L=1;  // changable!
		
		maxNumR= size_t(2.0 * n * X_max * pow2((1 - 1 / e) * sqrt(log(6.0 / delta)) + sqrt((1 - 1.0 / e) * (logcnk(n - rumorSet.size(), targetSize) + log(6.0 / delta)))) / OPT_L / pow2(epsilon)) + 1;

		numIter = (size_t)log2(maxNumR / numRbase) + 1;
		const double a1 = log(numIter * 3.0 / delta);
		const double a2 = log(numIter * 3.0 / delta);

		double x_max=(double)U(arg.T-1);
		for (auto idx = 0; idx < numIter; idx++)
		{
			auto numR = numRbase << idx;
			cout<<"numR is: "<<numR<<endl;
			//R1
			while(Lhyperedges_1.size()<numR) WS_TRW_R1(arg);
			//R2
			while(Lhyperedges_2.size()<numR) WS_TRW_R2(arg);
			double cov_ub;
			double cov_r2;
			double lowerSelect;
			double upperOPT;

			//LOWER BOUND
			cov_ub = buildSeedSet(Lhyperedges_1, node_Lhyperedges_1, node_Lhyperedge_weight_1, arg.k, LB_seedSet);
			// cout<<"cov_ub: "<<cov_ub<<endl;
			cov_r2=compute_cov_r2(Lhyperedges_2,node_Lhyperedges_2,node_Lhyperedge_weight_2, arg.k, LB_seedSet);
			lowerSelect = pow2(sqrt(cov_r2/x_max + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
			upperOPT = pow2(sqrt(cov_ub/x_max + a2 / 2.0) + sqrt(a2 / 2.0));
			auto approx_LB = lowerSelect / upperOPT;


			//UPPER BOUND
			cov_ub = buildSeedSet(Uhyperedges_1, node_Uhyperedges_1, node_Uhyperedge_weight_1, arg.k, UB_seedSet);
			// cout<<"cov_ub: "<<cov_ub<<endl;
			cov_r2=compute_cov_r2(Uhyperedges_2,node_Uhyperedges_2,node_Uhyperedge_weight_2, arg.k, UB_seedSet);
			lowerSelect = pow2(sqrt(cov_r2/x_max + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
			upperOPT = pow2(sqrt(cov_ub/x_max + a2 / 2.0) + sqrt(a2 / 2.0));
			auto approx_UB = lowerSelect / upperOPT;
			cout<<"LB approximation ratio: "<<approx_LB<<endl;
			cout<<"UB approximation ratio: "<<approx_UB<<endl;
			if (approx_LB >= approx - epsilon && approx_UB >= approx - epsilon)
			{
				double inf_after_block_LB = stopping_rule(LB_seedSet, T, gamma, delta);
				std::cout << "inf_after_block_LB: " << inf_after_block_LB << std::endl;

				double inf_after_block_UB = stopping_rule(UB_seedSet, T, gamma, delta);
				std::cout << "inf_after_block_UB: " << inf_after_block_UB << std::endl;
			
				if(inf_after_block_LB<inf_after_block_UB) 
				{
					seedSet=LB_seedSet;
					return inf_after_block_LB;
				}
				else
				{
					seedSet=UB_seedSet;
					return inf_after_block_UB;
				}
			}

		}
		//i=i_max
		double inf_after_block_LB = stopping_rule(LB_seedSet, T, gamma, delta);
		double inf_after_block_UB = stopping_rule(UB_seedSet, T, gamma, delta);
		if(inf_after_block_LB<inf_after_block_UB) 
		{
			seedSet=LB_seedSet;
			return inf_after_block_LB;
		}
		else
		{
			seedSet=UB_seedSet;
			return inf_after_block_UB;
		}

	}


	
};