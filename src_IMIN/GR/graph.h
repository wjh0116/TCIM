#define HEAD_INFO
#include "sfmt/SFMT.h"
#include "head.h"

#include "MeasureM.h"
#include "memoryusage.h"

using namespace std;
typedef double (*pf)(int, int);
void handle_error(const char* msg);

class Graph
{
public:    	
	unsigned int n, m;

    vector<vector<int>> gT;
    vector<vector<int>> gT_reverse;//寻找节点的入节点
	vector<vector<double>> probT;
    vector<vector<double>> probT2;//用于复原probT
    vector<vector<int>> edge_visit;
    vector<vector<int>> edge_sample;

    


    enum InfluModel {IC, LT};
    InfluModel influModel;
    void setInfuModel(InfluModel p)
    {
        influModel = p;
        TRACE(influModel == IC);        
        TRACE(influModel == LT);        
    }

    string folder;
    string graph_file;
    void readNM()
    {
        ifstream cin((folder + "attribute.txt").c_str());
        ASSERT(!cin == false);
        string s;
        while (cin >> s)
        {
            if (s.substr(0, 2) == "n=")
            {
                n = atoi(s.substr(2).c_str());
                continue;
            }
            if (s.substr(0, 2) == "m=")
            {
                m = atoi(s.substr(2).c_str());
                continue;
            }
            ASSERT(false);
        }
        TRACE(n, m );
        cin.close();
    }

    void readGraph()
    {
		size_t length;
		int fd = open((graph_file).c_str(), O_RDWR);
		if (fd == -1)
			handle_error("open");
		struct stat sb;
		int rc = fstat(fd, &sb);
		if (rc == -1)
			handle_error("fstat");

		length = sb.st_size;
		auto ptr = static_cast<char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));  //byte by byte
		auto f = ptr;

		int gap = 2 * sizeof(int) + sizeof(double);
        //ASSERT(fin != false);
        unsigned int readCnt = 0;
        for (unsigned int i = 0; i < m; i++)
        {
            readCnt ++;
            unsigned int a, b;
            double p;
			memcpy(&a, f, sizeof(int));
			memcpy(&b, f + sizeof(int), sizeof(int));
			memcpy(&p, f + 2 * sizeof(int), sizeof(double));
			f += gap;

            ASSERT( a < n );
            ASSERT( b < n );
			
			
			probT[a].push_back(p);
            probT2[a].push_back(p);
            edge_visit[a].push_back(0);
            edge_sample[a].push_back(0);
			gT[a].push_back(b);
            gT_reverse[b].push_back(a);
        }        

        ASSERT(readCnt == m);
		rc = munmap(ptr, length);
		close(fd);
        
    }

    Graph(string folder, string graph_file): folder(folder), graph_file(graph_file)
    {		
		readNM();

		gT = vector<vector<int>>(n+2, vector<int>());
        gT_reverse = vector<vector<int>>(n+2, vector<int>());
		probT = vector<vector<double>>(n+2, vector<double>());
        probT2 = vector<vector<double>>(n+2, vector<double>());
        edge_visit = vector<vector<int>>(n+2, vector<int>());
        edge_sample = vector<vector<int>>(n+2, vector<int>());
		readGraph();
    }
};

double sqr(double t)
{
    return t * t;
}

void handle_error(const char* msg) 
{
	perror(msg);
	exit(255);
}

