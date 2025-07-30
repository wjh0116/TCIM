#define HEAD_INFO

#include "sfmt/SFMT.h"
#include "head.h"

class Argument {
public:
    int k;
    unsigned int Rumor_num;
    string dataset;
    string model;
    string res;
    string seedfile;
    string algo;
    int time;
};

#include "graph.h"
#include "infgraph.h"
#include "AG.h"

//static unsigned int rr_num=0;
void OutputSeedSetToFile(vector<int> seed_set, string seedfile)
{
    ofstream of(seedfile, ios::app);
    //of.open(seedfile);
    for (int seed : seed_set)
    {
        of << seed << endl;
    }
    of << endl;
    of.close();
}

void run_with_parameter(InfGraph& g, const Argument& arg)
{
    CP::CP_based_AMR(g, arg);

    //INFO(g.seedSet);
    //OutputSeedSetToFile(g.seedSet, arg.seedfile);
}
void Run(int argn, char** argv)
{
    Argument arg;

    arg.k = 0;

    for (int i = 0; i < argn; i++)
    {
        if (argv[i] == string("-help") || argv[i] == string("--help") || argn == 1)
        {
            cout << "./imm -dataset *** -epsilon *** -k ***  -model IC|LT -seedfile *** -time *** -batch ***" << endl;
            return;
        }
        if (argv[i] == string("-k"))
            arg.k = atoi(argv[i + 1]);
        if (argv[i] == string("-rumorNum"))
            arg.Rumor_num = atoi(argv[i + 1]);
        if (argv[i] == string("-dataset"))
            arg.dataset = argv[i + 1];
        if (argv[i] == string("-algo"))
            arg.algo = argv[i + 1];
        if (argv[i] == string("-model"))
            arg.model = argv[i + 1];
    }
    ASSERT(arg.dataset != "");
    ASSERT(arg.model == "IC" || arg.model == "LT");
    string tempname = arg.dataset.substr(arg.dataset.find_last_of("/") + 1);
    arg.res = "results/res_" + tempname + "_|S|=" + to_string(arg.Rumor_num)+ "_K=" + to_string(arg.k) +"_"+arg.algo;
    arg.dataset = arg.dataset + "/";
    string graph_file;
    if (arg.model == "IC")
        graph_file = arg.dataset + "graph_ic.inf";
    else if (arg.model == "LT")
        graph_file = arg.dataset + "graph_lt.inf";
    else
        ASSERT(false);

    InfGraph g(arg.dataset, graph_file);

    if (arg.model == "IC")
        g.setInfuModel(InfGraph::IC);
    else if (arg.model == "LT")
        g.setInfuModel(InfGraph::LT);
    else
        ASSERT(false);

    run_with_parameter(g, arg);
}


int main(int argn, char** argv)
{
    __head_version = "v1";
    OutputInfo info(argn, argv);

    Run(argn, argv);
}

