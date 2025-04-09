#define HEAD_INFO

#include "sfmt/SFMT.h"
#include "head.h"

class Argument {
public:
    unsigned int k;
    string dataset;
    string res;
    int Rumor_num;
    int T; // deadline
    double gamma;
    double epsilon;
    string algo;
    int ws_flag;
    int pr_flag;
};

#include "graph.h"
#include "infgraph.h"
#include "Sandwich.h"

void run_with_parameter(InfGraph& g,  Argument& arg)
{
    Sand::SandTCIM(g, arg);
}
void Run(int argn, char** argv)
{
    Argument arg;
    arg.k = 0;

    for (int i = 0; i < argn; i++)
    {
        if (argv[i] == string("-k"))
            arg.k = atoi(argv[i + 1]);
        if (argv[i] == string("-T"))
            arg.T = atoi(argv[i + 1]);
        if (argv[i] == string("-dataset"))
            arg.dataset = argv[i + 1];
        if (argv[i] == string("-algo"))
            arg.algo = argv[i + 1];
        if (argv[i] == string("-rumorNum"))
            arg.Rumor_num = atoi(argv[i + 1]);
        if (argv[i] == string("-gamma"))
            arg.gamma = stod(argv[i + 1]);
        if (argv[i] == string("-epsilon"))
            arg.epsilon = stod(argv[i + 1]);
    }
    if(arg.algo=="SandTCIM-WS")
    {
        arg.pr_flag=1;
        arg.ws_flag=0;
    }
    if(arg.algo=="SandTCIM-PR")
    {
        arg.pr_flag=0;
        arg.ws_flag=1;
    }
    if(arg.algo=="SandTCIM-")
    {
        arg.pr_flag=0;
        arg.ws_flag=0;
    }
    if(arg.algo=="SandTCIM")
    {
        arg.pr_flag=1;
        arg.ws_flag=1;
    }
    ASSERT(arg.dataset != "");
    string temp_name = arg.dataset.substr(arg.dataset.find_last_of("/") + 1);
    arg.res = "results/res_" + temp_name +"_algo=" + arg.algo + "_|S|=" + to_string(arg.Rumor_num) + "_K=" + to_string(arg.k) + "_T=" + to_string(arg.T) 
    + "_epsilon=" + to_string(arg.epsilon) + "_gamma=" + to_string(arg.gamma);
    arg.dataset = arg.dataset + "/";
    string graph_file;
    graph_file = arg.dataset + "graph_ic.inf";
    InfGraph g(arg.dataset, graph_file);// graph reading
    run_with_parameter(g, arg);
}


int main(int argn, char** argv)
{
    __head_version = "v1";
    OutputInfo info(argn, argv);
    Run(argn, argv);
}

