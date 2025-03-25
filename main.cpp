#include "Graph.h"
#include "Utility.h"

using namespace std;

// Usage: [0]exe [1]input_graph [2]k
int main(int argc, const char *argv[])
{
    if (argc < 2)
    {
        cout << "\t Usage: [0]exe [1]input_graph [2]k" << endl;
        exit(1);
    }

    string file_path = argv[1];
    ui K = atoi(argv[2]);

    if (K < 1)
    {
        cout << "Error: k must be a positive integer." << endl;
        exit(1);
    }

    // load graph
    Graph *graph = new Graph(K);
    graph->load_graph(file_path);

    // heuristic signed kplex
    graph->heu_signed_kplex();

    // find_signed_kplex
    graph->find_signed_kplex();

    // print result
    graph->print_result(true);

    delete graph;
    return 0;
}