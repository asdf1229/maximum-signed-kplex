#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"

class Graph
{
private:
    vector<ui> kplex;
    int K;

    ui n;   // number of vertices
    ept m;  // number of edges
    ept pm; // number of positive edges
    ept nm; // number of negative edges

    /*Store edges in linear arrays*/
    ept *pstart; // start of edge number of a point
    ept *pend;   // end of edge number of a point
    ui *edges;   // edges

    ept *p_pstart; // start of positive edge number of a point
    ept *p_pend;   // end of positive edge number of a point
    ui *p_edges;   // positive edges

    ept *n_pstart; // start of negative edge number of a point
    ept *n_pend;   // end of negative edge number of a point
    ui *n_edges;   // negative edges

    ui *degree;   // degree of a point
    ui *p_degree; // positive degree of a point
    ui *n_degree; // negative degree of a point
    ept *tri_cnt;

    ui *v_rid;
    ui *vis;

    bool *v_del;
    bool *e_del;

    int lb, ub;
    int s_n;

public:
    Graph(const int _k);
    ~Graph();

    void load_graph(string input_graph);
    void find_signed_kplex();
    void check_is_kplex(ui ids_n, ui *ids);
    void print_result(bool print_solution);
    void heu_signed_kplex();

private:
    ui degen(ui *dorder);
    void get_degree();
    void get_tricnt();
    ui extract_subgraph(ui u, vector<Edge> &vp, ui &ids_n);
    ui get_g(ui u, vector<Edge> &vp, ui &ids_n);
    void rebuild_graph(bool *v_del, bool *e_del);
    void CTCP(int tv, int te, int del_v);
};

#endif