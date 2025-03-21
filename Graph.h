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

    /**
     * @brief 从文件加载图结构
     * @param input_graph 输入文件路径
     * 支持读取标准图文件格式，初始化邻接表等数据结构
     */
    void load_graph(string input_graph);
    /**
     * @brief 主算法入口
     * 执行带符号k-plex挖掘的核心算法流程
     */
    void find_signed_kplex();
    /**
     * @brief 检查给定的顶点集是否构成一个有效的带符号k-plex
     *
     * @param ids_n 顶点集的大小
     * @param ids 顶点集数组,存储顶点的ID
     */
    void check_is_kplex(ui ids_n, ui *ids);
    void print_result(bool print_solution);
    void heu_signed_kplex();

private:
    /**
     * @brief 计算图的退化序
     * @param dorder 输出参数，存储计算得到的退化序
     * @return 图的退化序值
     */
    ui degen(ui *dorder);
    /**
     * @brief 计算图中所有顶点的度数
     * 同时计算正负边的度数分布
     */
    void get_degree();
    /**
     * @brief 计算图中每个顶点的三角形数量
     * 用于后续剪枝优化
     */
    void get_tricnt();
    /**
     * @brief 提取以顶点u为中心的k-plex子图
     * @param u 中心顶点ID
     * @param vp 边集合的引用
     * @param ids_n 子图顶点数量
     * @return 子图的规模
     */
    ui extract_subgraph(ui u, vector<Edge> &vp, ui &ids_n);
    /**
     * @brief 计算顶点u的g值（k-plex核心度）
     * @param u 目标顶点ID
     * @param vp 边集合引用
     * @param ids_n 候选顶点数量
     * @return 顶点u的g值
     */
    ui get_g(ui u, vector<Edge> &vp, ui &ids_n);
    /**
     * @brief 重建图结构
     * @param v_del 顶点删除标记数组
     * @param e_del 边删除标记数组
     * 根据删除标记压缩图结构，优化存储空间
     */
    void rebuild_graph(bool *v_del, bool *e_del);
    /**
     * @brief 核心三角形剪枝算法（Core Triangle-based Pruning）
     * @param tv 顶点删除阈值
     * @param te 边删除阈值
     * @param del_v 顶点删除标记
     * 通过三角形计数进行有效的分支限界剪枝
     */
    void CTCP(int tv, int te, int del_v);
};

#endif