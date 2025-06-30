#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"
#include "MyBitset.h"

class Graph
{
private:
    vector<ui> kplex;
    int K;

    ui N, M, PM, NM;
    ui n;   // number of vertices
    ept m;  // number of edges
    ept pm; // number of positive edges
    ept nm; // number of negative edges

    /*Store edges in linear arrays*/
    ept *pstart; // start of edge number of a point
    ept *pend;   // end of edge number of a point
    ui *edges;   // edges
    ui *rev_edges;
    int *esign;

    ui *degree;   // degree of a point
    ept *tri_cnt;

    ui *vis;

    int lb, ub;

    MyBitset v_sel;

public:
    Graph(const int _k);
    ~Graph();

    /**
     * @brief 从文件加载图结构
     *
     * @param input_graph 输入文件路径
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
    /**
     * @brief 打印算法运行结果
     *
     * @param print_solution 是否打印具体的k-plex顶点集
     */
    void print_result(bool print_solution);
    /**
     * @brief 启发式算法寻找signed k-plex
     *
     * 该函数使用启发式方法在图中搜索signed k-plex。
     * 主要步骤包括:
     * 1. 使用度数排序选择初始顶点
     * 2. 从初始顶点开始贪心扩展子图
     * 3. 检查和维护符号平衡性约束
     * 4. 更新当前找到的最大k-plex
     */
    void heu_signed_kplex();

private:
    /**
     * @brief 计算图的退化序列和上界
     *
     * @param dorder 用于存储退化序列的数组
     * @return 返回图的上界值
     *
     * @details 得到了图的上界，得到了dorder数组，包含图的dorder序
     */
    ui degen(ui *dorder);
    /**
     * @brief 计算图中所有顶点的度数
     *
     * @details 根据正确的邻接表，统计图中所有顶点的度数
     */
    void get_degree();
    /**
     * @brief 计算图中所有边的三角形数量
     *
     * @details 根据正确的邻接表，统计图中每条边参与构成的三角形数量
     */
    void get_tricnt();
    /**
     * @brief 提取以顶点u为中心的k-plex子图
     * @param u 中心顶点ID
     * @param vp 边集合的引用
     * @param ids_n 子图顶点数量
     * @return 子图中u的id
     */
    ui extract_subgraph_two_hop(ui u, vector<Edge> &vp, ui &ids_n, ui *ids, ui *rid, ui *Q, ui *nei_degree, ui *mark);
    /**
     * @brief 提取以顶点u为中心的k-plex子图
     * @param u 中心顶点ID
     * @param vp 边集合的引用
     * @param ids_n 子图顶点数量
     * @return 子图中u的id
     */
    ui extract_subgraph_one_hop(ui u, vector<Edge> &vp, ui &ids_n, ui *ids, ui *rid, ui *Q, ui *nei_degree, ui *mark);
    /**
      * @brief 提取顶点u的2-hop子图,不进行剪枝
      *
      * @param u 需要提取子图的中心顶点
      * @param vp 用于存储子图边的vector
      * @param ids_n 返回子图的顶点数
      * @param sub_v_rid 顶点映射表，维护原图顶点到子图顶点的映射
      * @return ui 返回中心顶点的映射后的新ID
      */
    ui extract_graph_without_prune(ui u, vector<Edge> &vp, ui &ids_n, ui *ids);
    /**
     * @brief 重建图结构
     *
     * @param v_del 顶点删除标记数组
     * @param e_del 边删除标记数组
     *
     */
    void rebuild_graph(MyBitset &v_del, MyBitset &e_del);
    /**
     * @brief core-truss co-pruning
     *
     * @param tv 顶点度数阈值,度数小于tv的顶点将被删除
     * @param te 三角形数量阈值,三角形数小于te的边将被删除
     * @param del_v 指定要删除的顶点,默认为-1表示不指定
     *
     */
    void CTCP(int lb, int del_v = -1);
    /**
     * @brief core-truss co-pruning
     *
     * @param tv 顶点度数阈值,度数小于tv的顶点将被删除
     * @param te 三角形数量阈值,三角形数小于te的边将被删除
     * @param del_v 指定要删除的顶点,默认为-1表示不指定
     *
     */
    void CTCP_new(int lb, int del_v = -1);
};

#endif