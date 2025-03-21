/**
 * @file Graph.cpp
 * @brief 图数据结构的实现文件
 * @author Trae
 * @date 2025-03-21
 */
#include "Graph.h"
#include "Utility.h"
#include "SignedKplex.h"

static long long REBUILD_TIME = 0;
static long long HEU_TIME = 0;
static long long CTCP_TIME = 0;
static long long BNB_TIME = 0;
static long long TOT_TIME = 0;

Graph::Graph(const int _k)
{
    K = _k;
    kplex.clear();
    n = m = pm = nm = 0;
    lb = ub = 0;

    edges = p_edges = n_edges = nullptr;
    pstart = p_pstart = n_pstart = nullptr;
    pend = p_pend = n_pend = nullptr;
    degree = v_rid = vis = nullptr;
    tri_cnt = nullptr;
}

Graph::~Graph()
{
    kplex.clear();

    delete[] pstart;
    delete[] pend;
    delete[] edges;
    delete[] p_pstart;
    delete[] p_pend;
    delete[] p_edges;
    delete[] n_pstart;
    delete[] n_pend;
    delete[] n_edges;
    delete[] degree;
    delete[] tri_cnt;
    delete[] v_rid;
    delete[] vis;

    pstart = pend = nullptr;
    edges = p_edges = n_edges = nullptr;
    p_pstart = p_pend = nullptr;
    n_pstart = n_pend = nullptr;
    degree = v_rid = vis = nullptr;
    tri_cnt = nullptr;
}

/**
 * @brief 检查给定的顶点集是否构成一个有效的带符号k-plex
 *
 * @param ids_n 顶点集的大小
 * @param ids 顶点集数组,存储顶点的ID
 */
void Graph::check_is_kplex(ui ids_n, ui *ids)
{
    // 记录每个顶点在子图中的度数
    ui *deg = new ui[ids_n];
    memset(deg, 0, sizeof(ui) * ids_n);

    // 标记数组,用于快速判断顶点是否在ids中
    ui *mark = new ui[n];
    memset(mark, 0, sizeof(ui) * n);

    // 邻接矩阵存储边的符号
    int *mat = new int[ids_n * ids_n];
    memset(mat, 0, sizeof(int) * ids_n * ids_n);

    // 子图建图
    for (ui i = 0; i < ids_n; i++) mark[ids[i]] = i + 1;
    for (ui i = 0; i < ids_n; i++) {
        ui u = ids[i];
        ui ids_u = mark[u] - 1;

        // 处理正边
        for (ept j = p_pstart[u]; j < p_pend[u]; j++) {
            ui v = p_edges[j];
            if (mark[v]) {
                ui ids_v = mark[v] - 1;
                deg[ids_u]++; deg[ids_v]++;
                mat[ids_u * ids_n + ids_v] = mat[ids_v * ids_n + ids_u] = 1;
            }
        }

        // 处理负边
        for (ept j = n_pstart[u]; j < n_pend[u]; j++) {
            ui v = n_edges[j];
            if (mark[v]) {
                ui ids_v = mark[v] - 1;
                deg[ids_u]++; deg[ids_v]++;
                mat[ids_u * ids_n + ids_v] = mat[ids_v * ids_n + ids_u] = -1;
            }
        }
    }

    // 检查k-plex条件
    for (ui i = 0; i < ids_n; i++) {
        if (deg[i] + K < ids_n) {
            printf("not a signed kplex!!!kplex\n");
            delete[] deg;
            delete[] mark;
            delete[] mat;
            return;
        }
    }

    // 检查三角形的符号平衡性
    for (ui i = 0; i < ids_n; i++) {
        for (ui j = i + 1; j < ids_n; j++) {
            for (ui k = j + 1; k < ids_n; k++) {
                if (mat[i * ids_n + j] && mat[i * ids_n + k] && mat[j * ids_n + k]) {
                    int tri_cn = mat[i * ids_n + j] + mat[i * ids_n + k] + mat[j * ids_n + k];
                    if (tri_cn == 1 || tri_cn == -3) {
                        printf("not a signed kplex!!!signed\n");
                        delete[] mat;
                        delete[] deg;
                        delete[] mark;
                        return;
                    }
                }
            }
        }
    }

    printf("is a signed kplex\n");

    delete[] mat;
    delete[] deg;
    delete[] mark;
}
/**
 * @brief 从文件加载图结构
 *
 * @param input_graph 输入文件路径
 */
void Graph::load_graph(string input_graph)
{
    Timer t;
    t.restart();

    ifstream input_file(input_graph, ios::in);
    if (!input_file.is_open()) {
        cout << "cannot open file : " << input_graph << endl;
        fflush(stdout);
        exit(1);
    }

    input_file >> n >> m;
    map<ui, int> *s_G = new map<ui, int>[n];
    int u, v, flag;
    while (input_file >> u >> v >> flag) {
        if (u == v) continue;
        assert(u >= 0 && u < n);
        assert(v >= 0 && v < n);
        assert(flag == 1 || flag == -1);
        s_G[u].insert(make_pair(v, flag));
        s_G[v].insert(make_pair(u, flag));
    }

    m = pm = nm = 0;
    for (ui i = 0; i < n; i++) {
        for (auto e : s_G[i]) {
            ++m;
            if (e.second == 1)
                ++pm;
            else
                ++nm;
        }
    }
    assert(m == pm + nm);
    assert(m % 2 == 0 && pm % 2 == 0 && nm % 2 == 0);

    input_file.close();

    edges = new ui[m];
    p_edges = new ui[m];
    n_edges = new ui[m];
    tri_cnt = new ept[m];
    pstart = new ept[n + 1];
    p_pstart = new ept[n + 1];
    n_pstart = new ept[n + 1];
    pend = new ept[n];
    p_pend = new ept[n];
    n_pend = new ept[n];
    degree = new ui[n];
    v_rid = new ui[n];
    vis = new ui[n];

    // construct edges
    ui idx = 0, p_idx = 0, n_idx = 0;
    for (ui u = 0; u < n; u++) {
        pstart[u] = idx;
        p_pstart[u] = p_idx;
        n_pstart[u] = n_idx;
        for (auto e : s_G[u]) {
            edges[idx++] = e.first;
            if (e.second == 1)
                p_edges[p_idx++] = e.first;
            if (e.second == -1)
                n_edges[n_idx++] = e.first;
        }
        pend[u] = idx;
        p_pend[u] = p_idx;
        n_pend[u] = n_idx;
    }
    pstart[n] = idx;
    p_pstart[n] = p_idx;
    n_pstart[n] = n_idx;
    pm = p_idx, nm = n_idx;
    assert(idx == m && m == pm + nm);
    assert(m % 2 == 0 && pm % 2 == 0 && nm % 2 == 0);

    delete[] s_G;

    for (ui i = 0; i < n; i++) {
        sort(edges + pstart[i], edges + pend[i]);
        sort(p_edges + p_pstart[i], p_edges + p_pend[i]);
        sort(n_edges + n_pstart[i], n_edges + n_pend[i]);
    }

    for (ui i = 0; i < n; i++) v_rid[i] = i;

    lb = 0, ub = n;

    TOT_TIME += t.elapsed();
#ifndef NPRINT
    cout << "\t -------------------load_graph-------------------" << endl;
    cout << "\t Graph: " << input_graph << endl;
    cout << "\t G : n = " << n << ", m = " << m / 2 << ", pm = " << pm / 2 << ", nm = " << nm / 2 << endl;
    cout << "\t load_graph: time cost = " << integer_to_string(t.elapsed()) << endl;
#endif
}
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
void Graph::heu_signed_kplex()
{
    Timer t;
    t.restart();
    lb = max((int)kplex.size(), 2 * K - 1);
    CTCP(lb + 1 - K, lb + 1 - 2 * K);
    ui *dorder = new ui[n];
    ub = degen(dorder);

    ui *pn = new ui[n];
    ui *in_pn = new ui[n];
    ui *pn_rid = new ui[n];
    ui *heu_kplex = new ui[n];
    ui p_end = 0, n_end = 0, heu_n = 0;
    ui *degree_in_S = new ui[n];
    ui *neighbor = new ui[n];
    memset(in_pn, 0, sizeof(ui) * n);
    memset(neighbor, 0, sizeof(ui) * n);
    memset(degree_in_S, 0, sizeof(ui) * n);
    memset(vis, 0, sizeof(ui) * n);

    ui num = min(n, (ui)10);
    for (ui i = 1; i <= num; i++) {
        // ListLinearHeap* heap = new ListLinearHeap(n, n - 1);
        ui u = dorder[n - i];
        p_end = n_end = 0;

        pn_rid[u] = p_end;
        pn[p_end++] = u;
        vis[u] = 3;
        heu_n = 0;
        heu_kplex[heu_n++] = u;

        for (ept j = p_pstart[u]; j < p_pend[u]; j++) {
            ui v = p_edges[j];
            in_pn[v] = u;
            degree_in_S[p_end] = 1;
            pn_rid[v] = p_end;
            pn[p_end++] = v;
            vis[v] = 1;
        }
        n_end = p_end;
        for (ept j = n_pstart[u]; j < n_pend[u]; j++) {
            ui v = n_edges[j];
            in_pn[v] = u;
            vis[v] = 2;
            degree_in_S[n_end] = 1;
            pn_rid[v] = n_end;
            pn[n_end++] = v;
        }
        printf("p_end = %d, n_end = %d, heu_n = %d\n", p_end, n_end, heu_n);
        while (1) {
#ifndef NDEBUG
            for (ui j = p_end; j < n_end; j++) {
                if (vis[pn[j]] == 3) {
                    assert(degree_in_S[j] + K >= heu_n + 1);
                }
            }
            for (ui j = p_end; j < n_end; j++) {
                if (vis[pn[j]] == 1 || vis[pn[j]] == 2) {
                    assert(degree_in_S[j] + K >= heu_n + 1);
                }
            }
            // ui tot = 0;
            // for (ui j = 0; j < n_end; j++) {
            //     if (vis[pn[j]] == 1 || vis[pn[j]] == 2) tot++;
            // }
            // printf("tot = %d\n", tot);

            for (ui j = 0; j < n_end; j++)
                assert(neighbor[j] == 0);
            for (ui j = 0; j < n_end; j++)
                assert(pn_rid[pn[j]] == j);
#endif

            ui nxt_u = n;
            for (ui j = 0; j < n_end; j++) {
                if (vis[pn[j]] == 1 || vis[pn[j]] == 2) {
                    if (nxt_u == n)
                        nxt_u = j;
                    else if (degree_in_S[nxt_u] < degree_in_S[j])
                        nxt_u = j;
                }
            }
            if (nxt_u == n)
                break;

            nxt_u = pn[nxt_u];
            assert(vis[nxt_u] == 1 || vis[nxt_u] == 2);
            ui u_pn = vis[nxt_u];
            vis[nxt_u] = 3;
            // printf("heu_n = %d\n", heu_n);
            heu_kplex[heu_n++] = nxt_u;

            for (ept j = p_pstart[nxt_u]; j < p_pend[nxt_u]; j++) {
                ui v = p_edges[j];
                if (vis[v] == 1 || vis[v] == 2) {
                    if (in_pn[v] == u) {
                        assert(pn_rid[v] >= 0 && pn_rid[v] < n_end);
                        neighbor[pn_rid[v]] = 1;
                        if ((vis[v] == 1 && u_pn == 1) || (vis[v] == 2 && u_pn == 2)) {
                            degree_in_S[pn_rid[v]]++;
                        }
                        else {
                            vis[v] = 0;
                        }
                    }
                }
                else if (vis[v] == 3) {
                    if (in_pn[v] == u) {
                        assert(pn_rid[v] >= 0 && pn_rid[v] < n_end);
                        neighbor[pn_rid[v]] = 1;
                        degree_in_S[pn_rid[v]]++;
                    }
                }
            }
            for (ept j = n_pstart[nxt_u]; j < n_pend[nxt_u]; j++) {
                ui v = n_edges[j];
                if (vis[v] == 1 || vis[v] == 2 && pn_rid[v] < n_end) {
                    if (in_pn[v] == u) {
                        assert(pn_rid[v] >= 0 && pn_rid[v] < n_end);
                        neighbor[pn_rid[v]] = 1;
                        if ((vis[v] == 1 && u_pn == 2) || (vis[v] == 2 && u_pn == 1)) {
                            degree_in_S[pn_rid[v]]++;
                        }
                        else {
                            vis[v] = 0;
                        }
                    }
                }
                else if (vis[v] == 3) {
                    if (in_pn[v] == u) {
                        assert(pn_rid[v] >= 0 && pn_rid[v] < n_end);
                        neighbor[pn_rid[v]] = 1;
                        degree_in_S[pn_rid[v]]++;
                    }
                }
            }

            bool ub_mark = false;
            for (ui j = 0; j < p_end; j++) {
                if (neighbor[j] == 0) {
                    if (vis[pn[j]] == 3) {
                        if (degree_in_S[j] + K <= heu_n) {
                            ub_mark = true;
                        }
                    }
                }
                neighbor[j] = 0;
            }
            for (ui j = p_end; j < n_end; j++) {
                if (neighbor[j] == 0) {
                    if (vis[pn[j]] == 1 || vis[pn[j]] == 2) {
                        if (degree_in_S[j] + K <= heu_n + 1)
                            vis[pn[j]] = 0;
                    }
                }
                neighbor[j] = 0;
            }
            if (ub_mark)
                break;
        }

#ifndef NDEBUG
        ui vis1Num = 0, vis2Num = 0, vis3Num = 0;
        for (ui j = 0; j < n; j++) {
            if (vis[j] == 1)
                vis1Num++;
            if (vis[j] == 2)
                vis2Num++;
            if (vis[j] == 3)
                vis3Num++;
        }
        printf("vis1 = %d, vis2 = %d, vis3 = %d\n", vis1Num, vis2Num, vis3Num);
#endif

        printf("heu_n = %d\n", heu_n);

        // ui addNum = min(n, (ui)100);
        // // Add as many nodes as possible
        // for (ui j = 1; j <= addNum; j++) {
        //     ui v = dorder[n - j];
        //     if (heu_add_check(heu_n, heu_kplex, v)) {
        //         heu_kplex[heu_n++] = v;
        //     }
        // }
        // printf("heu_n(add nodes) = %d\n", heu_n);

        if (heu_n > kplex.size()) {
            kplex.clear();
            for (ui j = 0; j < heu_n; j++)
                kplex.push_back(heu_kplex[j]);
        }
    }

    if (lb < kplex.size()) {
        lb = kplex.size();
        CTCP(lb + 1 - K, lb + 1 - 2 * K);
    }

    delete[] dorder;
    delete[] pn;
    delete[] in_pn;
    delete[] pn_rid;
    delete[] heu_kplex;
    delete[] degree_in_S;
    delete[] neighbor;
    HEU_TIME = t.elapsed();
    TOT_TIME += HEU_TIME;
#ifndef NPRINT
    cout << "\t -------------------heu_find_kplex-------------------" << endl;
    cout << "\t heu_kplex_size = " << kplex.size() << ",\t time cost = " << integer_to_string(t.elapsed()) << endl;
    cout << "\t after graph reduction, G: n = " << n / 2 << ", m = " << m / 2 << ", pm = " << pm / 2 << ", nm = " << nm << endl;
#endif
}
/**
 * @brief find_signed_kplex
 *
 */
void Graph::find_signed_kplex()
{
    Timer t;
    t.restart();
    if (kplex.size() >= ub) return;

    SIGNED_KPLEX *signed_kplex_solver = new SIGNED_KPLEX();
    signed_kplex_solver->allocateMemory(n, m);

    vector<Edge> vp;
    vp.reserve(m);

    ui *ids = new ui[n];
    ui *mark = new ui[n];
    memset(mark, 0, sizeof(ui) * n);
    ui *Q = new ui[n];
    ui *nei_degree = new ui[n];
    ui *rid = new ui[n];

    while (n > lb) {
        printf("n = %d, lb = %d\n", n, lb);
        // get u
        get_degree();
        ui u = 0;
        for (ui i = 1; i < n; i++) if (degree[u] > degree[i]) u = i;
        // printf("deg = %d, u = %d\n", degree[u], u);

        // get g
        ui s_n = 0;
        vp.clear();
        ui rid_u = 0;
        // rid_u = extract_graph_without_prune(u, vp, s_n, ids);
        extract_subgraph(u, vp, s_n, ids, rid, Q, nei_degree, mark);

        // //输出vp和s_n的大小信息
        // printf("vp size = %d, s_n = %d\n", vp.size(), s_n);
        // for (auto &v : vp)
        //     printf("(%d, %d, %d)\n", v.a, v.b, v.c);

        if (s_n > lb) {
            // bnb
            Timer t_bnb;
            t_bnb.restart();
            signed_kplex_solver->load_graph(s_n, vp);
            signed_kplex_solver->kPlex(K, kplex, (s_n == n) ? -1 : rid_u);
            BNB_TIME += t_bnb.elapsed();
        }

        if (kplex.size() > lb) {
            lb = kplex.size();
            for (auto &v : kplex)
                v = v_rid[ids[v]];
        }

        if (s_n == n) break;

        // CTCP
        CTCP(lb + 1 - K, lb + 1 - 2 * K, u);
    }

    delete signed_kplex_solver;
    TOT_TIME += t.elapsed();
}
/**
 * @brief 打印算法运行结果
 *
 * @param print_solution 是否打印具体的k-plex顶点集
 */
void Graph::print_result(bool print_solution)
{
    cout << "\t -------------------print_result-------------------" << endl;
    cout << "\t CTCP_TIME: " << integer_to_string(CTCP_TIME);
    cout << "\t REBUILD_TIME: " << integer_to_string(REBUILD_TIME);
    cout << "\t HEU_TIME: " << integer_to_string(HEU_TIME);
    cout << ",\t BNB_TIME: " << integer_to_string(BNB_TIME);
    cout << "\t TOT_TIME: " << integer_to_string(TOT_TIME) << endl;
    if (kplex.size() < 2 * K - 1)
        printf("***We can't find a plex larger than 2k-2!! The following is a heuristic solution.\n");
    cout << "\t kplex_size: " << kplex.size() << endl;
    if (print_solution) {
        cout << "\t -------------------maximum k-plex-------------------" << endl;
        for (auto u : kplex) printf("%d ", u);
        printf("\n");
    }
    printf("----------------------------------------------------\n");
}
/**
 * @brief 计算图中所有顶点的度数
 *
 * @details 根据正确的邻接表，统计图中所有顶点的度数
 */
void Graph::get_degree()
{
    for (ui i = 0; i < n; i++) {
        degree[i] = pstart[i + 1] - pstart[i];
    }
}
/**
 * @brief 计算图中所有边的三角形数量
 *
 * @details 根据正确的邻接表，统计图中每条边参与构成的三角形数量
 */
void Graph::get_tricnt()
{
    ui *mark = vis;
    memset(mark, 0, sizeof(ui) * n);
    memset(tri_cnt, 0, sizeof(ept) * m);

    for (ui u = 0; u < n; u++) {
        for (ept i = pstart[u]; i < pend[u]; i++) mark[edges[i]] = i + 1;

        for (ept i = pstart[u]; i < pend[u]; i++) {
            ui v = edges[i];
            if (u < v) {
                for (ept j = pstart[v]; j < pend[v]; j++) {
                    ui w = edges[j];
                    if (mark[w] && v < w) {
                        ept id_uv = i;
                        ept id_vu = pstart[v] + find(edges + pstart[v], edges + pend[v], u);
                        ept id_vw = j;
                        ept id_wv = pstart[w] + find(edges + pstart[w], edges + pend[w], v);
                        ept id_uw = mark[w] - 1;
                        ept id_wu = pstart[w] + find(edges + pstart[w], edges + pend[w], u);
#ifndef NDEBUG
                        ept id_uv1 = pstart[u] + find(edges + pstart[u], edges + pend[u], v);
                        ept id_vw1 = pstart[v] + find(edges + pstart[v], edges + pend[v], w);
                        ept id_uw1 = pstart[u] + find(edges + pstart[u], edges + pend[u], w);
                        assert(id_uv == id_uv1);
                        assert(id_vw == id_vw1);
                        assert(id_uw == id_uw1);
#endif
                        tri_cnt[id_uv]++;
                        tri_cnt[id_vu]++;
                        tri_cnt[id_vw]++;
                        tri_cnt[id_wv]++;
                        tri_cnt[id_uw]++;
                        tri_cnt[id_wu]++;
                    }
                }
            }
        }

        for (ept i = pstart[u]; i < pend[u]; i++) mark[edges[i]] = 0;
    }
}
/**
  * @brief 提取顶点u的2-hop子图,并进行剪枝
  *
  * @param u 需要提取子图的中心顶点
  * @param vp 用于存储子图边的vector
  * @param ids_n 返回子图的顶点数
  * @param sub_v_rid 顶点映射表，维护原图顶点到子图顶点的映射
  * @return ui 返回中心顶点的映射后的新ID
  *
  * @details
  * 中心顶点标记为1，邻居标记为2，2-hop邻居标记为3
  * 该函数提取以u为中心的2-hop子图,并进行剪枝优化:
  * 1. 首先提取u的直接邻居
  * 2. 计算邻居节点在子图中的度数,删除度数过小的节点
  * 3. 提取2-hop邻居并继续剪枝
  * 4. 重新编号并构建子图的边集
  */
ui Graph::extract_subgraph(ui u, vector<Edge> &vp, ui &ids_n, ui *ids, ui *rid, ui *Q, ui *nei_degree, ui *mark)
{
#ifndef NDEBUG
    for (ui i = 0; i < n; i++) assert(mark[i] == 0);
#endif
    vp.clear();

    // 添加中心顶点u
    ids_n = 0;
    ids[ids_n++] = u;
    mark[u] = 1;

    // 处理u的直接邻居
    for (ept i = pstart[u]; i < pend[u]; i++) {
        ui v = edges[i];
        mark[v] = 2;
        ids[ids_n++] = v;
    }

    // Any two adjacent vertices must have at least l−2k common neighbors
    ui Q_n = 0;
    for (ui i = 1; i < ids_n; i++) {
        ui v = ids[i];
        nei_degree[v] = 0;
        for (ept j = pstart[v]; j < pend[v]; j++) if (mark[edges[j]] == 2) ++nei_degree[v];
        if (nei_degree[v] + 2 * K <= kplex.size()) Q[Q_n++] = v;
    }
    for (ui i = 0; i < Q_n; i++) {
        ui v = Q[i];
        mark[v] = 10; // deleted
        for (ept j = pstart[v]; j < pend[v]; j++) if (mark[edges[j]] == 2) {
            if ((nei_degree[edges[j]]--) + 2 * K == kplex.size() + 1) {
                Q[Q_n++] = edges[j];
            }
        }
    }
    assert(Q_n <= ids_n - 1);

    // 中心顶点的上界
    if (ids_n - 1 - Q_n + K <= kplex.size()) {
        for (ui i = 0; i < ids_n; i++) mark[ids[i]] = 0;
        ids_n = 0;
        return 0;
    }

    // 处理2-hop邻居
    ui nei_size = ids_n;
    for (ui i = 1; i < nei_size; i++) if (mark[ids[i]] == 2) {
        ui v = ids[i];
        for (ept j = pstart[v]; j < pend[v]; j++) {
            ui w = edges[j];
            if (!mark[w]) {
                ids[ids_n++] = w;
                mark[w] = 3;
                nei_degree[w] = 1;
            }
            else if (mark[w] == 3) nei_degree[w]++;
        }
    }

    ui new_size = 1;
    for (ui i = 1; i < nei_size; i++) {
        if (mark[ids[i]] == 10) mark[ids[i]] = 0;
        else ids[new_size++] = ids[i];
    }

    assert(new_size + Q_n == nei_size);
    ui old_nei_size = nei_size;
    nei_size = new_size;
    // Any two non-adjacent vertices must have at least l−2k+2 common neighbors
    for (ui i = old_nei_size; i < ids_n; i++) {
        if (nei_degree[ids[i]] + 2 * K <= kplex.size() + 2) mark[ids[i]] = 0;
        else ids[new_size++] = ids[i];
    }
    ids_n = new_size;

#ifndef NDEBUG
    assert(mark[ids[0]] == 1);
    for (ui i = 1; i < nei_size; i++)
        assert(mark[ids[i]] == 2);
    for (ui i = nei_size; i < ids_n; i++)
        assert(mark[ids[i]] == 3);
#endif

    // 重新编号
    for (ui i = 0; i < ids_n; i++) {
        assert(mark[ids[i]]);
        rid[ids[i]] = i;
    }

    for (ui i = 0; i < ids_n; i++) {
        u = ids[i];
        for (ept i = p_pstart[u]; i < p_pend[u]; i++) {
            ui v = p_edges[i];
            if (mark[v] && u < v) {
                vp.push_back(Edge(rid[u], rid[v], 1));
            }
        }
        for (ept i = n_pstart[u]; i < n_pend[u]; i++) {
            ui v = n_edges[i];
            if (mark[v] && u < v) {
                vp.push_back(Edge(rid[u], rid[v], -1));
            }
        }
    }

    for (ui i = 0; i < ids_n; i++) mark[ids[i]] = 0;
#ifndef NDEBUG
    for (ui i = 0; i < n; i++) assert(mark[i] == 0);
#endif
    return 0;
}
/**
  * @brief 提取顶点u的2-hop子图,不进行剪枝
  *
  * @param u 需要提取子图的中心顶点
  * @param vp 用于存储子图边的vector
  * @param ids_n 返回子图的顶点数
  * @param sub_v_rid 顶点映射表，维护原图顶点到子图顶点的映射
  * @return ui 返回中心顶点的映射后的新ID
  */
ui Graph::extract_graph_without_prune(ui u, vector<Edge> &vp, ui &ids_n, ui *ids)
{
    Timer t;
    t.restart();

    bool *v_sel = new bool[n];
    memset(v_sel, 0, sizeof(bool) * n);

    v_sel[u] = 1;
    for (ept i = pstart[u]; i < pend[u]; i++) {
        ui v = edges[i];
        v_sel[v] = 1;
        for (ept j = pstart[v]; j < pend[v]; j++) {
            ui w = edges[j];
            v_sel[w] = 1;
        }
    }

    // 构建新旧顶点ID的映射
    ui rid_u = n;
    ui *rid = new ui[n];
    ui cnt = 0;

    for (ui i = 0; i < n; i++) if (v_sel[i]) {
        if (u == i) rid_u = cnt;
        ids[cnt] = i;
        rid[i] = cnt++;
    }
    ids_n = cnt;
    printf("cnt = %d\n", cnt);

    for (ui u = 0; u < n; u++) if (v_sel[u]) {
        for (ept i = p_pstart[u]; i < p_pend[u]; i++) {
            ui v = p_edges[i];
            if (v_sel[v]) {
                if (u < v) {
                    vp.push_back(Edge(rid[u], rid[v], 1));
                }
            }
        }
        for (ept i = n_pstart[u]; i < n_pend[u]; i++) {
            ui v = n_edges[i];
            if (v_sel[v]) {
                if (u < v) {
                    vp.push_back(Edge(rid[u], rid[v], -1));
                }
            }
        }
    }

    delete[] rid;
    delete[] v_sel;

#ifndef NDEBUG
    cout << "\t extract_graph_without_prune, T : " << integer_to_string(t.elapsed()) << ",\t n=" << ids_n << endl;
#endif
    assert(rid_u != n);
    return rid_u;
}
/**
 * @brief 计算图的退化序列和上界
 *
 * @param dorder 用于存储退化序列的数组
 * @return 返回图的上界值
 *
 * @details 得到了图的上界，得到了dorder数组，包含图的dorder序
 */
ui Graph::degen(ui *dorder)
{
    Timer t;
    t.restart();
    ui threshold = kplex.size() + 1 > K ? kplex.size() + 1 - K : 0;

    for (ui i = 0; i < n; i++) degree[i] = pstart[i + 1] - pstart[i];

    ui dorder_n = 0, new_size = 0;

    for (ui i = 0; i < n; i++) if (degree[i] < threshold) dorder[dorder_n++] = i;
    for (ui i = 0; i < dorder_n; i++) {
        ui u = dorder[i];
        degree[u] = 0;
        for (ept j = pstart[u]; j < pstart[u + 1]; j++) if (degree[edges[j]] > 0) {
            if ((degree[edges[j]]--) == threshold) dorder[dorder_n++] = i;
        }
    }

    ui UB = n;
    if (dorder_n == n) UB = kplex.size();

    memset(vis, 0, sizeof(ui) * n);
    for (ui i = 0; i < n; i++) {
        if (degree[i] >= threshold) dorder[dorder_n + (new_size++)] = i;
        else vis[i] = 1;
    }
    assert(dorder_n + new_size == n);

    ListLinearHeap *heap = new ListLinearHeap(n, n - 1);

    if (new_size != 0) {
        heap->init(new_size, new_size - 1, dorder + dorder_n, degree);
        ui max_core = 0;
        UB = 0;
        for (ui i = 0; i < new_size; i++) {
            ui u, key;
            heap->pop_min(u, key);
            if (key > max_core) max_core = key;
            dorder[dorder_n + i] = u;

            ui t_UB = max_core + K;
            if (new_size - i < t_UB) t_UB = new_size - i;
            if (t_UB > UB) UB = t_UB;

            vis[u] = 1;

            for (ept j = pstart[u]; j < pend[u]; j++) if (vis[edges[j]] == 0) {
                heap->decrement(edges[j], 1);
            }
        }

        printf("*** Degen(unsigned): max_core: %u, UB: %u, Time: %s (microseconds)\n", max_core, UB, integer_to_string(t.elapsed()).c_str());
    }
    delete heap;

    return UB;
}

/**
 * @brief 重建图结构
 *
 * @param v_del 顶点删除标记数组
 * @param e_del 边删除标记数组
 *
 */
void Graph::rebuild_graph(bool *v_del, bool *e_del)
{
    Timer t;
    t.restart();

    ui *rid = new ui[n];
    ui new_n = 0;
    for (ui i = 0; i < n; i++) {
        if (!v_del[i]) {
            v_rid[new_n] = v_rid[i];
            rid[i] = new_n++;
        }
    }

    new_n = 0;
    ept pos = 0, p_pos = 0, n_pos = 0;
    pstart[0] = p_pstart[0] = n_pstart[0] = 0;
    for (ui u = 0; u < n; u++) if (!v_del[u]) {
        ept p_pointer = p_pstart[u], n_pointer = n_pstart[u];
        for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) {
            ui v = edges[i];
            if (!v_del[v]) {
                edges[pos++] = rid[v];

                // 寻找边(u,v)在正负边表中的位置
                while (p_pointer < p_pend[u] && p_edges[p_pointer] < v) p_pointer++;
                while (n_pointer < n_pend[u] && n_edges[n_pointer] < v) n_pointer++;

                if (p_pointer < p_pend[u] && p_edges[p_pointer] == v) {
                    p_edges[p_pos++] = rid[v];
                    p_pointer++;
                }
                else if (n_pointer < n_pend[u] && n_edges[n_pointer] == v) {
                    n_edges[n_pos++] = rid[v];
                    n_pointer++;
                }
                else {
                    throw std::runtime_error("rebuild_graph error");
                }
            }
        }
        pend[new_n] = pos;
        p_pend[new_n] = p_pos;
        n_pend[new_n] = n_pos;
        new_n++;
    }

    assert(pos % 2 == 0 && p_pos % 2 == 0 && n_pos % 2 == 0);
    n = new_n;
    m = pos;
    pm = p_pos;
    nm = n_pos;

    for (ui u = 1; u <= n; u++) {
        pstart[u] = pend[u - 1];
        p_pstart[u] = p_pend[u - 1];
        n_pstart[u] = n_pend[u - 1];
    }

    delete[] rid;

    REBUILD_TIME += t.elapsed();
}

/**
 * @brief core-truss co-pruning
 *
 * @param tv 顶点度数阈值,度数小于tv的顶点将被删除
 * @param te 三角形数量阈值,三角形数小于te的边将被删除
 * @param del_v 指定要删除的顶点,默认为-1表示不指定
 *
 */
void Graph::CTCP(int tv, int te, int del_v)
{
    static int last_tv = 0;
    Timer t;
    t.restart();

    tv = max(0, tv);
    te = max(0, te);

    queue<ui> qv;
    queue<pair<ui, ept>> qe; // from, idx
    for (ui i = 0; i < n; i++)
        degree[i] = pstart[i + 1] - pstart[i];
    get_tricnt();

    bool *v_del = new bool[n];
    bool *e_del = new bool[m];
    ui *mark = new ui[n];
    memset(v_del, 0, sizeof(bool) * n);
    memset(e_del, 0, sizeof(bool) * m);
    memset(mark, 0, sizeof(ui) * n);
    if (del_v != -1)
        qv.push((ui)del_v);
    if (last_tv < tv) {
        for (ui u = 0; u < n; u++) {
            if (degree[u] < tv)
                qv.push(u);
            for (ept i = pstart[u]; i < pend[u]; i++) {
                ui v = edges[i];
                if (u < v && tri_cnt[i] < te)
                    qe.push(make_pair(u, i));
            }
        }
    }
    last_tv = tv;
    while (!qv.empty() || !qe.empty()) {
        while (!qe.empty()) {
            auto ue = qe.front();
            qe.pop();
            ui u = ue.first;
            ept id_uv = ue.second;
            ui v = edges[id_uv];
            ept id_vu = pstart[v] + find(edges + pstart[v], edges + pend[v], u);
            assert(e_del[id_uv] == e_del[id_vu]);
            if (v_del[u] || v_del[v] || e_del[id_uv]) continue;
            e_del[id_uv] = 1;
            e_del[id_vu] = 1;

            if ((degree[u]--) == tv) qv.push(u);
            if ((degree[v]--) == tv) qv.push(v);

            for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) mark[edges[i]] = i + 1;

            for (ept j = pstart[v]; j < pend[v]; j++) if (!e_del[j]) {
                ui w = edges[j];
                if (mark[w]) { // triangle count--
                    ept id_uw = mark[w] - 1;
                    ept id_wu = pstart[w] + find(edges + pstart[w], edges + pend[w], u);
                    if ((tri_cnt[id_uw]--) == te) qe.push(make_pair(u, id_uw));
                    tri_cnt[id_wu]--;

                    ept id_vw = j;
                    ept id_wv = pstart[w] + find(edges + pstart[w], edges + pend[w], v);
                    if ((tri_cnt[id_vw]--) == te) qe.push(make_pair(v, id_vw));
                    tri_cnt[id_wv]--;
                }
            }

            for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) mark[edges[i]] = 0;
        }
        if (!qv.empty()) {
            ui u = qv.front();
            qv.pop();
            if (v_del[u]) continue;
            v_del[u] = 1;

            for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) mark[edges[i]] = i + 1;

            for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) {
                ui v = edges[i];
                for (ept j = pstart[v]; j < pend[v]; j++) if (!e_del[j]) {
                    ui w = edges[j];
                    if (mark[w] && v < w) { // triangle count--
                        ept id_vw = j;
                        ept id_wv = pstart[w] + find(edges + pstart[w], edges + pend[w], v);
                        if ((tri_cnt[id_vw]--) == te) qe.push(make_pair(v, id_vw));
                        tri_cnt[id_wv]--;
                    }
                }
            }

            for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) mark[edges[i]] = 0;

            for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) {
                ui v = edges[i];
                if ((degree[v]--) == tv) qv.push(v);
                ept id_uv = i;
                ept id_vu = pstart[v] + find(edges + pstart[v], edges + pend[v], u);
                e_del[id_uv] = 1;
                e_del[id_vu] = 1;
            }
        }
    }

    rebuild_graph(v_del, e_del);

    delete[] v_del;
    delete[] e_del;
    delete[] mark;
#ifndef NDEBUG
    get_degree();
    get_tricnt();
    for (ui u = 0; u < n; u++) {
        assert(degree[u] >= tv);
        for (ept i = pstart[u]; i < pend[u]; i++) {
            assert(tri_cnt[i] >= te);
        }
    }
#endif
    CTCP_TIME += t.elapsed();
#ifndef NDEBUG
    cout << "\t CTCP, T : " << integer_to_string(t.elapsed()) << ",\t n = " << n << ", m = " << m / 2 << endl;
#endif
}