/**
 * @file Graph.cpp
 * @brief 图数据结构的实现文件
 * @author asdf1229
 * @date 2025-03-21
 */
#include "Graph.h"
#include "Utility.h"
#include "SignedKplex.h"
#include "SignedKplex_new.h"

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
    edges = nullptr;
    pstart = pend = nullptr;
    esign = nullptr;
    degree = v_rid = vis = nullptr;
    tri_cnt = nullptr;
}

Graph::~Graph()
{
    kplex.clear();

    delete[] pstart;
    delete[] pend;
    delete[] edges;
    delete[] esign;
    delete[] degree;
    delete[] tri_cnt;
    delete[] v_rid;
    delete[] vis;

    edges = nullptr;
    pstart = pend = nullptr;
    esign = nullptr;
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

        for (ept j = pstart[u]; j < pend[u]; j++) {
            ui v = edges[j];
            if (mark[v]) {
                ui ids_v = mark[v] - 1;
                int sign_uv = esign[j];
                deg[ids_u]++; deg[ids_v]++;
                mat[ids_u * ids_n + ids_v] = mat[ids_v * ids_n + ids_u] = sign_uv;
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

    m = 0;
    ept pm = 0, nm = 0;
    for (ui i = 0; i < n; i++) {
        for (auto e : s_G[i]) {
            ++m;
            if (e.second == 1) ++pm;
            else ++nm;
        }
    }
    assert(m == pm + nm);
    assert(m % 2 == 0 && pm % 2 == 0 && nm % 2 == 0);

    input_file.close();

    edges = new ui[m];
    pstart = new ept[n + 1];
    pend = new ept[n];
    esign = new int[m];
    tri_cnt = new ept[m];
    degree = new ui[n];
    v_rid = new ui[n];
    vis = new ui[n];

    // construct edges
    ui idx = 0;
    for (ui u = 0; u < n; u++) {
        pstart[u] = idx;
        for (auto e : s_G[u]) {
            edges[idx] = e.first;
            esign[idx] = e.second;
            idx++;
        }
        pend[u] = idx;
    }
    pstart[n] = idx;
    // 统计正负边数量
    pm = nm = 0;
    for (ui i = 0; i < m; i++) {
        if (esign[i] == 1) pm++;
        else nm++;
    }
    assert(idx == m && m == pm + nm);
    assert(m % 2 == 0 && pm % 2 == 0 && nm % 2 == 0);
    delete[] s_G;

    // 对每个顶点的邻接表按照边和符号一起排序
    for (ui i = 0; i < n; i++) {
        vector<pair<ui, int>> temp;
        for (ept j = pstart[i]; j < pend[i]; j++) temp.push_back({ edges[j], esign[j] });
        sort(temp.begin(), temp.end());
        for (ept j = pstart[i]; j < pend[i]; j++) {
            edges[j] = temp[j - pstart[i]].first;
            esign[j] = temp[j - pstart[i]].second;
        }
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
    lb = max((int)kplex.size(), 2 * K - 2);
    CTCP(lb + 1 - K, lb + 1 - 2 * K);
    ui *dorder = new ui[n];
    ub = degen(dorder);

    SIGNED_KPLEX *signed_kplex_solver = new SIGNED_KPLEX();
    signed_kplex_solver->allocateMemory(n, m);
    SIGNED_KPLEX_BITSET *signed_kplex_solver_bitset = new SIGNED_KPLEX_BITSET();
    signed_kplex_solver_bitset->allocateMemory(n);

    vector<ui> kplex_bitset;

    vector<Edge> vp;
    vp.reserve(m);

    ui *ids = new ui[n];
    ui *mark = new ui[n];
    memset(mark, 0, sizeof(ui) * n);
    ui *Q = new ui[n];
    ui *nei_degree = new ui[n];
    ui *rid = new ui[n];

    // 1-hop heuristic
    ui cur = 0;
    ui num1 = min(n, (ui)20);
    for (ui i = 1; i <= num1; i++) {
        if (n < cur + 1) break;
        ui u = dorder[n - (++cur)];
        if (degree[u] + 1 <= lb) continue;

        // get g
        ui s_n = 0;
        ui rid_u = 0;
        extract_subgraph_one_hop(u, vp, s_n, ids, rid, Q, nei_degree, mark);

        if (s_n > lb) {
            signed_kplex_solver->load_graph(s_n, vp);
            signed_kplex_solver->heu_kPlex(K, kplex, rid_u);

            signed_kplex_solver_bitset->load_graph(s_n, vp);
            signed_kplex_solver_bitset->heu_kPlex(K, kplex_bitset, rid_u);
#ifndef NDEBUG
            printf("std_heu_kplex_size = %lu, bitset_heu_kplex_size = %lu\n", kplex.size(), kplex_bitset.size());
#endif
        }

        if (kplex.size() > lb) {
            for (auto &v : kplex) v = v_rid[ids[v]];
            lb = kplex.size();
            CTCP(lb + 1 - K, lb + 1 - 2 * K);
            ub = min(ub, (int)degen(dorder));
            cur = 0;
        }
    }

    // 2-hop heuristic
    ui num2 = min(n, (ui)20);
    for (ui i = 1; i <= num2; i++) {
        if (n < cur + 1) break;
        ui u = dorder[n - (++cur)];
        if (degree[u] + K <= lb) continue;

        // get g
        ui s_n = 0;
        ui rid_u = 0;
        extract_subgraph_two_hop(u, vp, s_n, ids, rid, Q, nei_degree, mark);
        if (s_n > lb) {
            signed_kplex_solver->load_graph(s_n, vp);
            signed_kplex_solver->heu_kPlex(K, kplex, rid_u);

            signed_kplex_solver_bitset->load_graph(s_n, vp);
            signed_kplex_solver_bitset->heu_kPlex(K, kplex_bitset, rid_u);
#ifndef NDEBUG
            printf("std_heu_kplex_size = %lu, bitset_heu_kplex_size = %lu\n", kplex.size(), kplex_bitset.size());
#endif
        }

        if (kplex.size() > lb) {
            for (auto &v : kplex) v = v_rid[ids[v]];
            lb = kplex.size();
            CTCP(lb + 1 - K, lb + 1 - 2 * K);
            ub = min(ub, (int)degen(dorder));
            cur = 0;
        }
    }

    // kplex.push_back(1000000);
    // printf("kplex_size = %lu, lb = %d\n", kplex.size(), lb);
    // fflush(stdout);
    // if (kplex.size() > lb) {
    //     // for (auto &v : kplex) v = v_rid[ids[v]];
    //     lb = kplex.size();
    //     CTCP(lb + 1 - K, lb + 1 - 2 * K);
    //     ub = min(ub, (int)degen(dorder));
    //     cur = 0;
    // }

    delete[] dorder;
    delete[] ids;
    delete[] mark;
    delete[] Q;
    delete[] nei_degree;
    delete[] rid;
    delete signed_kplex_solver;
    delete signed_kplex_solver_bitset;

    HEU_TIME = t.elapsed();
    TOT_TIME += HEU_TIME;

    cout << "-------------------heu_find_kplex-------------------" << endl;
    printf("std_heu_kplex_size = %lu, bitset_heu_kplex_size = %lu\n", kplex.size(), kplex_bitset.size());
    cout << "\theu_kplex_size = " << kplex.size() << ",\t time cost = " << integer_to_string(t.elapsed()) << endl;
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
    SIGNED_KPLEX_BITSET *signed_kplex_solver_bitset = new SIGNED_KPLEX_BITSET();
    signed_kplex_solver_bitset->allocateMemory(n);

    vector<Edge> vp;
    vp.reserve(m);

    ui *ids = new ui[n];
    memset(ids, 0, sizeof(ui) * n);
    ui *mark = new ui[n];
    memset(mark, 0, sizeof(ui) * n);
    ui *Q = new ui[n];
    ui *nei_degree = new ui[n];
    ui *rid = new ui[n];
    vector<ui> kplex_bitset = kplex;

    long long BNB_STD_TIME = 0;
    long long BNB_BITSET_TIME = 0;

    while (n > lb) {
        // get u
        get_degree();
        ui u = 0;
        for (ui i = 1; i < n; i++) if (degree[u] > degree[i]) u = i;

        // get g
        ui s_n = 0;
        vp.clear();
        ui rid_u = 0;
        extract_subgraph_two_hop(u, vp, s_n, ids, rid, Q, nei_degree, mark);

        if (s_n > lb) {
            // printf("start bnb: s_n = %d, s_m = %d, u = %d\n", s_n, vp.size(), v_rid[u]);
            // bnb
            Timer t_bnb;
            t_bnb.restart();

            // standard
            Timer t_standard;
            t_standard.restart();
            // cout << "standard bnb start" << endl;
            signed_kplex_solver->load_graph(s_n, vp);
            signed_kplex_solver->kPlex(K, kplex, (s_n == n) ? -1 : rid_u);
            BNB_STD_TIME = t_standard.elapsed();
            // cout << "standard time cost = " << integer_to_string(BNB_STD_TIME) << endl;

            // bitset
            Timer t_bitset;
            t_bitset.restart();
            // cout << "bitset bnb start" << endl;
            signed_kplex_solver_bitset->load_graph(s_n, vp);
            signed_kplex_solver_bitset->kPlex(K, kplex_bitset, (s_n == n) ? -1 : rid_u);
            BNB_BITSET_TIME = t_bitset.elapsed();
            // cout << "bitset time cost = " << integer_to_string(BNB_BITSET_TIME) << endl;

            // printf("**end bnb: kplex.size() = %d, kplex_bitset.size() = %d\n", kplex.size(), kplex_bitset.size());
            if (kplex.size() != kplex_bitset.size()) {
                printf("kplex: ");
                for (auto v : kplex) printf("%d ", v);
                // for (auto v : kplex) printf("%d ", v_rid[ids[v]]);
                printf("\n");
                printf("kplex_bitset: ");
                for (auto v : kplex_bitset) printf("%d ", v);
                // for (auto v : kplex_bitset) printf("%d ", v_rid[ids[v]]);
                // printf("\n");
                printf("vid:\n");
                for (ui i = 0; i < n; i++) printf("ids[%d] = %d\n", i, v_rid[ids[i]]);
                // printf("\n");
                // printf("vp:\n");
                // for (auto e : vp) printf("(%d, %d, %d)\n", v_rid[ids[e.a]], v_rid[ids[e.b]], e.c);
                // 终止程序
                cout << "ERROR: kplex.size() != kplex_bitset.size()" << endl;
                exit(1);
            }

            BNB_TIME += t_bnb.elapsed();
        }

        if (kplex.size() > lb) {
            lb = kplex.size();
            for (auto &v : kplex) v = v_rid[ids[v]];
        }

        if (s_n == n) break;

        // CTCP
        CTCP(lb + 1 - K, lb + 1 - 2 * K, u);
    }

    // signed_kplex_solver->print_dfs_cnt();
    delete signed_kplex_solver;
    delete signed_kplex_solver_bitset;
    delete[] ids;
    delete[] mark;
    delete[] Q;
    delete[] nei_degree;
    delete[] rid;
    TOT_TIME += t.elapsed();

    printf("BNB_STD_TIME = %lld, BNB_BITSET_TIME = %lld\n", BNB_STD_TIME, BNB_BITSET_TIME);
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
                        // // 判断是否是平衡三角形
                        // int sign_uv = esign[i];
                        // int sign_vw = esign[j];
                        // int sign_uw = esign[mark[w] - 1];
                        // int tri_cn = sign_uv + sign_vw + sign_uw;
                        // assert(tri_cn == 3 || tri_cn == 1 || tri_cn == -1 || tri_cn == -3);
                        // if (tri_cn == 1 || tri_cn == -3) continue;

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
 * @brief 提取并剪枝顶点u的2-hop子图
 *
 * @param u 中心顶点
 * @param vp 存储子图边的容器
 * @param ids_n 子图的顶点数量
 * @param ids 存储子图顶点ID的数组
 * @param rid 顶点重编号映射数组
 * @param Q 用于BFS的队列数组
 * @param nei_degree 存储顶点在子图中度数的数组
 * @param mark 顶点标记数组(1:中心顶点, 2:直接邻居, 3:2-hop邻居)
 * @return ui 返回中心顶点在子图中的新ID
 *
 * @details
 * 该函数执行以下步骤:
 * 1. 将中心顶点u加入子图并标记为1
 * 2. 添加u的直接邻居并标记为2
 * 3. 计算邻居节点度数,剪枝度数不满足条件的节点
 * 4. 添加2-hop邻居并标记为3,继续剪枝
 * 5. 对保留的顶点重新编号,构建子图的边集
 */
ui Graph::extract_subgraph_two_hop(ui u, vector<Edge> &vp, ui &ids_n, ui *ids, ui *rid, ui *Q, ui *nei_degree, ui *mark)
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
        for (ept j = pstart[u]; j < pend[u]; j++) {
            ui v = edges[j];
            if (mark[v] && u < v) vp.push_back(Edge(rid[u], rid[v], esign[j]));
        }
    }

    for (ui i = 0; i < ids_n; i++) mark[ids[i]] = 0;
#ifndef NDEBUG
    for (ui i = 0; i < n; i++) assert(mark[i] == 0);
#endif
    return 0;
}
/**
 * @brief 提取并剪枝顶点u的1-hop子图
 *
 */
ui Graph::extract_subgraph_one_hop(ui u, vector<Edge> &vp, ui &ids_n, ui *ids, ui *rid, ui *Q, ui *nei_degree, ui *mark)
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
        if (nei_degree[v] + K <= kplex.size()) Q[Q_n++] = v;
    }
    for (ui i = 0; i < Q_n; i++) {
        ui v = Q[i];
        mark[v] = 10; // deleted
        for (ept j = pstart[v]; j < pend[v]; j++) if (mark[edges[j]] == 2) {
            if ((nei_degree[edges[j]]--) + K == kplex.size() + 1) {
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

    ui nei_size = ids_n;
    ui new_size = 1;
    for (ui i = 1; i < nei_size; i++) {
        if (mark[ids[i]] == 10) mark[ids[i]] = 0;
        else ids[new_size++] = ids[i];
    }
    ids_n = new_size;

    // 重新编号
    for (ui i = 0; i < ids_n; i++) {
        assert(mark[ids[i]]);
        rid[ids[i]] = i;
    }

    for (ui i = 0; i < ids_n; i++) {
        u = ids[i];
        for (ept j = pstart[u]; j < pend[u]; j++) {
            ui v = edges[j];
            if (mark[v] && u < v) vp.push_back(Edge(rid[u], rid[v], esign[j]));
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
        for (ept j = pstart[u]; j < pend[u]; j++) {
            ui v = edges[j];
            if (v_sel[v] && u < v) vp.push_back(Edge(rid[u], rid[v], esign[j]));
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
    ept pos = 0;
    pstart[0] = 0;
    for (ui u = 0; u < n; u++) if (!v_del[u]) {
        for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) {
            ui v = edges[i];
            if (!v_del[v]) {
                edges[pos] = rid[v];
                esign[pos] = esign[i];
                pos++;
            }
        }
        pend[new_n] = pos;
        new_n++;
    }

    assert(pos % 2 == 0);
    n = new_n;
    m = pos;

    // 统计正负边数量
    pm = nm = 0;
    for (ui i = 0; i < m; i++) {
        if (esign[i] == 1) pm++;
        else nm++;
    }

    for (ui u = 1; u <= n; u++)  pstart[u] = pend[u - 1];

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
    if (del_v != -1) qv.push((ui)del_v);
    if (last_tv < tv) {
        for (ui u = 0; u < n; u++) {
            if (degree[u] < tv) qv.push(u);
            for (ept i = pstart[u]; i < pend[u]; i++) {
                ui v = edges[i];
                if (u < v && tri_cnt[i] < te) qe.push(make_pair(u, i));
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
    // #ifndef NDEBUG
    cout << "\tCTCP, T : " << integer_to_string(t.elapsed()) << ",\t n = " << n << ", m = " << m / 2 << endl;
    // #endif
}