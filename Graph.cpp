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
#include "MyBitset.h"

static long long REBUILD_TIME = 0;
static long long HEU_TIME = 0;
static long long CTCP_TIME = 0;
static long long BNB_TIME = 0;
static long long TOT_TIME = 0;

Graph::Graph(const int _k)
{
    K = _k;
    kplex.clear();
    N = M = PM = NM = 0;
    n = m = pm = nm = 0;
    lb = ub = 0;
    edges = nullptr;
    rev_edges = nullptr;
    pstart = pend = nullptr;
    esign = nullptr;
    degree = vis = nullptr;
    tri_cnt = nullptr;
}

Graph::~Graph()
{
    kplex.clear();

    delete[] pstart;
    delete[] pend;
    delete[] edges;
    delete[] rev_edges;
    delete[] esign;
    delete[] degree;
    delete[] tri_cnt;
    delete[] vis;

    edges = nullptr;
    rev_edges = nullptr;
    pstart = pend = nullptr;
    esign = nullptr;
    degree = vis = nullptr;
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
    ui *mark = new ui[N];
    memset(mark, 0, sizeof(ui) * N);

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

    N = n;
    M = m;
    PM = pm;
    NM = nm;
    edges = new ui[M];
    rev_edges = new ui[M];
    pstart = new ept[N + 1];
    pend = new ept[N];
    esign = new int[M];
    tri_cnt = new ept[M];
    degree = new ui[N];
    vis = new ui[N];
    v_sel.resize(N); v_sel.flip();

    // construct edges
    ui idx = 0;
    for (ui u = 0; u < N; u++) {
        pstart[u] = idx;
        for (auto e : s_G[u]) {
            edges[idx] = e.first;
            esign[idx] = e.second;
            idx++;
        }
        pend[u] = idx;
    }
    pstart[N] = idx;
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
    for (ui i = 0; i < N; i++) {
        vector<pair<ui, int>> temp;
        for (ept j = pstart[i]; j < pend[i]; j++) temp.push_back({ edges[j], esign[j] });
        sort(temp.begin(), temp.end());
        for (ept j = pstart[i]; j < pend[i]; j++) {
            edges[j] = temp[j - pstart[i]].first;
            esign[j] = temp[j - pstart[i]].second;
        }
    }

    // 计算每条边的反向边
    for (ui u = 0; u < N; u++) pend[u] = pstart[u];
    for (ui u = 0; u < N; u++) {
        for (ept i = pstart[u]; i < pstart[u + 1]; i++) {
            rev_edges[i] = pend[edges[i]]++;
        }
    }
    for (ui u = 0; u < N; u++) {
        if (pend[u] != pstart[u + 1]) {
            printf("u = %d, pend[u] = %ld, pstart[u + 1] = %ld\n", u, pend[u], pstart[u + 1]);
            assert(false);
            exit(1);
        }
    }

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
    CTCP(lb);
    ui *dorder = new ui[N];
    ub = degen(dorder);

    SIGNED_KPLEX_BITSET *signed_kplex_solver_bitset = new SIGNED_KPLEX_BITSET();
    signed_kplex_solver_bitset->allocateMemory(N);

    vector<Edge> vp;
    vp.reserve(m);

    ui *ids = new ui[N];
    ui *mark = new ui[N];
    memset(mark, 0, sizeof(ui) * N);
    ui *Q = new ui[N];
    ui *nei_degree = new ui[N];
    ui *rid = new ui[N];

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
            signed_kplex_solver_bitset->load_graph(s_n, vp);
            signed_kplex_solver_bitset->heu_kPlex(K, kplex, rid_u);
        }

        if (kplex.size() > lb) {
            for (auto &v : kplex) v = ids[v];
            lb = kplex.size();
            CTCP(lb);
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
            signed_kplex_solver_bitset->load_graph(s_n, vp);
            signed_kplex_solver_bitset->heu_kPlex(K, kplex, rid_u);
        }

        if (kplex.size() > lb) {
            for (auto &v : kplex) v = ids[v];
            lb = kplex.size();
            CTCP(lb);
            ub = min(ub, (int)degen(dorder));
            cur = 0;
        }
    }

    delete[] dorder;
    delete[] ids;
    delete[] mark;
    delete[] Q;
    delete[] nei_degree;
    delete[] rid;
    delete signed_kplex_solver_bitset;

    HEU_TIME = t.elapsed();
    TOT_TIME += HEU_TIME;

    cout << "-------------------heu_find_kplex-------------------" << endl;
    cout << "\theu_kplex_size = " << kplex.size() << ",\t time cost = " << integer_to_string(t.elapsed()) << endl;
}
/**
 * @brief find_signed_kplex
 *
 */
void Graph::find_signed_kplex()
{
    dfs_cnt = 0;
    dfs_cnt_1 = 0;
    dfs_cnt_2 = 0;
    dfs_cnt_after_prune = 0;
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
    ui *mark = new ui[N];
    memset(mark, 0, sizeof(ui) * N);
    ui *Q = new ui[n];
    ui *nei_degree = new ui[N];
    ui *rid = new ui[N];
    vector<ui> kplex_bitset = kplex;

    long long BNB_STD_TIME = 0;
    long long BNB_BITSET_TIME = 0;
    while (n > lb) {
        // get u
        get_degree();
        int u = -1;
        for (ui i : v_sel) if (u == -1 || degree[u] > degree[i]) u = i;

        // get g
        ui s_n = 0;
        vp.clear();
        ui rid_u = 0;
        extract_subgraph_two_hop(u, vp, s_n, ids, rid, Q, nei_degree, mark);

        if (s_n > lb) {
            // bnb
            Timer t_bnb;
            t_bnb.restart();

            // standard
            Timer t_standard;
            t_standard.restart();
            // cout << "standard bnb start" << endl;
            signed_kplex_solver->load_graph(s_n, vp);
            signed_kplex_solver->kPlex(K, kplex, (s_n == n) ? -1 : rid_u);
            BNB_STD_TIME += t_standard.elapsed();
            // cout << "standard time cost = " << integer_to_string(BNB_STD_TIME) << endl;

            // bitset
            Timer t_bitset;
            t_bitset.restart();
            // cout << "bitset bnb start" << endl;
            signed_kplex_solver_bitset->load_graph(s_n, vp);
            signed_kplex_solver_bitset->kPlex(K, kplex_bitset, (s_n == n) ? -1 : rid_u);
            BNB_BITSET_TIME += t_bitset.elapsed();
            // cout << "bitset time cost = " << integer_to_string(BNB_BITSET_TIME) << endl;

            // printf("**end bnb: kplex.size() = %d, kplex_bitset.size() = %d\n", kplex.size(), kplex_bitset.size());
            if (kplex.size() != kplex_bitset.size()) {
                printf("kplex.size() = %ld, kplex_bitset.size() = %ld\n", kplex.size(), kplex_bitset.size());
                printf("kplex: ");
                for (auto v : kplex) printf("%d ", v);
                printf("\n");
                printf("kplex_bitset: ");
                for (auto v : kplex_bitset) printf("%d ", v);
                printf("vid:\n");
                for (ui i = 0; i < n; i++) printf("ids[%d] = %d\n", i, ids[i]);
                cout << "ERROR: kplex.size() != kplex_bitset.size()" << endl;
                exit(1);
            }

            // cout << "bnb time cost = " << integer_to_string(t_bnb.elapsed()) << endl;

            BNB_TIME += t_bnb.elapsed();
        }

        if (kplex.size() > lb) {
            lb = kplex.size();
            for (auto &v : kplex) v = ids[v];
        }

        if (s_n == n) break;

        // CTCP
        CTCP(lb, u);
    }

    printf("dfs_cnt = %lld, dfs_cnt_1 = %lld, dfs_cnt_2 = %lld, dfs_cnt_after_prune = %lld\n", dfs_cnt, dfs_cnt_1, dfs_cnt_2, dfs_cnt_after_prune);

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
    for (ui u : v_sel) degree[u] = pend[u] - pstart[u];
}
/**
 * @brief 计算图中所有边的三角形数量
 *
 * @details 根据正确的邻接表，统计图中每条边参与构成的三角形数量
 */
void Graph::get_tricnt()
{
    ui *mark = vis;
    memset(mark, 0, sizeof(ui) * N);
    memset(tri_cnt, 0, sizeof(ept) * M);

    for (auto u : v_sel) {
        for (ept i = pstart[u]; i < pend[u]; i++) mark[edges[i]] = i + 1;

        for (ept i = pstart[u]; i < pend[u]; i++) {
            ui v = edges[i];
            if (u < v) {
                for (ept j = pstart[v]; j < pend[v]; j++) {
                    ui w = edges[j];
                    if (mark[w] && v < w) {
                        ept id_uv = i;
                        ept id_vu = rev_edges[i];
                        ept id_vw = j;
                        ept id_wv = rev_edges[j];
                        ept id_uw = mark[w] - 1;
                        ept id_wu = rev_edges[mark[w] - 1];

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
    for (ui u : v_sel) assert(mark[u] == 0);
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
    for (ui u : v_sel) assert(mark[u] == 0);
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
    for (ui u : v_sel) assert(mark[u] == 0);
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
    for (ui u : v_sel) assert(mark[u] == 0);
#endif
    return 0;
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

    get_degree();

    ui dorder_n = 0, new_size = 0;

    for (ui u : v_sel) if (degree[u] < threshold) dorder[dorder_n++] = u;
    for (ui i = 0; i < dorder_n; i++) {
        ui u = dorder[i];
        degree[u] = 0;
        for (ept j = pstart[u]; j < pend[u]; j++) if (degree[edges[j]] > 0) {
            if ((degree[edges[j]]--) == threshold) dorder[dorder_n++] = i;
        }
    }

    ui UB = n;
    if (dorder_n == n) UB = kplex.size();

    memset(vis, 0, sizeof(ui) * N);
    for (ui u : v_sel) {
        if (degree[u] >= threshold) dorder[dorder_n + (new_size++)] = u;
        else vis[u] = 1;
    }
    assert(dorder_n + new_size == n);

    ListLinearHeap *heap = new ListLinearHeap(N, N - 1);
    if (new_size != 0) {
        heap->init(new_size, N - 1, dorder + dorder_n, degree);
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
void Graph::rebuild_graph(MyBitset &v_del, MyBitset &e_del)
{
    Timer t;
    t.restart();

    v_sel ^= v_del;
    n = v_sel.size();
    m = pm = nm = 0;
    for (ui u : v_sel) {
        ept pos = pstart[u];
        for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) {
            ui v = edges[i];
            if (!v_del[v]) {
                edges[pos] = v;
                esign[pos] = esign[i];
                if (esign[pos] == 1) pm++;
                else nm++;
                pos++;
            }
        }
        pend[u] = pos;
        m += pend[u] - pstart[u];
    }

    // 计算每条边的反向边
    ui *tmp = vis;
    for (ui u : v_sel) tmp[u] = pstart[u];
    for (ui u : v_sel) {
        for (ept i = pstart[u]; i < pend[u]; i++) {
            rev_edges[i] = tmp[edges[i]]++;
        }
    }

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
void Graph::CTCP(int lb, int del_v)
{
    static int last_lb = 0;
    Timer t;
    t.restart();

    int tv = max(0, lb + 1 - K);
    int te = max(0, lb + 1 - 2 * K);

    queue<ui> qv;
    queue<pair<ui, ept>> qe; // from, idx
    get_degree();
    get_tricnt();

    MyBitset v_del(N);
    MyBitset e_del(M);
    // MyBitset in_qe(M);
    ui *mark = vis;
    memset(mark, 0, sizeof(ui) * N);
    if (del_v != -1) qv.push((ui)del_v);
    if (last_lb < lb) {
        for (ui u : v_sel) {
            if (degree[u] < tv) qv.push(u);
            for (ept i = pstart[u]; i < pend[u]; i++) {
                ui v = edges[i];
                if (u < v && tri_cnt[i] < te) qe.push(make_pair(u, i));
            }
        }
    }
    last_lb = lb;
    while (!qv.empty() || !qe.empty()) {
        while (!qe.empty()) {
            auto ue = qe.front(); qe.pop();
            ui u = ue.first;
            ept id_uv = ue.second;
            ui v = edges[id_uv];
            ept id_vu = rev_edges[id_uv];
            assert(e_del[id_uv] == e_del[id_vu]);
            if (v_del[u] || v_del[v] || e_del[id_uv]) continue;
            e_del.insert(id_uv);
            e_del.insert(id_vu);

            if ((degree[u]--) == tv) qv.push(u);
            if ((degree[v]--) == tv) qv.push(v);

            for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) mark[edges[i]] = i + 1;

            for (ept j = pstart[v]; j < pend[v]; j++) if (!e_del[j]) {
                ui w = edges[j];
                if (mark[w]) { // triangle count--
                    ept id_uw = mark[w] - 1;
                    ept id_wu = rev_edges[id_uw];
                    if ((tri_cnt[id_uw]--) == te) qe.push(make_pair(u, id_uw));
                    tri_cnt[id_wu]--;

                    ept id_vw = j;
                    ept id_wv = rev_edges[id_vw];
                    if ((tri_cnt[id_vw]--) == te) qe.push(make_pair(v, id_vw));
                    tri_cnt[id_wv]--;
                }
            }

            for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) mark[edges[i]] = 0;
        }
        if (!qv.empty()) {
            ui u = qv.front(); qv.pop();
            if (v_del[u]) continue; v_del.insert(u);

            for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) mark[edges[i]] = i + 1;

            for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) {
                ui v = edges[i];
                for (ept j = pstart[v]; j < pend[v]; j++) if (!e_del[j]) {
                    ui w = edges[j];
                    if (mark[w] && v < w) { // triangle count--
                        ept id_vw = j;
                        ept id_wv = rev_edges[id_vw];
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
                ept id_vu = rev_edges[id_uv];
                e_del.insert(id_uv);
                e_del.insert(id_vu);
            }
        }
    }
    rebuild_graph(v_del, e_del);
#ifndef NDEBUG
    get_degree();
    get_tricnt();
    for (ui u : v_sel) {
        assert(degree[u] >= tv);
        for (ept i = pstart[u]; i < pend[u]; i++) {
            if (tri_cnt[i] < te) {
                printf("tri_cnt[%u] = %lu, te = %u\n", i, tri_cnt[i], te);
                fflush(stdout);
                assert(false);
            }
        }
    }
#endif
    CTCP_TIME += t.elapsed();
    // #ifndef NDEBUG
    cout << "\tCTCP, T : " << integer_to_string(t.elapsed()) << ",\t n = " << n << ", m = " << m / 2 << endl;
    // #endif
}
/**
 * @brief core-truss co-pruning
 *
 * @param lb 核心阈值
 * @param del_v 指定要删除的顶点,默认为-1表示不指定
 *
 * @details //TODO 删除一个点会删除所有与该点相连的边，从这一点入手进行优化
 * 先不计算边对应的三角形变化，待点缩减完成后再计算
 */
void Graph::CTCP_new(int lb, int del_v)
{
    static int last_lb = 0;
    Timer t;
    t.restart();

    int tv = max(0, lb + 1 - K);
    int te = max(0, lb + 1 - 2 * K);

    queue<ui> qv; // 等待计算三角形计数的点
    queue<ui> qv_2; // 等待缩减的点
    queue<pair<ui, ept>> qe; // from, idx
    get_degree();
    get_tricnt();

    MyBitset v_del(N);
    MyBitset e_del(M);
    ui *mark = vis;
    memset(mark, 0, sizeof(ui) * N);

    if (del_v != -1) qv_2.push((ui)del_v);
    if (last_lb < lb) {
        for (ui u : v_sel) {
            if (degree[u] < tv) qv_2.push(u);
            for (ept i = pstart[u]; i < pend[u]; i++) {
                ui v = edges[i];
                if (u < v && tri_cnt[i] < te) qe.push(make_pair(u, i));
            }
        }
    }
    last_lb = lb;

    while (!qv.empty() || !qe.empty() || !qv_2.empty()) {
        // 直接删除点
        while (!qv_2.empty()) {
            ui u = qv_2.front(); qv_2.pop();
            if (v_del[u]) continue; v_del.insert(u);
            qv.push(u);
            for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) {
                ui v = edges[i];
                if ((degree[v]--) == tv) qv_2.push(v);
            }
        }
        // 删边同时更新三角形计数
        while (!qv.empty()) {
            ui u = qv.front(); qv.pop();

            for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) {
                ui v = edges[i];
                if (v_del[v]) continue;
                for (ept j = i + 1; j < pend[u]; j++) if (!e_del[j]) {
                    ui w = edges[j];
                    if (v_del[w]) continue;
                    ept id_vw = j;
                    ept id_wv = rev_edges[id_vw];
                    if ((tri_cnt[id_vw]--) == te) qe.push(make_pair(v, id_vw));
                    tri_cnt[id_wv]--;
                }
            }

            for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) {
                ept id_uv = i;
                ept id_vu = rev_edges[id_uv];
                e_del.insert(id_uv);
                e_del.insert(id_vu);
            }
        }
        while (!qe.empty()) {
            auto ue = qe.front(); qe.pop();
            ui u = ue.first; ept id_uv = ue.second;
            ui v = edges[id_uv]; ept id_vu = rev_edges[id_uv];
            assert(e_del[id_uv] == e_del[id_vu]);
            if (v_del[u] || v_del[v] || e_del[id_uv]) continue;
            e_del.insert(id_uv);
            e_del.insert(id_vu);

            bool flag = false;
            if ((degree[u]--) == tv) qv_2.push(u), flag = true;
            if ((degree[v]--) == tv) qv_2.push(v), flag = true;

            for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) mark[edges[i]] = i + 1;

            for (ept j = pstart[v]; j < pend[v]; j++) if (!e_del[j]) {
                ui w = edges[j];
                if (mark[w]) { // triangle count--
                    ept id_uw = mark[w] - 1;
                    ept id_wu = rev_edges[id_uw];
                    if ((tri_cnt[id_uw]--) == te) qe.push(make_pair(u, id_uw));
                    tri_cnt[id_wu]--;

                    ept id_vw = j;
                    ept id_wv = rev_edges[id_vw];
                    if ((tri_cnt[id_vw]--) == te) qe.push(make_pair(v, id_vw));
                    tri_cnt[id_wv]--;
                }
            }

            for (ept i = pstart[u]; i < pend[u]; i++) if (!e_del[i]) mark[edges[i]] = 0;

            if (flag) break; // 如果有删点，则暂停进行三角形计数
        }
    }
    rebuild_graph(v_del, e_del);
#ifndef NDEBUG
    get_degree();
    get_tricnt();
    for (ui u : v_sel) {
        assert(degree[u] >= tv);
        for (ept i = pstart[u]; i < pend[u]; i++) {
            if (tri_cnt[i] < te) {
                printf("tri_cnt[%u] = %lu, te = %u\n", i, tri_cnt[i], te);
                fflush(stdout);
                assert(false);
            }
        }
    }
#endif
    CTCP_TIME += t.elapsed();
    // #ifndef NDEBUG
    cout << "\tCTCP, T : " << integer_to_string(t.elapsed()) << ",\t n = " << n << ", m = " << m / 2 << endl;
    // #endif
}