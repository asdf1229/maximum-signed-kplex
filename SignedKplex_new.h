#ifndef _SIGNED_KPLEX_BITSET_
#define _SIGNED_KPLEX_BITSET_

#include "Utility.h"
#include "Timer.h"
#include "MyBitset.h"
#include "LinearHeap.h"

#define _SECOND_ORDER_PRUNING_
using Set = MyBitset;

#define PRINT_SET(x)  print_set((x), #x)

class SIGNED_KPLEX_BITSET
{
private:
    ui K;
    ui n;
    ept m;
    AdjacentMatrix matrix, non_matrix;
    AdjacentMatrix p_matrix, n_matrix;
    AdjacentMatrix non_p_matrix, non_n_matrix;
    ui N; // matrix size
    int ori_u;

    ui best_solution_size;
    Set best_solution;
    vector<ui> degree;
    vector<ui> degree_in_S;
    vector<ui> array_n_1;
    vector<ui> array_n_2;

    int addedVertex;

public:
    SIGNED_KPLEX_BITSET()
    {
        N = K = n = m = 0;
        best_solution_size = 0;
        addedVertex = -1;
    }

    ~SIGNED_KPLEX_BITSET()
    {
    }

    void allocateMemory(ui _N)
    {
        N = 0;
        best_solution = Set(_N);
    }

    void load_graph(ui _n, const vector<Edge> &vp)
    {
        n = _n;
        m = vp.size();

        matrix.resize(n);
        p_matrix.resize(n);
        n_matrix.resize(n);

        non_matrix.resize(n);
        non_p_matrix.resize(n);
        non_n_matrix.resize(n);
        non_matrix.flip();
        non_p_matrix.flip();
        non_n_matrix.flip();

        degree.resize(n);
        degree_in_S.resize(n);
        array_n_1.resize(n);
        array_n_2.resize(n);

        for (ui i = 0; i < m; i++) {
            int a = vp[i].a, b = vp[i].b, c = vp[i].c;
            assert(a < n);
            assert(b < n);
            assert(c == 1 || c == -1);
            matrix.add_edge(a, b);
            matrix.add_edge(b, a);
            non_matrix.remove_edge(a, b);
            non_matrix.remove_edge(b, a);
            if (c == 1) {
                p_matrix.add_edge(a, b);
                p_matrix.add_edge(b, a);
                non_p_matrix.remove_edge(a, b);
                non_p_matrix.remove_edge(b, a);
            }
            if (c == -1) {
                n_matrix.add_edge(a, b);
                n_matrix.add_edge(b, a);
                non_n_matrix.remove_edge(a, b);
                non_n_matrix.remove_edge(b, a);
            }
        }

#ifndef NDEBUG
        printf("matrix kplex(bitset) load_graph: n=%u, m=%lu\n", n, m);
#endif
    }

    void kPlex(ui K_, vector<ui> &kplex, int choose_u)
    {
        K = K_;
        best_solution_size = kplex.size();
        ori_u = choose_u;
        Set S(n), C(n);
        init_bnb(S, C, choose_u);
        bnb_search_two_hop(S, C);
        if (best_solution_size > kplex.size()) {
            kplex.clear();
            for (int v : best_solution) kplex.push_back(v);
        }
        // printf("SIGNED_KPLEX_BITSET: best_solution_size = %d\n", best_solution_size);
    }

    void heu_kPlex(ui K_, vector<ui> &kplex, int choose_u)
    {
        K = K_;
        best_solution_size = kplex.size();
#ifndef NDEBUG
        printf("bitset heu_kPlex begin: n = %d, K = %u, kplex.size = %lu\n", n, K, kplex.size());
#endif
        // 无符号图中的heu_plex
        heu_kPlex_in_unsigned_graph();

        // 有符号图中的heu_plex
        heu_kPlex_in_signed_graph();

        if (best_solution_size > kplex.size()) {
            kplex.clear();
            for (int v : best_solution) kplex.push_back(v);
        }
    }
private:
    void heu_kPlex_in_unsigned_graph()
    {
#ifndef NDEBUG
        printf("heu_kPlex_in_unsigned_graph begin:\n");
#endif
        // 1. 破坏所有的不平衡三角形，对于每个不平衡三角形，选择度数最小的点，将其删去
        for (ui u = 0; u < n; u++) degree[u] = matrix[u].size();
        Set vis(n);
        for (ui i = 0; i < n; i++) if (!vis[i]) {
            for (ui j = i + 1; j < n; j++) if (!vis[j]) {
                if (!matrix[i][j]) continue;
                for (ui k = j + 1; k < n; k++) if (!vis[k]) {
                    if (!matrix[i][k]) continue;
                    if (!matrix[j][k]) continue;
                    int tri_sum = 0;
                    tri_sum += p_matrix[i][j] + p_matrix[j][k] + p_matrix[i][k];
                    tri_sum -= n_matrix[i][j] + n_matrix[j][k] + n_matrix[i][k];
                    if (tri_sum == 1 || tri_sum == -3) {
                        ui min_degree_vertex = i;
                        if (i == 0) min_degree_vertex = degree[j] < degree[k] ? j : k;
                        else {
                            if (degree[j] < degree[min_degree_vertex]) min_degree_vertex = j;
                            if (degree[k] < degree[min_degree_vertex]) min_degree_vertex = k;
                        }
                        vis.insert(min_degree_vertex);
                        for (auto v : matrix[min_degree_vertex]) degree[v]--;

                        assert(i != j && j != k && i != k);
                        if (min_degree_vertex == i) break;
                        if (min_degree_vertex == j) break;
                        if (min_degree_vertex == k) continue;
                    }
                }
                if (vis[i]) break;
            }
        }

        vis.flip();
        Set in = vis;
        vis.clear(); vis.flip();
        ui new_n = in.size();
#ifndef NDEBUG
        printf("after remove unbalanced triangle: n = %u\n", new_n);
#endif
        // 2. 对于剩下的图，计算无符号图的heu_plex
        // 做degeneracy peeling
        for (auto u : in) degree[u] = in.intersect(matrix[u]);

        vector<ui> peel_sequence;
        ui idx = new_n;
        for (ui i = 0; i < new_n; i++) {
            ui u, min_degree = new_n;
            auto cand = in & vis;
            for (auto j : cand) if (degree[j] < min_degree) {
                min_degree = degree[j];
                u = j;
            }

            peel_sequence.push_back(u);
            vis.remove(u);

            if (idx == new_n && min_degree + K >= new_n - i) idx = i;

            cand &= matrix[u];
            for (auto j : cand) degree[j]--;
        }

        // 目前得到的启发式解是degeneracy中idx到new_n-1的后缀
        Set heu_solution(n);

        for (ui i = idx; i < new_n; i++) heu_solution.insert(peel_sequence[i]);
#ifndef NDEBUG
        printf("\tdegeneracy order find a kplex of size %u\n", heu_solution.size());
        for (ui u : heu_solution) {
            ui degree_in_heu = heu_solution.intersect(matrix[u]);
            assert(degree_in_heu + K >= heu_solution.size());
        }
#endif
        // 从idx-1开始，依次减到0，尝试将其加入到heu_solution中
        for (ui cur = 0; cur + 1 <= idx; cur++) {
            ui u = peel_sequence[idx - 1 - cur];
            bool can_add_to_heu = true;

            heu_solution.insert(u);
            for (auto v : heu_solution) {
                ui degree_in_heu = heu_solution.intersect(matrix[v]);
                if (degree_in_heu + K < heu_solution.size()) {
                    can_add_to_heu = false;
                    break;
                }
            }

            if (!can_add_to_heu) heu_solution.remove(u);
        }
#ifndef NDEBUG
        printf("\tadd vertices to heu_solution find a kplex of size %u\n", heu_solution.size());
        for (ui u : heu_solution) {
            ui degree_in_heu = heu_solution.intersect(matrix[u]);
            assert(degree_in_heu + K >= heu_solution.size());
        }
#endif

        if (heu_solution.size() > best_solution_size) {
            best_solution_size = heu_solution.size();
            best_solution = heu_solution;
#ifndef NDEBUG
            printf("\tubsigned: find a heu solution of size %u\n", best_solution_size);
#endif
        }
    }

    void heu_kPlex_in_signed_graph()
    {
    }

    void init_bnb(Set &S, Set &C, int choose_u)
    {
        // initialize S and C
        C.flip();
        addedVertex = choose_u;
        if (choose_u != -1) {
            S.insert(choose_u);
            C.remove(choose_u);
        }
    }

    void bnb_search_two_hop(Set &S, Set &C)
    {
        // ui level = S.size();
        bool pruned = false;
        update_SC(S, C, pruned);
        if (pruned) return;

        ui ub = upper_bound_two_hop(S, C);
        if (ub <= best_solution_size) return;
#ifndef NDEBUG
        // 确定C中每个点都能加入到S中
        for (ui v : C) {
            if (S.intersect(matrix[v]) + K <= S.size()) {
                printf("v = %d, S.size = %d, K = %d, intersect_sum = %d\n",
                    v, S.size(), K, S.intersect(matrix[v]));
                assert(false);
            }
        }
#endif

        //         if (ori_u != -1 && C.intersect(non_matrix[ori_u]) == 0) {
        // #ifndef NDEBUG
        //             printf("ori_u = %d, C.size = %d, non_matrix[ori_u].size = %d\n",
        //                 ori_u, C.size(), non_matrix[ori_u].size());
        // #endif
        //             Set S0 = S, SL(n), SR(n);
        //             Set CL = C & matrix[ori_u], CR = C & matrix[ori_u];
        //             bnb_search_one_hop(S0, SL, SR, CL, CR);
        //             return;
        //         }
        // choose branching vertex
        ui u = choose_branch_vertex_two_hop(S, C);
        assert(C.contains(u));
        // generate two branches
        // the first branch includes u into S
        {
            auto new_S = S, new_C = C;
            new_S.insert(u);
            addedVertex = u;
            new_C.remove(u);
            bnb_search_two_hop(new_S, new_C);
        }
        // the second branch excludes u from S
        {
            auto new_S = S, new_C = C;
            new_C.remove(u);
            addedVertex = -1;
            bnb_search_two_hop(new_S, new_C);
        }
    }
    /**
     * 当待选集中只剩下初始点u的邻居时，调用bnb_search_one_hop
     */
    void bnb_search_one_hop(Set &S0, Set &SL, Set &SR, Set &CL, Set &CR)
    {

    }
    /**
     * @brief reduce the candidate set C
     * @param S: the current set S
     * @param C: the current candidate set C
     * @return the reduced candidate set C
     *
     * 判断S本身是否是signed k-plex
     * 保证C中的点都能直接加入到S中
     *
     */
    void update_SC(Set &S, Set &C, bool &pruned)
    {
        pruned = false;
        for (ui u : S) {
            degree_in_S[u] = S.intersect(matrix[u]);
            degree[u] = degree_in_S[u] + C.intersect(matrix[u]);
        }
        for (ui u : C) {
            degree_in_S[u] = S.intersect(matrix[u]);
            degree[u] = degree_in_S[u] + C.intersect(matrix[u]);
        }
        // S中有新增的点，判断S是否是signed k-plex，并更新degree_in_S
        if (addedVertex != -1) {
            assert(addedVertex < n);
            for (ui u : S) {
                if (degree_in_S[u] + K < S.size()) {
                    pruned = true;
                    return;
                }
                else if (degree_in_S[u] + K == S.size()) {
                    C &= matrix[u];
                }
            }

            // reduce C
            for (ui u : C) if (degree_in_S[u] + K <= S.size()) C.remove(u);

            // update lower bound
            if (S.size() > best_solution_size) {
                best_solution_size = S.size();
                best_solution = S;
            }

            // check balance
            // TODO 结合一下这里的计算量，是否能利用上
            Set p_neighbors_in_S = S & p_matrix[addedVertex];
            Set n_neighbors_in_S = S & n_matrix[addedVertex];
            Set p_neighbors_in_C = C & p_matrix[addedVertex];
            Set n_neighbors_in_C = C & n_matrix[addedVertex];
            Set unbalanced_vertices = Set(n);
            for (ui u : p_neighbors_in_S) {
                unbalanced_vertices |= p_neighbors_in_C & n_matrix[u];
                unbalanced_vertices |= n_neighbors_in_C & p_matrix[u];
            }
            for (ui u : n_neighbors_in_S) {
                unbalanced_vertices |= n_neighbors_in_C & n_matrix[u];
                unbalanced_vertices |= p_neighbors_in_C & p_matrix[u];
            }
            C ^= unbalanced_vertices;
        }
        // 计算每个点的上界
        queue<ui> q;
        while (!q.empty()) q.pop();
        for (ui u : S) if (degree[u] + K <= best_solution_size) {
            pruned = true;
            return;
        }
        for (ui u : C) if (degree[u] + K <= best_solution_size) {
            q.push(u);
        }

        while (!q.empty()) {
            ui u = q.front(); q.pop();
            assert(C.contains(u));
            C.remove(u);
            Set neighbors_in_S = (p_matrix[u] | n_matrix[u]) & S;
            Set neighbors_in_C = (p_matrix[u] | n_matrix[u]) & C;
            for (ui v : neighbors_in_S) if ((--degree[v]) + K == best_solution_size) {
                pruned = true;
                return;
            }
            for (ui v : neighbors_in_C) if ((--degree[v]) + K == best_solution_size) {
                q.push(v);
            }
        }
        // kPlexT中的上界
        // ui C_size_old = C.size();
        reduce_kPlexT(S, C);
        // printf("reduce_kPlexT END: C_size: %d -> %d\n", C_size_old, C.size());
    }

    ui upper_bound_two_hop(Set &S, Set &C)
    {
        ui ub = S.size() + C.size();
        ub = min(ub, upperbound_based_partition_two_hop(S, C));
        return ub;
    }

    ui upperbound_based_partition_two_hop(Set &S, Set &C)
    {
        // TODO：reduce pi_0
        ui ub = S.size();
        // pi0: 与S全连接的C中的点
        Set S_copy = S;
        Set C_copy = C;

        // 计算S中每个顶点在C中的非邻居数量和可缺失的边数
        vector<ui> &cost = array_n_1; // S中每个顶点还能支持的缺失边数
        for (ui v : S_copy) {
            assert(K >= S_copy.size() - degree_in_S[v]);
            cost[v] = K - (S_copy.size() - degree_in_S[v]);
        }

        while (S_copy.size()) {
            int u = -1;
            int u_pi = 0, u_cost = 0;

            for (ui v : S_copy) {
                int v_pi = C_copy.intersect(non_matrix[v]);
                int v_cost = cost[v];
                if (v_pi <= v_cost) {
                    S_copy.remove(v);
                    continue;
                }
                if (v_cost == 0) {
                    u = v;
                    break;
                }
                // 计算v在C中的缺失边数 v_pi
                // u_pi / u_cost --- v_pi / cost[v]
                else if (u == -1 || u_pi * v_cost < v_pi * u_cost) {
                    u_pi = v_pi;
                    u_cost = v_cost;
                    u = v;
                }
                else if (u_pi * v_cost == v_pi * u_cost) {
                    if (v_cost < u_cost) {
                        u_pi = v_pi;
                        u_cost = v_cost;
                        u = v;
                    }
                }
            }
            if (u == -1) break;
            if (u_pi <= u_cost) {
                printf("u_pi = %d, u_cost = %d\n", u_pi, u_cost);
                fflush(stdout);
            }
            assert(u_pi > u_cost);

            S_copy.remove(u);
            ub += u_cost;
            C_copy = C_copy & matrix[u];
        }

        // C中剩余的点是pi_0
        ub += C_copy.size();

        return ub;
    }
    // 计算C中每个点v的上界
    void reduce_kPlexT(Set &S, Set &C)
    {
        if (S.size() <= 1) return;

        Set S2(n);
        Set SC = S | C;
        ui SC_size = S.size() + C.size();

        // 检查S中的点
        for (ui u : S) if (SC_size > degree[u] + K) S2.insert(u);
        if (S2.size() == 0) return;

        int S2_support = 0; // S2中所有点还能增加的缺失边数
        for (ui u : S2) S2_support += K - (S.size() - degree_in_S[u]);

        vector<pair<ui, ui>> vertex_loss_sorted; // (loss_S2_v, v)
        if (S2.size() == S.size()) for (ui v : C) {
            ui loss_S2_v = S.size() - degree_in_S[v];
            vertex_loss_sorted.push_back({ loss_S2_v, v });
        }
        else for (ui v : C) {  // S2.size() < S.size()
            // 计算C中点v在S2中的缺失边
            ui loss_S2_v = S2.intersect(non_matrix[v]);
            vertex_loss_sorted.push_back({ loss_S2_v, v });
        }

        sort(vertex_loss_sorted.begin(), vertex_loss_sorted.end());

        // 计算C中每个点v的上界
        for (ui u : C) {
            ui u_ub = S.size() + 1; // add u to S
            int u_sup = S2_support - S2.intersect(non_matrix[u]);
            // u在S中最多有K-1个非邻居
            ui non_neighbor_of_u_in_S = S.size() - degree_in_S[u];
            // 贪心地加入缺失边最少的邻居
            for (auto p : vertex_loss_sorted) {
                ui v = p.second;
                ui loss_S2_v = p.first;
                if (u_ub > best_solution_size) break;
                if (loss_S2_v > u_sup) break;
                if (u == v) continue;
                if (!C.contains(v)) continue;
                if (non_matrix[u][v]) {
                    if (non_neighbor_of_u_in_S == K - 1) continue;
                    non_neighbor_of_u_in_S++;
                }
                u_sup -= loss_S2_v;
                u_ub++;
            }

            if (u_ub <= best_solution_size) C.remove(u);
        }
    }

    ui choose_branch_vertex_two_hop(Set &S, Set &C)
    {
        ui u = n;
        ui min_degree_in_S = n;

        Set cand = C;
        if (ori_u != -1 && cand.intersect(non_matrix[ori_u]) > 0) cand &= non_matrix[ori_u];

        for (ui v : cand) if (degree_in_S[v] < min_degree_in_S) {
            u = v;
            min_degree_in_S = degree_in_S[v];
        }
        assert(u != n);
        return u;
    }

    // 输出set中的每个数
    void print_set(Set &S, string name)
    {
        cout << name << ": size = " << S.size() << "; ";
        for (ui v : S) printf("%u ", v);
        printf("\n");
    }
};

#endif /* _SIGNED_KPLEX_BITSET_ */