/*
This file contains code from the Maximum-kPlex project, which is licensed under the MIT License.
The original code and license can be found at: https://github.com/LijunChang/Maximum-kPlex
*/

#ifndef _SIGNED_KPLEX_
#define _SIGNED_KPLEX_

#include "Utility.h"
#include "Timer.h"

#define _SECOND_ORDER_PRUNING_
class SIGNED_KPLEX
{
private:
    ui K;
    ui n;
    ept m;
    int *matrix;
    long long matrix_size;

#ifdef _SECOND_ORDER_PRUNING_
    ui *cn;
    std::queue<Edge> Qe;
    std::vector<Edge> removed_edges;
    long long removed_edges_n;
#endif

    ui best_solution_size;
    ui *best_solution;
    ui *degree;
    ui *degree_in_S;
    ui *neighbors;
    ui *nonneighbors;
    ui *SC;
    ui *SC_rid;
    ui *level_id;
    ui *vis;
    queue<ui> Qv;

    int dfs_cnt = 0;
    int dfs_cnt_after_ub = 0;
    vector<int> level_cnt;
    vector<int> level_cnt_after_ub;

public:
    SIGNED_KPLEX()
    {
        n = m = 0;
        matrix = NULL;
        matrix_size = 0;
        best_solution_size = 0;
        best_solution = NULL;
#ifdef _SECOND_ORDER_PRUNING_
        cn = NULL;
        removed_edges_n = 0;
#endif
        degree = NULL;
        degree_in_S = NULL;
        neighbors = NULL;
        nonneighbors = NULL;
        SC = NULL;
        SC_rid = NULL;
        level_id = NULL;
        vis = NULL;
    }

    ~SIGNED_KPLEX()
    {
        if (matrix != NULL) {
            delete[] matrix;
            matrix = NULL;
        }
        if (degree != NULL) {
            delete[] degree;
            degree = NULL;
        }
#ifdef _SECOND_ORDER_PRUNING_
        if (cn != NULL) {
            delete[] cn;
            cn = NULL;
        }
#endif
        if (best_solution != NULL) {
            delete[] best_solution;
            best_solution = NULL;
        }
        if (degree_in_S != NULL) {
            delete[] degree_in_S;
            degree_in_S = NULL;
        }
        if (neighbors != NULL) {
            delete[] neighbors;
            neighbors = NULL;
        }
        if (nonneighbors != NULL) {
            delete[] nonneighbors;
            nonneighbors = NULL;
        }
        if (SC != NULL) {
            delete[] SC;
            SC = NULL;
        }
        if (SC_rid != NULL) {
            delete[] SC_rid;
            SC_rid = NULL;
        }
        if (level_id != NULL) {
            delete[] level_id;
            level_id = NULL;
        }
        if (vis != NULL) {
            delete[] vis;
            vis = NULL;
        }
    }

    void allocateMemory(ui n, ui m)
    {
        matrix_size = m * 2;
        matrix = new int[matrix_size];
#ifdef _SECOND_ORDER_PRUNING_
        cn = new ui[matrix_size];
#endif
        best_solution = new ui[n];
        degree = new ui[n];
        degree_in_S = new ui[n];
        neighbors = new ui[n];
        nonneighbors = new ui[n];
        SC = new ui[n];
        SC_rid = new ui[n];
        level_id = new ui[n];
        vis = new ui[n];
        level_cnt.resize(n, 0);
        level_cnt_after_ub.resize(n, 0);
    }

    void load_graph(ui _n, const vector<Edge> &vp)
    {
        n = _n;
        if (1ll * n * n > matrix_size) {
            while (1ll * n * n > matrix_size) matrix_size *= 2;
            delete[] matrix; matrix = new int[matrix_size];
#ifdef _SECOND_ORDER_PRUNING_
            delete[] cn; cn = new ui[matrix_size];
#endif
        }

#ifdef _SECOND_ORDER_PRUNING_
        memset(cn, 0, sizeof(ui) * matrix_size);
#endif
        memset(matrix, 0, sizeof(int) * matrix_size);
        memset(degree, 0, sizeof(ui) * n);
        memset(degree_in_S, 0, sizeof(ui) * n);

        m = vp.size();
        for (ui i = 0; i < m; i++) {
            int a = vp[i].a, b = vp[i].b, c = vp[i].c;
            assert(a < n);
            assert(b < n);
            assert(c == 1 || c == -1);
            degree[a]++;
            degree[b]++;
            matrix[a * n + b] = matrix[b * n + a] = c;
        }

#ifndef NDEBUG
        printf("matrix kplex load_graph: n=%u, m=%lu\n", n, m);
#endif
    }

    void kPlex(ui K_, vector<ui> &kplex, int choose_u)
    {
        K = K_;
        best_solution_size = kplex.size();
        ui S_end = 0, C_end = 0;
        vector<ui> pivot_set;
        pivot_set.clear();
        init(S_end, C_end, choose_u, pivot_set);
        if (C_end) kplex_search(S_end, C_end, 1, pivot_set);
        if (best_solution_size > kplex.size()) {
            kplex.clear();
            for (int i = 0; i < best_solution_size; i++)
                kplex.push_back(best_solution[i]);
        }
    }

    void heu_kPlex(ui K_, vector<ui> &kplex, int choose_u)
    {
        K = K_;
        best_solution_size = kplex.size();
#ifndef NDEBUG
        printf("heu_kplex: n = %u\n", n);
#endif

        // 对于每一个不平衡三角形，选择度数最小的点，将其删去
        memset(vis, 0, sizeof(ui) * n);
        for (ui i = 0; i < n; i++) if (!vis[i]) {
            for (ui j = i + 1; j < n; j++) if (!vis[j]) {
                if (!matrix[i * n + j]) continue;
                for (ui k = j + 1; k < n; k++) if (!vis[k]) {
                    if (!matrix[i * n + k] || !matrix[j * n + k]) continue;
                    int tri_sum = matrix[i * n + j] + matrix[j * n + k] + matrix[i * n + k];
                    if (tri_sum == 1 || tri_sum == -3) {
                        ui min_degree_vertex = i;
                        if (i == 0) min_degree_vertex = degree[j] < degree[k] ? j : k;
                        else {
                            if (degree[j] < degree[min_degree_vertex]) min_degree_vertex = j;
                            if (degree[k] < degree[min_degree_vertex]) min_degree_vertex = k;
                        }
                        vis[min_degree_vertex] = 1;
                        for (ui v = 0; v < n; v++) if (matrix[min_degree_vertex * n + v]) degree[v]--;
                    }
                }
            }
        }
        assert(choose_u == 0);
        assert(vis[choose_u] == 0);

        vector<ui> mapping;
        ui N = n;
        for (ui i = 0; i < N; i++) if (!vis[i]) mapping.push_back(i);
        n = mapping.size();

        // 重新计算matrix
        if (n < N) {
            for (ui i = 0; i < n; i++) for (ui j = i + 1; j < n; j++) matrix[i * n + j] = matrix[mapping[i] * N + mapping[j]];
            for (ui i = 0; i < n; i++) for (ui j = i + 1; j < n; j++) matrix[j * n + i] = matrix[i * n + j];
            for (ui i = 0; i < n; i++) matrix[i * n + i] = 0;
            // 重新计算每个点的度数
            for (ui i = 0; i < n; i++) {
                vis[i] = degree[i] = 0;
                for (ui j = 0; j < n; j++) if (matrix[i * n + j]) degree[i]++, m++;
            }
        }
#ifndef NDEBUG
        printf("after remove unbalanced triangle: n = %u, m = %lu\n", n, m);
#endif
        // 做无符号图的启发式kplex
        ui *peel_sequence = neighbors;
        ui *core = nonneighbors;
        ui max_core = 0, UB = 0, idx = n;
        memset(vis, 0, sizeof(ui) * n);
        for (ui i = 0; i < n; i++) {
            ui u, min_degree = n;
            for (ui j = 0; j < n; j++) if (!vis[j] && degree[j] < min_degree) {
                u = j;
                min_degree = degree[j];
            }
            if (min_degree > max_core) max_core = min_degree;
            core[u] = max_core;
            peel_sequence[i] = u;
            vis[u] = 1;

            ui t_UB = core[u] + K;
            if (n - i < t_UB) t_UB = n - i;
            if (t_UB > UB) UB = t_UB;

            if (idx == n && min_degree + K >= n - i) idx = i;

            for (ui j = 0; j < n; j++) if (!vis[j] && matrix[u * n + j]) --degree[j];
        }

        // 目前得到的启发式解是idx到n-1的一个后缀
        vector<ui> heu_solution_in_map;
        heu_solution_in_map.clear();
        vector<ui> degree_in_heu;
        degree_in_heu.clear();

        for (ui i = idx; i < n; i++) heu_solution_in_map.push_back(peel_sequence[i]);
        for (ui i = 0; i < heu_solution_in_map.size(); i++) {
            degree_in_heu.push_back(0);
            ui u = heu_solution_in_map[i];
            for (ui j = 0; j < heu_solution_in_map.size(); j++) {
                ui v = heu_solution_in_map[j];
                if (matrix[u * n + v]) degree_in_heu[i]++;
            }
        }

        // 从idx-1开始，依次减到0，尝试将其加入到S中
        for (ui cur = 0; cur + 1 <= idx; cur++) {
            ui u = peel_sequence[idx - 1 - cur];
            bool can_add_to_heu = true;

            vector<ui> neighbors_vec, nonneighbors_vec;
            for (ui i = 0; i < heu_solution_in_map.size(); i++) {
                if (matrix[u * n + heu_solution_in_map[i]]) neighbors_vec.push_back(i);
                else nonneighbors_vec.push_back(i);
            }

            for (ui i = 0; i < neighbors_vec.size(); i++) {
                ui v = heu_solution_in_map[neighbors_vec[i]];
                assert(degree_in_heu[neighbors_vec[i]] + K >= heu_solution_in_map.size());
                if (degree_in_heu[neighbors_vec[i]] + K == heu_solution_in_map.size()) {
                    can_add_to_heu = false;
                    break;
                }
            }
            if (neighbors_vec.size() + K <= heu_solution_in_map.size()) can_add_to_heu = false;

            if (can_add_to_heu) {
                for (ui i = 0; i < neighbors_vec.size(); i++)  degree_in_heu[neighbors_vec[i]]++;
                heu_solution_in_map.push_back(u);
                degree_in_heu.push_back(neighbors_vec.size());
            }
        }

        if (heu_solution_in_map.size() > best_solution_size) {
            best_solution_size = heu_solution_in_map.size();
            for (int i = 0; i < best_solution_size; i++) best_solution[i] = mapping[heu_solution_in_map[i]];
            printf("Degen find a solution of size %u\n", best_solution_size);
        }
        if (best_solution_size > kplex.size()) {
            kplex.clear();
            for (int i = 0; i < best_solution_size; i++) kplex.push_back(best_solution[i]);
        }
    }

private:
    // init S, C
    void init(ui &S_end, ui &C_end, int choose_u, vector<ui> &pivot_set)
    {
        S_end = 0, C_end = 0;
        // k-core
        queue<ui> q;
        memset(vis, 0, sizeof(ui) * n);
        for (ui i = 0; i < n; i++) if (degree[i] + K <= best_solution_size) q.push(i);
        while (!q.empty()) {
            ui u = q.front(); q.pop();
            if (vis[u]) continue; vis[u] = 1;
            for (ui v = 0; v < n; v++) if (matrix[u * n + v]) {
                if ((degree[v]--) + K == best_solution_size + 1) q.push(v);
            }
        }

        if (choose_u != -1 && vis[choose_u]) return;

        for (ui i = 0; i < n; i++) SC_rid[i] = n;

        for (ui i = 0; i < n; i++) if (!vis[i]) {
            SC[C_end] = i;
            SC_rid[i] = C_end++;
        }

        memset(degree, 0, sizeof(ui) * n);
        for (ui i = 0; i < C_end; i++) {
            ui u = SC[i];
            assert(SC_rid[u] == i);
            for (ui j = 0; j < C_end; j++) if (matrix[u * n + SC[j]]) degree[u]++;
        }

        memset(level_id, 0, sizeof(ui) * n);
        for (ui i = 0; i < C_end; i++) level_id[SC[i]] = n;

        while (!Qv.empty()) Qv.pop();
#ifdef _SECOND_ORDER_PRUNING_
        while (!Qe.empty()) Qe.pop();
#endif

#ifdef _SECOND_ORDER_PRUNING_
        for (ui i = 0; i < C_end; i++) {
            ui neighbors_n = 0;
            int *t_matrix = matrix + SC[i] * n;
            for (ui j = 0; j < C_end; j++) if (t_matrix[SC[j]]) neighbors[neighbors_n++] = SC[j];
            for (ui j = 0; j < neighbors_n; j++) for (ui k = j + 1; k < neighbors_n; k++) {
                ++cn[neighbors[j] * n + neighbors[k]];
                ++cn[neighbors[k] * n + neighbors[j]];
            }
        }

        assert(Qe.empty());
        for (ui i = 0; i < C_end; i++) for (ui j = i + 1; j < C_end; j++) {
            if (matrix[SC[i] * n + SC[j]] && upper_bound_based_prune(S_end, SC[i], SC[j])) {
                Qe.push(Edge(SC[i], SC[j], matrix[SC[i] * n + SC[j]]));
            }
        }
        removed_edges_n = 0;
#endif
        if (remove_vertices(S_end, C_end, 0)) {
            C_end = 0;
            return;
        }
        if (SC_rid[choose_u] >= C_end) {
            C_end = 0;
            return;
        }
        // 将choose_u加入S
        if (choose_u != -1) {
            bool pruned = moveu_C_to_S(S_end, C_end, 0, choose_u);
            if (pruned) {
                S_end = 0, C_end = 0;
                return;
            }
            assert(choose_u == SC[0]);
            assert(SC_rid[choose_u] == 0);
            assert(S_end == 1);
            update_pivot_set(S_end, C_end, choose_u, pivot_set);
        }
    }

    void kplex_search(ui S_end, ui C_end, ui level, vector<ui> pivot_set)
    {
        // printf("S_end = %d, C_end = %d, level = %d\n", S_end, C_end, level);
        if (S_end > best_solution_size) {
            best_solution_size = S_end;
            for (ui i = 0; i < best_solution_size; i++) best_solution[i] = SC[i];
#ifndef NDEBUG
            // 每个点在S中的缺失边少于K
            for (ui i = 0; i < best_solution_size; i++)
                assert(degree_in_S[best_solution[i]] + K >= best_solution_size);
            printf("Find a k-plex of size: %u\n", best_solution_size);
#endif
        }
        if (C_end <= best_solution_size) return;

#ifndef NDEBUG
#ifdef _SECOND_ORDER_PRUNING_
        for (ui i = 0; i < C_end; i++) for (ui j = i + 1; j < C_end; j++) {
            ui v = SC[i], w = SC[j];
            ui common_neighbors = 0;
            for (ui k = S_end; k < C_end; k++) if (matrix[SC[k] * n + v] && matrix[SC[k] * n + w]) ++common_neighbors;
            assert(cn[v * n + w] == common_neighbors);
            assert(cn[w * n + v] == common_neighbors);
        }
#endif
        for (ui i = 0; i < C_end; i++) assert(degree[SC[i]] + K > best_solution_size);
        for (ui i = 0; i < S_end; i++) assert(degree_in_S[SC[i]] + K >= S_end);
        for (ui i = S_end; i < C_end; i++) assert(degree_in_S[SC[i]] + K > S_end);
        for (ui i = 0; i < C_end; i++) assert(level_id[SC[i]] == n);
        for (ui i = 0; i < C_end; i++) assert(check_balance(S_end, SC[i]));
#endif

        dfs_cnt++;
        level_cnt[level]++;

        // upper bound
        ui ub = upper_bound(S_end, C_end);
        if (ub <= best_solution_size) return;

        dfs_cnt_after_ub++;
        level_cnt_after_ub[level]++;

        ui old_kplex_size = best_solution_size, old_C_end = C_end;
        ui old_removed_edges_n = 0;
#ifdef  _SECOND_ORDER_PRUNING_
        old_removed_edges_n = removed_edges_n;
#endif

        // choose branching vertex
        ui u = choose_branch_vertex_with_min_degree(S_end, C_end, pivot_set);
        assert(u != n);
        assert(SC[SC_rid[u]] == u && SC_rid[u] >= S_end && SC_rid[u] < C_end);
        assert(degree[u] + K > best_solution_size);
        assert(check_balance(S_end, u));
        assert(Qv.empty());
#ifdef _SECOND_ORDER_PRUNING_
        assert(Qe.empty());
#endif
        // the first branch includes u into S
        // printf("level = %d, u: C->S, u = %d\n", level, u);
        vector<ui> new_pivot_set;
        bool pruned = moveu_C_to_S(S_end, C_end, level, u);
        // 加入点u时，向kplex_search函数传入空的pivot_set
        if (!pruned) {
            new_pivot_set.clear();
            kplex_search(S_end, C_end, level + 1, new_pivot_set);
        }

        restore_C(S_end, C_end, old_C_end, old_removed_edges_n, level);
        moveu_S_to_C(S_end, C_end, level);
        assert(C_end == old_C_end);
        assert(u == SC[S_end]);
#ifdef  _SECOND_ORDER_PRUNING_
        assert(removed_edges_n == old_removed_edges_n);
#endif

        // the second branch exclude u from S
        // printf("level = %d, u: C->X, u = %d\n", level, u);
        assert(Qv.empty());
#ifdef _SECOND_ORDER_PRUNING_
        assert(Qe.empty());
#endif
        pruned = false;
        if (best_solution_size > old_kplex_size) pruned = reduce_SC_based_lb(S_end, C_end, level);
        // ub = upper_bound(S_end, C_end);
        // if (!pruned) pruned = (ub <= best_solution_size);
        if (!pruned) pruned = moveu_C_to_X(S_end, C_end, u, level);
        if (!pruned) {
            new_pivot_set.clear();
            update_pivot_set(S_end, C_end, u, new_pivot_set);
            kplex_search(S_end, C_end, level + 1, new_pivot_set);
        }
        restore_C(S_end, C_end, old_C_end, old_removed_edges_n, level);
#ifdef  _SECOND_ORDER_PRUNING_
        assert(removed_edges_n == old_removed_edges_n);
#endif
    }
    /**
     * 选择分支顶点
     *
     * 选择策略是找到在集合S中度数最小的顶点
     *
     * 20250331：加入pivot点
     */
    ui choose_branch_vertex_with_min_degree(ui S_end, ui C_end, vector<ui> pivot_set)
    {
        // 找到degree_in_S最小的点
        ui u = n;
        ui min_degree_in_S = n;

        if (!pivot_set.empty()) {
            for (ui i = 0; i < pivot_set.size(); i++) {
                ui v = pivot_set[i];
                assert(SC_rid[v] >= S_end && SC_rid[v] < C_end);
                assert(SC[SC_rid[v]] == v);
                if (degree_in_S[v] < min_degree_in_S) {
                    u = v;
                    min_degree_in_S = degree_in_S[v];
                }
            }
        }
        else {
            for (ui i = S_end; i < C_end; i++) {
                ui v = SC[i];
                if (degree_in_S[v] < min_degree_in_S) {
                    u = v;
                    min_degree_in_S = degree_in_S[v];
                }
            }
        }
        assert(u != n);
        return u;
    }
    /**
     * 将顶点u从候选集C移动到S中，并更新集合SC
     *
     * 主要步骤:
     * 1. 将u移到S集合末尾
     * 2. 更新u的邻居在S中的度数
     * 3. 根据k-plex约束进行剪枝
     * 4. 检查平衡性约束
     * 5. SECOND_ORDER_PRUNING
     */
    bool moveu_C_to_S(ui &S_end, ui &C_end, ui level, ui u)
    {
        assert(Qv.empty());
#ifdef _SECOND_ORDER_PRUNING_
        assert(Qe.empty());
#endif
        // 1
        swap_pos(SC_rid[u], S_end++);
        assert(u == SC[S_end - 1]);

        // 2
        ui neighbors_n = 0, nonneighbors_n = 0;
        int *t_matrix = matrix + u * n;
        for (ui i = 0; i < C_end; i++) if (SC[i] != u) {
            if (t_matrix[SC[i]]) neighbors[neighbors_n++] = SC[i];
            else nonneighbors[nonneighbors_n++] = SC[i];
        }

        for (ui i = 0; i < neighbors_n; i++) degree_in_S[neighbors[i]]++;

        // 3
        if (degree_in_S[u] + K == S_end) {
            ui i = 0;
            while (i < nonneighbors_n && SC_rid[nonneighbors[i]] < S_end) i++;
            for (; i < nonneighbors_n; i++) if (level_id[nonneighbors[i]] == n) {
                level_id[nonneighbors[i]] = level;
                Qv.push(nonneighbors[i]);
            }
        }
        else {
            ui i = 0;
            while (i < nonneighbors_n && SC_rid[nonneighbors[i]] < S_end) i++;
            for (; i < nonneighbors_n; i++) if (level_id[nonneighbors[i]] == n) {
                if (degree_in_S[nonneighbors[i]] + K <= S_end) {
                    level_id[nonneighbors[i]] = level;
                    Qv.push(nonneighbors[i]);
                }
            }
        }
        for (ui i = 0; i < nonneighbors_n && SC_rid[nonneighbors[i]] < S_end; i++) {
            if (degree_in_S[nonneighbors[i]] + K == S_end) {
                int *tt_matrix = matrix + nonneighbors[i] * n;
                for (ui j = S_end; j < C_end; j++) {
                    if (level_id[SC[j]] == n && !tt_matrix[SC[j]]) {
                        level_id[SC[j]] = level;
                        Qv.push(SC[j]);
                    }
                }
            }
        }

        // 4.check balance
        {
            ui i = 0;
            while (i < neighbors_n && SC_rid[neighbors[i]] < S_end) i++;
            for (; i < neighbors_n; i++) {
                ui v = neighbors[i];
                if (level_id[neighbors[i]] == n && !check_balance(S_end, neighbors[i])) {
                    level_id[neighbors[i]] = level;
                    Qv.push(neighbors[i]);
                }
            }
        }

#ifdef _SECOND_ORDER_PRUNING_
        // update cn(.,.)
        for (ui i = 0; i < neighbors_n; i++) { // process common neighbors of u
            for (ui j = i + 1; j < neighbors_n; j++) {
                ui v = neighbors[i], w = neighbors[j];
                assert(cn[v * n + w]);
                --cn[v * n + w];
                --cn[w * n + v];
            }
        }

        // u和u的非邻居，组成的节点对，其upper_bound_based_prune会减小
        for (ui i = 0; i < nonneighbors_n; i++) {
            int v = nonneighbors[i];
            assert(!t_matrix[v]);
            if (SC_rid[v] < S_end || level_id[v] == level || t_matrix[v]) continue;
            if (upper_bound_based_prune(S_end, u, v)) {
                level_id[v] = level;
                Qv.push(v);
            }
        }

        // 如果有一个是u的邻居，upper_bound_based_prune不变
        int new_n = 0;
        for (ui i = 0; i < nonneighbors_n; i++) if (level_id[nonneighbors[i]] > level) nonneighbors[new_n++] = nonneighbors[i];
        nonneighbors_n = new_n;
        for (ui i = 1; i < nonneighbors_n; i++) { // process common non-neighbors of u
            ui w = nonneighbors[i];
            for (ui j = 0; j < i; j++) {
                ui v = nonneighbors[j];
                if (!upper_bound_based_prune(S_end, v, w)) continue;
                if (SC_rid[w] < S_end) return true; // v, w \in S --- UB2
                else if (SC_rid[v] >= S_end) { // v, w, \in R --- RR5
                    if (matrix[v * n + w]) Qe.push(Edge(v, w, matrix[v * n + w]));
                }
                else { // RR4
                    assert(level_id[w] > level);
                    level_id[w] = level;
                    Qv.push(w);
                    break;
                }
            }
        }
#endif
        return remove_vertices(S_end, C_end, level);
    }
    /**
     * 从候选集C中移除顶点和边
     */
    bool remove_vertices(ui S_end, ui &C_end, ui level)
    {
#ifdef _SECOND_ORDER_PRUNING_
        while (!Qv.empty() || !Qe.empty())
#else
        while (!Qv.empty())
#endif
        {
            while (!Qv.empty()) {
                // remove u from C
                ui u = Qv.front(); Qv.pop();
                // printf("delete qv: u = %d\n", u);
                assert(level_id[u] == level);
                assert(SC[SC_rid[u]] == u);
                assert(SC_rid[u] >= S_end && SC_rid[u] < C_end);
                swap_pos(SC_rid[u], --C_end);
                assert(u == SC[C_end]);

                int *t_matrix = matrix + u * n;
                bool terminate = false;

                ui neighbors_n = 0;
                for (ui i = 0; i < C_end; i++) if (t_matrix[SC[i]]) {
                    ui v = SC[i];
                    neighbors[neighbors_n++] = v;
                    degree[v]--;
                    if (degree[v] + K <= best_solution_size) {
                        if (i < S_end) terminate = true;
                        else if (level_id[v] == n) {
                            level_id[v] = level;
                            Qv.push(v);
                        }
                    }
                }
#ifdef _SECOND_ORDER_PRUNING_
                for (ui i = 1; i < neighbors_n; i++) {
                    ui w = neighbors[i];
                    for (ui j = 0; j < i; j++) {
                        ui v = neighbors[j];
                        assert(cn[v * n + w]);
                        --cn[v * n + w];
                        --cn[w * n + v];

                        if (!upper_bound_based_prune(S_end, v, w)) continue;

                        if (SC_rid[w] < S_end) terminate = true; // v < w < S_end
                        else if (SC_rid[v] >= S_end) { // S_end < v < w
                            if (matrix[v * n + w]) Qe.push(Edge(v, w, matrix[v * n + w]));
                        }
                        else if (level_id[w] == n) { // v < S_end < w
                            level_id[w] = level;
                            Qv.push(w);
                        }
                    }
                }
#endif
                if (terminate) return true;
            }
#ifdef _SECOND_ORDER_PRUNING_
            // 从Qe中取出一条边，并删除
            if (Qe.empty()) break;
            ui v = Qe.front().a, w = Qe.front().b;
            int sign = Qe.front().c; Qe.pop();
            // printf("delete qe: %d, %d, %d\n", v, w, sign);
            if (level_id[v] <= level || level_id[w] <= level || !matrix[v * n + w]) continue;
            assert(SC_rid[v] >= S_end && SC_rid[v] < C_end && SC_rid[w] >= S_end && SC_rid[w] < C_end);

            assert(level_id[v] == n && level_id[w] == n);
            if (degree[v] + K - 1 <= best_solution_size) {
                level_id[v] = level;
                Qv.push(v);
            }
            if (degree[w] + K - 1 <= best_solution_size) {
                level_id[w] = level;
                Qv.push(w);
            }
            if (!Qv.empty()) continue;

            assert(matrix[v * n + w] == sign);
            matrix[v * n + w] = matrix[w * n + v] = 0;
            --degree[v]; --degree[w];

            if (removed_edges.size() == removed_edges_n) {
                removed_edges.push_back(Edge(v, w, sign));
                ++removed_edges_n;
            }
            else removed_edges[removed_edges_n++] = Edge(v, w, sign);

            int *t_matrix = matrix + v * n;
            for (ui i = 0; i < C_end; i++) if (t_matrix[SC[i]]) {
                --cn[w * n + SC[i]];
                --cn[SC[i] * n + w];
                if (!upper_bound_based_prune(S_end, w, SC[i])) continue;
                if (i < S_end) {
                    if (level_id[w] == n) {
                        level_id[w] = level;
                        Qv.push(w);
                    }
                }
                else if (matrix[w * n + SC[i]]) Qe.push(Edge(w, SC[i], matrix[w * n + SC[i]]));
            }
            t_matrix = matrix + w * n;
            for (ui i = 0; i < C_end; i++) if (t_matrix[SC[i]]) {
                --cn[v * n + SC[i]];
                --cn[SC[i] * n + v];
                if (!upper_bound_based_prune(S_end, v, SC[i])) continue;
                if (i < S_end) {
                    if (level_id[v] == n) {
                        level_id[v] = level;
                        Qv.push(v);
                    }
                }
                else if (matrix[v * n + SC[i]]) Qe.push(Edge(v, SC[i], matrix[v * n + SC[i]]));
            }
#endif
        }
        return false;
    }
    /**
     * 恢复C集合中的顶点和边
     */
    void restore_C(ui S_end, ui &C_end, ui old_C_end, ui old_removed_edges_n, ui level)
    {
        while (!Qv.empty()) {
            ui u = Qv.front(); Qv.pop();
            assert(level_id[u] == level);
            assert(SC_rid[u] < C_end);
            level_id[u] = n;
        }
#ifdef _SECOND_ORDER_PRUNING_
        while (!Qe.empty()) Qe.pop();
#endif
        // 恢复C集合中的顶点
        while (C_end < old_C_end) {
            ui u = SC[C_end];
            assert(level_id[u] == level && SC_rid[u] == C_end);
            level_id[u] = n;

            ui neighbors_n = 0;
            int *t_matrix = matrix + u * n;
            for (ui i = 0; i < C_end; i++) if (t_matrix[SC[i]]) {
                neighbors[neighbors_n++] = SC[i];
                degree[SC[i]]++;
            }
#ifdef _SECOND_ORDER_PRUNING_
            // update cn(.,.)
            for (ui i = 0; i < neighbors_n; i++) {
                ui v = neighbors[i];
                for (ui j = i + 1; j < neighbors_n; j++) {
                    ui w = neighbors[j];
                    cn[v * n + w]++;
                    cn[w * n + v]++;
                }
            }
            for (ui i = 0; i < C_end; i++) cn[u * n + SC[i]] = 0;
            for (ui i = 0; i < neighbors_n; i++) if (SC_rid[neighbors[i]] >= S_end) {
                ui v = neighbors[i];
                for (ui j = 0; j < C_end; j++) if (matrix[v * n + SC[j]]) cn[u * n + SC[j]]++;
            }
            for (ui i = 0; i < C_end; i++) cn[SC[i] * n + u] = cn[u * n + SC[i]];
#ifndef NDEBUG
            for (ui i = 0; i < C_end; i++) {
                ui common_neighbors = 0, v = SC[i], w = u;
                for (ui k = S_end; k < C_end; k++) if (matrix[SC[k] * n + v] && matrix[SC[k] * n + w]) ++common_neighbors;
                if (cn[u * n + SC[i]] != common_neighbors) printf("cn[u * n + SC[i]] = %u, comon_neighbors = %u\n", cn[u * n + SC[i]], common_neighbors);
                assert(cn[u * n + SC[i]] == common_neighbors);
            }
#endif
#endif
            C_end++;
        }
#ifdef _SECOND_ORDER_PRUNING_
        // 恢复删掉的边
        // printf("old_removed_edges_n = %lu, removed_edges_n = %lu\n", old_removed_edges_n, removed_edges_n);
        for (ui i = old_removed_edges_n; i < removed_edges_n; i++) {
            ui v = removed_edges[i].a, w = removed_edges[i].b;
            int sign = removed_edges[i].c;
            assert(SC_rid[v] >= S_end && SC_rid[v] < C_end && SC_rid[w] >= S_end && SC_rid[w] < C_end);
            assert(!matrix[v * n + w] && !matrix[w * n + v]);
            // if (matrix[v * n + w]) continue;

            matrix[v * n + w] = matrix[w * n + v] = sign;
            ++degree[v]; ++degree[w];

            for (ui i = 0; i < C_end; i++) {
                if (matrix[v * n + SC[i]]) {
                    ++cn[w * n + SC[i]];
                    ++cn[SC[i] * n + w];
                }
                if (matrix[w * n + SC[i]]) {
                    ++cn[v * n + SC[i]];
                    ++cn[SC[i] * n + v];
                }
            }
        }
        removed_edges_n = old_removed_edges_n;
#endif
    }
    /**
     * 将顶点u从集合S移动到集合C中
     *
     * 1. 将S末尾的顶点u移出S集合
     * 2. 更新与u相邻的顶点在S中的度数
     */
    void moveu_S_to_C(ui &S_end, ui C_end, ui level)
    {
        assert(S_end);
        ui u = SC[--S_end];
        int *t_matrix = matrix + u * n;
        for (ui i = 0; i < C_end; i++) if (t_matrix[SC[i]]) degree_in_S[SC[i]]--;

#ifdef _SECOND_ORDER_PRUNING_
        // update cn(.,.)
        ui neighbors_n = 0;
        for (ui i = 0; i < C_end; i++) if (t_matrix[SC[i]]) neighbors[neighbors_n++] = SC[i];
        for (ui i = 0; i < neighbors_n; i++) {
            ui v = neighbors[i];
            for (ui j = i + 1; j < neighbors_n; j++) {
                ui w = neighbors[j];
                ++cn[v * n + w];
                ++cn[w * n + v];
            }
        }
#endif
    }
    /**
     * 基于下界剪枝来化简SC集合
     */
    bool reduce_SC_based_lb(ui S_end, ui &C_end, ui level)
    {
        assert(Qv.empty());
        for (ui i = 0; i < S_end; i++) if (degree[SC[i]] + K <= best_solution_size) return true;
#ifdef _SECOND_ORDER_PRUNING_
        for (ui i = 0; i < S_end; i++) for (ui j = i + 1; j < S_end; j++) {
            if (upper_bound_based_prune(S_end, SC[i], SC[j])) return true;
        }
#endif
        for (ui i = S_end; i < C_end; i++) if (level_id[SC[i]] == n) {
            ui u = SC[i];
            if (degree[u] + K <= best_solution_size && level_id[SC[i]] == n) {
                level_id[u] = level;
                Qv.push(u);
                continue;
            }

#ifdef _SECOND_ORDER_PRUNING_
            for (ui j = 0; j < S_end; j++) {
                ui v = SC[j];
                if (degree_in_S[v] + K == S_end) assert(matrix[v * n + u]);
                if (upper_bound_based_prune(S_end, v, u) && level_id[SC[i]] == n) {
                    level_id[u] = level;
                    Qv.push(u);
                }
            }
#endif
        }
#ifdef _SECOND_ORDER_PRUNING_
        for (ui i = S_end; i < C_end; i++) if (level_id[SC[i]] == n) {
            for (ui j = i + 1; j < C_end; j++) if (level_id[SC[j]] == n && matrix[SC[i] * n + SC[j]]) {
                if (upper_bound_based_prune(S_end, SC[i], SC[j]))  Qe.push(Edge(SC[i], SC[j], matrix[SC[i] * n + SC[j]]));
            }
        }
#endif
        return remove_vertices(S_end, C_end, level);
    }
    /**
     * 将顶点u从集合C中删去
     */
    bool moveu_C_to_X(ui S_end, ui &C_end, ui u, ui level)
    {
        assert(Qv.empty());
        if (u != SC[S_end]) return false;
        assert(level_id[u] == n);
        Qv.push(u);
        level_id[u] = level;
        return remove_vertices(S_end, C_end, level);
    }
    /**
     * 更新pivot集合
     *
     * @param S_end S集合的末尾索引
     * @param C_end C集合的末尾索引
     * @param u 当前处理的顶点
     * @param new_pivot_set 用于存储新的pivot点的向量
     *
     * 该函数用于更新pivot集合,主要包含两类顶点:
     * 1. 与顶点u不相邻的顶点
     * 2. 与顶点u相邻且在S中有共同非邻居的顶点
     */
    void update_pivot_set(ui S_end, ui C_end, ui u, vector<ui> &new_pivot_set)
    {
        new_pivot_set.clear();
        int *t_matrix = matrix + u * n;
        for (ui i = S_end; i < C_end; i++) {
            ui v = SC[i];
            if (!t_matrix[v]) {
                new_pivot_set.push_back(v);
            }
            else {
                //如果v和u在S中有共同的非邻居，那么v就是一个pivot点
                bool has_common_nonneighbor = false;
                for (ui j = 0; j < S_end; j++) {
                    ui w = SC[j];
                    if (!matrix[v * n + w] && !matrix[u * n + w]) {
                        has_common_nonneighbor = true;
                        break;
                    }
                }
                if (has_common_nonneighbor) {
                    new_pivot_set.push_back(v);
                }
            }
        }
    }

    void swap_pos(ui i, ui j)
    {
        swap(SC[i], SC[j]);
        SC_rid[SC[i]] = i;
        SC_rid[SC[j]] = j;
    }

    ui upper_bound(ui S_end, ui C_end)
    {
        return C_end;
        // ui ub1 = upper_bound_based_partition_1(S_end, C_end);
        // ui ub2 = upper_bound_based_partition_2(S_end, C_end);
        // ui ub3 = upper_bound_based_partition_3(S_end, C_end);
        // // assert(ub1 <= ub2);
        // // assert(ub2 == ub3);

        // return min(ub1, min(ub2, ub3));
    }

    ui upper_bound_based_partition_1(ui S_end, ui C_end)
    {
        ui ub = 0, pi0 = C_end - S_end;
        ui *count = neighbors;
        ui *C_missing_edges = nonneighbors;
        memset(count, 0, sizeof(ui) * n);
        memset(vis, 0, sizeof(ui) * n);
        for (ui i = 0; i < S_end; i++) {
            assert(K >= S_end - degree_in_S[SC[i]]);
            C_missing_edges[i] = K - (S_end - degree_in_S[SC[i]]);
            for (ui j = S_end; j < C_end; j++)
                if (!matrix[SC[i] * n + SC[j]])
                    count[i]++;
            // printf("count[%d] = %d\n", SC[i], count[i]);
        }
        for (ui i = 0; i < S_end; i++) {
            if (count[i] <= C_missing_edges[i]) {
                vis[i] = 1;
            }
        }
        while (1) {
            // find max dise
            int uid = -1;
            double max_dise = 0;
            for (ui i = 0; i < S_end; i++)
                if (!vis[i]) {
                    double dise = (C_missing_edges[i] == 0) ? 0 : (1.0 * count[i] / C_missing_edges[i]);
                    if (max_dise < dise) {
                        max_dise = dise;
                        uid = i;
                    }
                }
            if (uid == -1 || max_dise <= 1)
                break;
            vis[uid] = 1;
            ui u = SC[uid];
            // printf("u = %d, C_missing_edges = %d, count = %d\n", u, C_missing_edges[uid], count[uid]);
            // printf("u = %d, contribution = %d, max_dise = %.6lf\n", SC[uid], min(C_missing_edges[uid], count[uid]), max_dise);
            ub = ub + min(C_missing_edges[uid], count[uid]);

            for (ui j = S_end; j < C_end; j++) {
                if (!vis[j] && !matrix[u * n + SC[j]]) {
                    vis[j] = 1;
                    for (ui k = 0; k < S_end; k++) {
                        if (!vis[k] && !matrix[SC[k] * n + SC[j]]) {
                            assert(count[k] > 0);
                            count[k]--;
                        }
                    }
                    pi0--;
                }
            }
        }

        int total = 0;
        for (ui i = S_end; i < C_end; i++) {
            if (!vis[i])
                total++;
        }
        // printf("C - total = %d\n", C_end - total);

        ui *ub2_vertices = neighbors;
        ui ub2_vertices_n = 0;
        ui *colors = nonneighbors;
        memset(colors, 0, sizeof(ui) * n);
        for (ui i = S_end; i < C_end; i++) {
            if (!vis[i])
                ub2_vertices[++ub2_vertices_n] = SC[i];
        }

        // printf("S_end = %d, pi0 = %d\n", S_end, pi0);
        ub = S_end + pi0 + ub;
        // printf("S_end = %u, C_end = %u, ub = %u, kplex = %u, total = %d\n", S_end, C_end, ub, best_solution_size, total);
        assert(ub <= C_end);
        return ub;
    }

    ui upper_bound_based_partition_2(ui S_end, ui C_end)
    {
        ui ub = 0, pi0 = C_end - S_end;
        ui *count = neighbors;
        ui *C_missing_edges = nonneighbors;
        memset(count, 0, sizeof(ui) * n);
        memset(vis, 0, sizeof(ui) * n);
        for (ui i = 0; i < S_end; i++) {
            C_missing_edges[i] = K - (S_end - degree_in_S[SC[i]]);
            for (ui j = S_end; j < C_end; j++)
                if (!matrix[SC[i] * n + SC[j]])
                    count[i]++;
            printf("count[%d] = %d\n", SC[i], count[i]);
        }
        ui co = 0;
        for (ui i = 0; i < S_end; i++) {
            // find max dise
            int uid = i;
            double max_dise = 0;
            vis[uid] = 1;
            ui u = SC[uid];
            printf("u = %d, C_missing_edges = %d, count = %d\n", u, C_missing_edges[uid], count[uid]);
            // printf("u = %d, contribution = %d, max_dise = %.6lf\n", SC[uid], min(C_missing_edges[uid], count[uid]), max_dise);
            ub = ub + min(C_missing_edges[uid], count[uid]);

            for (ui j = S_end; j < C_end; j++) {
                if (!vis[j] && !matrix[u * n + SC[j]]) {
                    vis[j] = 1;
                    for (ui k = 0; k < S_end; k++) {
                        if (!vis[k] && !matrix[SC[k] * n + SC[j]]) {
                            assert(count[k] > 0);
                            count[k]--;
                        }
                    }
                    pi0--;
                }
            }
        }
        printf("S_end = %d, pi0 = %d\n", S_end, pi0);
        ub = S_end + pi0 + ub;

        assert(ub <= C_end);
        return ub;
    }

    ui upper_bound_based_partition_3(ui S_end, ui C_end)
    {
        ui ub = 0, pi0 = C_end - S_end;
        ui *C_missing_edges = nonneighbors;
        memset(vis, 0, sizeof(ui) * n);

        for (ui i = 0; i < S_end; i++) {
            C_missing_edges[i] = S_end - degree_in_S[SC[i]] - 1;
        }
        for (ui i = 0; i < S_end; i++) {
            ui u = SC[i];
            ui pii = 0;
            int *t_matrix = matrix + u * n;
            for (ui j = S_end; j < C_end; j++) {
                if (!vis[j] && !t_matrix[SC[j]]) {
                    vis[j] = 1;
                    pii++;
                    pi0--;
                }
            }
            ub = ub + min(K - 1 - C_missing_edges[i], pii);
        }
        ub = S_end + pi0 + ub;
        // printf("S_end = %u, C_end = %u, ub = %u, kplex = %u, total = %d\n", S_end, C_end, ub, best_solution_size, total);
        assert(ub <= n);
        return ub;
    }
#ifdef _SECOND_ORDER_PRUNING_
    bool upper_bound_based_prune(ui S_end, ui u, ui v)
    {
        // ui ub = S_end + 2*K - (S_end - degree_in_S[u]) - (S_end - degree_in_S[v]) + cn[u*n + v];
        ui ub = 2 * K + degree_in_S[u] + degree_in_S[v] + cn[u * n + v] - S_end;
        if (SC_rid[u] >= S_end) {
            --ub; // S_end ++
            if (matrix[u * n + v]) ++ub; // degree_in_S[v] ++
        }
        if (SC_rid[v] >= S_end) {
            --ub;
            if (matrix[v * n + u]) ++ub;
        }
        return ub <= best_solution_size;
    }
#endif
    bool check_balance(ui S_end, ui u)
    {
        vector<ui> nei;
        int *t_matrix = matrix + u * n;
        for (ui i = 0; i < S_end; i++) if (t_matrix[SC[i]]) nei.push_back(SC[i]);
        ui neighbors_n = nei.size();
        for (ui i = 0; i < neighbors_n; i++) for (ui j = i + 1; j < neighbors_n; j++) {
            ui v = nei[i], w = nei[j];
            if (!matrix[v * n + w]) continue;
            ui tri_sum = matrix[v * n + w] + t_matrix[v] + t_matrix[w];
            assert(tri_sum == 3 || tri_sum == 1 || tri_sum == -1 || tri_sum == -3);
            if (tri_sum == 1 || tri_sum == -3) return false;
        }
        return true;
    }
};

#endif /* _SIGNED_KPLEX_ */