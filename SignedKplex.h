/*
This file contains code from the Maximum-kPlex project, which is licensed under the MIT License.
The original code and license can be found at: https://github.com/LijunChang/Maximum-kPlex
*/

#ifndef _SIGNED_KPLEX_
#define _SIGNED_KPLEX_

#include "Utility.h"
#include "Timer.h"

class SIGNED_KPLEX
{
private:
    ui K;
    ui n;
    ept m;
    int *matrix;
    long long matrix_size;
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
    bool *s_matrix_flag;

public:
    SIGNED_KPLEX()
    {
        n = m = 0;
        matrix = NULL;
        matrix_size = 0;
        best_solution_size = 0;
        best_solution = NULL;

        degree = NULL;
        degree_in_S = NULL;
        neighbors = NULL;
        nonneighbors = NULL;
        SC = NULL;
        SC_rid = NULL;
        level_id = NULL;
        vis = NULL;
        s_matrix_flag = NULL;
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
        if (s_matrix_flag != NULL) {
            delete[] s_matrix_flag;
            s_matrix_flag = NULL;
        }
    }

    void allocateMemory(ui n, ui m)
    {
        matrix_size = m * 2;
        matrix = new int[matrix_size];
        best_solution = new ui[n];
        degree = new ui[n];
        degree_in_S = new ui[n];
        neighbors = new ui[n];
        nonneighbors = new ui[n];
        SC = new ui[n];
        SC_rid = new ui[n];
        level_id = new ui[n];
        vis = new ui[n];
        s_matrix_flag = new bool[matrix_size];
        level_cnt.resize(n, 0);
        level_cnt_after_ub.resize(n, 0);
    }

    void load_graph(ui _n, const vector<Edge> &vp)
    {
        n = _n;
        if (1ll * n * n > matrix_size) {
            while (1ll * n * n > matrix_size) matrix_size *= 2;
            delete[] matrix;
            matrix = new int[matrix_size];
            delete[] s_matrix_flag;
            s_matrix_flag = new bool[matrix_size];
        }

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
                if ((degree[v]--) + K == best_solution_size) q.push(v);
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

        assert(Qv.empty());

        // 将choose_u加入S
        if (choose_u != -1) {
            bool pruned = moveu_C_to_S(S_end, C_end, 0, choose_u);
            assert(pruned == false);
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
        for (ui i = 0; i < C_end; i++) assert(degree[SC[i]] + K > best_solution_size);
        for (ui i = 0; i < S_end; i++) assert(degree_in_S[SC[i]] + K >= S_end);
        for (ui i = S_end; i < C_end; i++) assert(degree_in_S[SC[i]] + K > S_end);
        for (ui i = 0; i < C_end; i++) assert(level_id[SC[i]] == n);
        for (ui i = 0; i < C_end; i++) assert(check_balance(S_end, SC[i]));
#endif

        dfs_cnt++;
        level_cnt[level]++;

        // upper bound
        // ui ub = upper_bound(S_end, C_end);
        // if (ub <= best_solution_size) return;

        dfs_cnt_after_ub++;
        level_cnt_after_ub[level]++;

        ui old_kplex_size = best_solution_size, old_C_end = C_end;

        // choose branching vertex
        ui u = choose_branch_vertex_with_min_degree(S_end, C_end, pivot_set);
        assert(u != n);
        assert(SC[SC_rid[u]] == u && SC_rid[u] >= S_end && SC_rid[u] < C_end);
        assert(degree[u] + K > best_solution_size);
        assert(check_balance(S_end, u));

        // the first branch includes u into S
        vector<ui> new_pivot_set;
        bool pruned = moveu_C_to_S(S_end, C_end, level, u);
        // 加入点u时，向kplex_search函数传入空的pivot_set
        if (!pruned) {
            new_pivot_set.clear();
            kplex_search(S_end, C_end, level + 1, new_pivot_set);
        }

        restore_C(S_end, C_end, old_C_end, level);
        moveu_S_to_C(S_end, C_end, level);
        assert(C_end == old_C_end);
        assert(Qv.empty());
        assert(u == SC[S_end]);

        // the second branch exclude u from S
        pruned = false;
        if (best_solution_size > old_kplex_size) {
            if (C_end <= best_solution_size) return;
            // ub = upper_bound(S_end, C_end);
            // if(ub <= best_solution_size) return;
            pruned = reduce_SC_based_lb(S_end, C_end, level);
        }
        if (!pruned) pruned = moveu_C_to_X(S_end, C_end, u, level);
        if (!pruned) {
            new_pivot_set.clear();
            update_pivot_set(S_end, C_end, u, new_pivot_set);
            kplex_search(S_end, C_end, level + 1, new_pivot_set);
        }
        restore_C(S_end, C_end, old_C_end, level);
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
     */
    bool moveu_C_to_S(ui &S_end, ui &C_end, ui level, ui u)
    {
        // 1
        swap_pos(SC_rid[u], S_end++);
        assert(u == SC[S_end - 1]);

        // 2
        ui neighbors_n = 0, nonneighbors_n = 0;
        for (ui i = 0; i < C_end; i++) if (SC[i] != u) {
            if (matrix[u * n + SC[i]]) neighbors[neighbors_n++] = SC[i];
            else nonneighbors[nonneighbors_n++] = SC[i];
        }
        assert(neighbors_n + nonneighbors_n == C_end - 1);
        for (ui i = 0; i < neighbors_n; i++) degree_in_S[neighbors[i]]++;

        // 3
        assert(Qv.empty());
        if (degree_in_S[u] + K == S_end) {
            ui i = 0;
            while (i < nonneighbors_n && SC_rid[nonneighbors[i]] < S_end) i++;
            for (; i < nonneighbors_n; i++) {
                level_id[nonneighbors[i]] = level;
                Qv.push(nonneighbors[i]);
            }
        }
        else {
            ui i = 0;
            while (i < nonneighbors_n && SC_rid[nonneighbors[i]] < S_end) i++;
            for (; i < nonneighbors_n; i++) if (degree_in_S[nonneighbors[i]] + K <= S_end) {
                level_id[nonneighbors[i]] = level;
                Qv.push(nonneighbors[i]);
            }
        }
        for (ui i = 0; i < nonneighbors_n && SC_rid[nonneighbors[i]] < S_end; i++) {
            if (degree_in_S[nonneighbors[i]] + K == S_end) {
                int *t_matrix = matrix + nonneighbors[i] * n;
                for (ui j = S_end; j < C_end; j++)
                    if (level_id[SC[j]] == n && !t_matrix[SC[j]]) {
                        level_id[SC[j]] = level;
                        Qv.push(SC[j]);
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
                    assert(neighbors[i] != u);
                    Qv.push(neighbors[i]);
                }
            }
        }
        return remove_vertices(S_end, C_end, level);
    }
    /**
     * 从候选集C中移除顶点
     */
    bool remove_vertices(ui S_end, ui &C_end, ui level)
    {
        while (!Qv.empty()) {
            // remove u from C
            ui u = Qv.front(); Qv.pop();
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
                neighbors[neighbors_n++] = SC[i];
                degree[v]--;
                if (degree[v] + K <= best_solution_size) {
                    if (i < S_end) terminate = true;
                    else if (level_id[v] == n) {
                        level_id[v] = level;
                        Qv.push(v);
                    }
                }
            }
            if (terminate) return true;
        }
        return false;
    }
    /**
     * 恢复C集合中的顶点
     */
    void restore_C(ui S_end, ui &C_end, ui old_C_end, ui level)
    {
        while (!Qv.empty()) {
            ui u = Qv.front(); Qv.pop();
            assert(level_id[u] == level && SC_rid[u] < C_end);
            level_id[u] = n;
        }
        while (C_end < old_C_end) {
            ui u = SC[C_end];
            assert(level_id[u] == level && SC_rid[u] == C_end);
            level_id[u] = n;
            for (ui i = 0; i < C_end; i++) if (matrix[u * n + SC[i]]) degree[SC[i]]++;
            C_end++;
        }
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
    }
    /**
     * 基于下界剪枝来化简SC集合
     */
    bool reduce_SC_based_lb(ui S_end, ui &C_end, ui level)
    {
        assert(Qv.empty());
        for (ui i = 0; i < S_end; i++) if (degree[SC[i]] + K <= best_solution_size) return true;
        for (ui i = S_end; i < C_end; i++) if (degree[SC[i]] + K <= best_solution_size) {
            level_id[SC[i]] = level;
            Qv.push(SC[i]);
        }
        return remove_vertices(S_end, C_end, level);
    }
    /**
     * 将顶点u从集合SC中删去
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