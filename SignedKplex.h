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

    int total = 0;
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
        ui C_end = 0;
        init(C_end, choose_u);
        if (C_end)
            kplex_search(0, C_end, 1, choose_u);
        if (best_solution_size > kplex.size()) {
            kplex.clear();
            for (int i = 0; i < best_solution_size; i++)
                kplex.push_back(best_solution[i]);
        }
    }

private:
    // init S, C
    void init(ui &C_end, int choose_u)
    {
        C_end = 0;
        // k-core
        queue<ui> q;
        memset(vis, 0, sizeof(ui) * n);
        for (ui i = 0; i < n; i++)
            if (degree[i] + K <= best_solution_size)
                q.push(i);
        while (!q.empty()) {
            ui u = q.front();
            q.pop();
            if (vis[u])
                continue;
            vis[u] = 1;
            for (ui v = 0; v < n; v++)
                if (matrix[u * n + v]) {
                    assert(degree[v]);
                    degree[v]--;
                    if (degree[v] + K == best_solution_size)
                        q.push(v);
                }
        }

        if (choose_u != -1 && vis[choose_u])
            return;

        for (ui i = 0; i < n; i++)
            SC_rid[i] = n;

        for (ui i = 0; i < n; i++)
            if (!vis[i]) {
                SC[C_end] = i;
                SC_rid[i] = C_end++;
            }

        memset(degree, 0, sizeof(ui) * n);
        for (ui i = 0; i < C_end; i++) {
            ui u = SC[i];
            assert(SC_rid[u] == i);
            for (ui j = 0; j < C_end; j++)
                if (matrix[u * n + SC[j]])
                    degree[u]++;
        }

        memset(level_id, 0, sizeof(ui) * n);
        for (ui i = 0; i < C_end; i++)
            level_id[SC[i]] = n;

        assert(Qv.empty());
    }

    void kplex_search(ui S_end, ui C_end, ui level, int choose_u)
    {
        // printf("S_end = %d, C_end = %d, level = %d\n", S_end, C_end, level);
        if (S_end > best_solution_size) {
            best_solution_size = S_end;
            for (ui i = 0; i < best_solution_size; i++)
                best_solution[i] = SC[i];
#ifndef NDEBUG
            for (ui i = 0; i < best_solution_size; i++)
                assert(degree_in_S[best_solution[i]] + K >= best_solution_size);
            printf("Find a k-plex of size: %u\n", best_solution_size);
#endif
        }
        if (C_end <= best_solution_size)
            return;

        // upper bound
        // ui ub1 = upper_bound_based_partition_1(S_end, C_end);
        // ui ub2 = upper_bound_based_partition_2(S_end, C_end);
        // ui ub3 = upper_bound_based_partition_3(S_end, C_end);
        // assert(!(ub1 <= best_solution_size && ub2 > best_solution_size));
        // assert(ub1 <= ub2);
        // assert(ub2 == ub3);
        // if(ub1 <= best_solution_size) return;
        // if(ub2 <= best_solution_size) return;
        // if (ub3 <= best_solution_size) return;

        ui pre_kplex_size = best_solution_size, old_C_end = C_end;
        total++;
        // printf("S_end = %d, C_end = %d, kplex = %d, total = %d\n", S_end, C_end, best_solution_size, total);

#ifndef NDEBUG
        for (ui i = 0; i < C_end; i++)
            assert(degree[SC[i]] + K > best_solution_size);
        for (ui i = 0; i < S_end; i++)
            assert(degree_in_S[SC[i]] + K >= S_end);
        for (ui i = S_end; i < C_end; i++)
            assert(degree_in_S[SC[i]] + K > S_end);
        for (ui i = 0; i < C_end; i++)
            assert(level_id[SC[i]] == n);
        for (ui i = 0; i < C_end; i++)
            assert(check_balance(S_end, SC[i]));
#endif
        // choose branching vertex
        bool must_choose = false;
        ui u = n;
        if (choose_u != -1) {
            assert(S_end == 0);
            assert(SC_rid[choose_u] < C_end);
            u = choose_u;
            must_choose = true;
        }
        else {
            u = choose_branch_vertex(S_end, C_end);
        }
        assert(u != n);
        assert(SC[SC_rid[u]] == u && SC_rid[u] >= S_end && SC_rid[u] < C_end);
        assert(degree[u] + K > best_solution_size);
        assert(check_balance(S_end, u));

        // the first branch includes u into S
        bool pruned = move_u_to_S(S_end, C_end, level, u);
        if (!pruned)
            kplex_search(S_end, C_end, level + 1, -1);
        restore_C(S_end, C_end, old_C_end, level);
        move_u_to_C(S_end, C_end, level);
        assert(C_end == old_C_end);
        assert(Qv.empty());
        assert(u == SC[S_end]);

        // the second branch exclude u from S
        if (must_choose)
            return;
        pruned = false;
        if (best_solution_size > pre_kplex_size) {
            if (C_end <= best_solution_size)
                return;
            // ub = upper_bound(S_end, C_end);
            // if(ub <= best_solution_size) return;
            pruned = collect_removable_vertices(S_end, C_end, level);
        }

        if (!pruned)
            pruned = remove_u_from_C(S_end, C_end, u, level);
        if (!pruned) {
#ifndef NDEBUG
            for (ui i = 0; i < C_end; i++)
                assert(degree[SC[i]] + K > best_solution_size);
            for (ui i = 0; i < S_end; i++)
                assert(degree_in_S[SC[i]] + K >= S_end);
            for (ui i = S_end; i < C_end; i++)
                assert(degree_in_S[SC[i]] + K > S_end);
            for (ui i = 0; i < C_end; i++)
                assert(level_id[SC[i]] == n);
#endif
            kplex_search(S_end, C_end, level + 1, -1);
        }
        restore_C(S_end, C_end, old_C_end, level);
    }

    ui choose_branch_vertex(ui S_end, ui C_end)
    {
        return SC[S_end];
        // ui u = n, min_degree_in_S = n;
        // for(ui i = 0; i < C_end; i++) {
        //     ui v = SC[i];
        //     // if(degree[v] + K >= C_end) continue;
        //     if(degree_in_S[v] < min_degree_in_S) {
        //         u = v;
        //         if(degree[v] + K >= C_end) {
        //             min_degree_in_S = degree_in_S[v];
        //         }
        //     }
        // }
        // assert(u != n);
        // if(min_degree_in_S == n) {
        //     total += 100000;
        // }
        // if(SC_rid[u] < S_end) {
        //     int *t_matrix = matrix+u*n;
        //     u = n;
        //     ui max_degree = 0;
        //     for(ui i = S_end; i < C_end; i++) if(!t_matrix[SC[i]]) {
        //         if(degree[SC[i]] > max_degree) {
        //             max_degree = degree[SC[i]];
        //             u = SC[i];
        //         }
        //     }
        // }
        // if(u == n) u = SC[S_end];
        // assert(u != n);
        // return u;
    }

    bool move_u_to_S(ui &S_end, ui &C_end, ui level, ui u)
    {
        swap_pos(SC_rid[u], S_end++);
        assert(u == SC[S_end - 1]);

        ui neighbors_n = 0, nonneighbors_n = 0;
        for (ui i = 0; i < C_end; i++)
            if (SC[i] != u) {
                if (matrix[u * n + SC[i]])
                    neighbors[neighbors_n++] = SC[i];
                else
                    nonneighbors[nonneighbors_n++] = SC[i];
            }
        assert(neighbors_n + nonneighbors_n == C_end - 1);
        for (ui i = 0; i < neighbors_n; i++)
            degree_in_S[neighbors[i]]++;

        assert(Qv.empty());
        // reduce
        if (degree_in_S[u] + K == S_end) {
            ui i = 0;
            while (i < nonneighbors_n && SC_rid[nonneighbors[i]] < S_end)
                i++;
            for (; i < nonneighbors_n; i++) {
                level_id[nonneighbors[i]] = level;
                Qv.push(nonneighbors[i]);
            }
        }
        else {
            ui i = 0;
            while (i < nonneighbors_n && SC_rid[nonneighbors[i]] < S_end)
                i++;
            for (; i < nonneighbors_n; i++)
                if (degree_in_S[nonneighbors[i]] + K <= S_end) {
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

        // check balance
        {
            ui i = 0;
            while (i < neighbors_n && SC_rid[neighbors[i]] < S_end)
                i++;
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

    // old_C_end->C_end
    bool remove_vertices(ui S_end, ui &C_end, ui level)
    {
        while (!Qv.empty()) {
            // remove u from C
            ui u = Qv.front();
            Qv.pop();
            assert(level_id[u] == level);
            assert(SC[SC_rid[u]] == u);
            assert(SC_rid[u] >= S_end && SC_rid[u] < C_end);
            swap_pos(SC_rid[u], --C_end);
            assert(u == SC[C_end]);

            int *t_matrix = matrix + u * n;
            bool terminate = false;

            ui neighbors_n = 0;
            for (ui i = 0; i < C_end; i++)
                if (t_matrix[SC[i]]) {
                    ui v = SC[i];
                    neighbors[neighbors_n++] = SC[i];
                    degree[v]--;
                    if (degree[v] + K <= best_solution_size) {
                        if (i < S_end)
                            terminate = true;
                        else if (level_id[v] == n) {
                            level_id[v] = level;
                            Qv.push(v);
                        }
                    }
                }
            if (terminate)
                return true;
        }
        return false;
    }

    // C_end->old_C_end
    void restore_C(ui S_end, ui &C_end, ui old_C_end, ui level)
    {
        while (!Qv.empty()) {
            ui u = Qv.front();
            Qv.pop();
            assert(level_id[u] == level && SC_rid[u] < C_end);
            level_id[u] = n;
        }
        while (C_end < old_C_end) {
            ui u = SC[C_end];
            assert(level_id[u] == level && SC_rid[u] == C_end);
            level_id[u] = n;
            for (ui i = 0; i < C_end; i++)
                if (matrix[u * n + SC[i]])
                    degree[SC[i]]++;
            C_end++;
        }
    }

    // u: S->C
    void move_u_to_C(ui &S_end, ui C_end, ui level)
    {
        assert(S_end);
        ui u = SC[--S_end];
        ui neighbors_n = 0;
        int *t_matrix = matrix + u * n;
        for (ui i = 0; i < C_end; i++)
            if (t_matrix[SC[i]])
                degree_in_S[SC[i]]--;
    }

    bool collect_removable_vertices(ui S_end, ui &C_end, ui level)
    {
        assert(Qv.empty());
        for (ui i = 0; i < S_end; i++)
            if (degree[SC[i]] + K <= best_solution_size)
                return true;
        for (ui i = S_end; i < C_end; i++)
            if (degree[SC[i]] + K <= best_solution_size) {
                level_id[SC[i]] = level;
                Qv.push(SC[i]);
            }
        return remove_vertices(S_end, C_end, level);
    }

    bool remove_u_from_C(ui S_end, ui &C_end, ui u, ui level)
    {
        // assert(Qv.empty());
        // if (u != SC[S_end]) return false;
        // swap_pos(S_end, --C_end);
        // assert(u == SC[C_end]);
        // assert(level_id[u] == n);
        // level_id[u] = level;

        // ui neighbors_n = 0, nonneighbors_n = 0;
        // for (ui i = 0; i < C_end; i++) {
        //     assert(SC[i] != u);
        //     if (matrix[u * n + SC[i]])
        //         neighbors[neighbors_n++] = SC[i];
        //     else
        //         nonneighbors[nonneighbors_n++] = SC[i];
        // }
        // assert(neighbors_n + nonneighbors_n == C_end);
        // for (ui i = 0; i < neighbors_n; i++)
        //     degree[neighbors[i]]--;

        // for (ui i = 0; i < neighbors_n; i++) {
        //     ui v = neighbors[i];
        //     if (degree[v] + K <= best_solution_size) {
        //         if (i < S_end)
        //             return true;
        //         else {
        //             assert(level_id[v] == n);
        //             level_id[v] = level;
        //             Qv.push(v);
        //         }
        //     }
        // }
        // return remove_vertices(S_end, C_end, level);
        assert(Qv.empty());
        if (u != SC[S_end])
            return false;
        assert(level_id[u] == n);
        Qv.push(u);
        level_id[u] = level;
        return remove_vertices(S_end, C_end, level);
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
        printf("S_end = %u, C_end = %u, ub = %u, kplex = %u, total = %d\n", S_end, C_end, ub, best_solution_size, total);
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
        ui neighbors_n = 0;
        ui *nei = vis;
        int *t_matrix = matrix + u * n;
        for (ui i = 0; i < S_end; i++)
            if (t_matrix[SC[i]])
                nei[neighbors_n++] = SC[i];
        for (ui i = 0; i < neighbors_n; i++)
            for (ui j = i + 1; j < neighbors_n; j++) {
                ui v = nei[i], w = nei[j];
                if (!matrix[v * n + w])
                    continue;
                ui tri = matrix[v * n + w] + t_matrix[v] + t_matrix[w];
                assert(tri == 3 || tri == 1 || tri == -1 || tri == -3);
                if (tri == 1 || tri == -3)
                    return false;
            }
        return true;
    }
};

#endif /* _SIGNED_KPLEX_ */

// /*
// This file contains code from the Maximum-kPlex project, which is licensed under the MIT License.
// The original code and license can be found at: https://github.com/LijunChang/Maximum-kPlex
// */

// #ifndef _SIGNED_KPLEX_
// #define _SIGNED_KPLEX_

// #include "Utility.h"
// #include "Timer.h"

// class SIGNED_KPLEX
// {
// private:
//     ui K;
//     ui n;
//     ept m;
//     int* matrix;
//     long long matrix_size;
//     ui best_solution_size;
//     ui* best_solution;
//     ui* degree;
//     ui* degree_in_S;
//     ui* neighbors;
//     ui* nonneighbors;
//     ui* SC;
//     ui* SC_rid;
//     ui* level_id;
//     ui* vis;
//     queue<ui> Qv;

//     int total = 0;
//     bool* s_matrix_flag;

// public:
//     SIGNED_KPLEX()
//     {
//         n = m = 0;
//         matrix = NULL;
//         matrix_size = 0;
//         best_solution_size = 0;
//         best_solution = NULL;

//         degree = NULL;
//         degree_in_S = NULL;
//         neighbors = NULL;
//         nonneighbors = NULL;
//         SC = NULL;
//         SC_rid = NULL;
//         level_id = NULL;
//         vis = NULL;
//         s_matrix_flag = NULL;
//     }

//     ~SIGNED_KPLEX()
//     {
//         if (matrix != NULL) { delete[] matrix; matrix = NULL; }
//         if (degree != NULL) { delete[] degree; degree = NULL; }
//         if (best_solution != NULL) { delete[] best_solution; best_solution = NULL; }
//         if (degree_in_S != NULL) { delete[] degree_in_S; degree_in_S = NULL; }
//         if (neighbors != NULL) { delete[] neighbors; neighbors = NULL; }
//         if (nonneighbors != NULL) { delete[] nonneighbors; nonneighbors = NULL; }
//         if (SC != NULL) { delete[] SC; SC = NULL; }
//         if (SC_rid != NULL) { delete[] SC_rid; SC_rid = NULL; }
//         if (level_id != NULL) { delete[] level_id; level_id = NULL; }
//         if (vis != NULL) { delete[] vis; vis = NULL; }
//         if (s_matrix_flag != NULL) { delete[] s_matrix_flag; s_matrix_flag = NULL; }
//     }

//     void allocateMemory(ui n, ui m)
//     {
//         matrix_size = m * 2;
//         matrix = new int[matrix_size];
//         best_solution = new ui[n];
//         degree = new ui[n];
//         degree_in_S = new ui[n];
//         neighbors = new ui[n];
//         nonneighbors = new ui[n];
//         SC = new ui[n];
//         SC_rid = new ui[n];
//         level_id = new ui[n];
//         vis = new ui[n];
//         s_matrix_flag = new bool[matrix_size];
//     }

//     void load_graph(ui _n, const vector<Edge>& vp)
//     {
//         n = _n;
//         if (1ll * n * n > matrix_size) {
//             while (1ll * n * n > matrix_size) matrix_size *= 2;
//             delete[] matrix;
//             matrix = new int[matrix_size];
//             delete[] s_matrix_flag;
//             s_matrix_flag = new bool[matrix_size];
//         }

//         memset(matrix, 0, sizeof(int) * matrix_size);
//         memset(degree, 0, sizeof(ui) * n);
//         memset(degree_in_S, 0, sizeof(ui) * n);

//         m = vp.size();
//         for (ui i = 0; i < m; i++) {
//             int a = vp[i].a, b = vp[i].b, c = vp[i].c;
//             assert(a < n);
//             assert(b < n);
//             assert(c == 1 || c == -1);
//             degree[a]++;
//             degree[b]++;
//             matrix[a * n + b] = matrix[b * n + a] = c;
//         }

// #ifndef NDEBUG
//         printf("matrix kplex load_graph: n=%u, m=%lu\n", n, m);
// #endif
//     }

//     void kPlex(ui K_, vector<ui>& kplex, int choose_u)
//     {
//         K = K_;
//         best_solution_size = kplex.size();
//         ui C_end = 0;
//         // remove unpromising nodes
//         init(C_end, choose_u);
//         // search
//         printf("\tn = %d, C_end = %d\n", n, C_end);
//         if (C_end) kplex_search(0, C_end, 1, choose_u);
//         if (best_solution_size > kplex.size()) {
//             kplex.clear();
//             for (int i = 0; i < best_solution_size; i++) kplex.push_back(best_solution[i]);
//         }
//     }

// private:
//     // init S, C
//     void init(ui& C_end, int choose_u)
//     {
//         C_end = 0;
//         // k-core
//         queue<ui> q;
//         memset(vis, 0, sizeof(ui) * n);
//         for (ui i = 0; i < n; i++) if (degree[i] + K <= best_solution_size) q.push(i);
//         while (!q.empty()) {
//             ui u = q.front(); q.pop();
//             if (vis[u]) continue;
//             vis[u] = 1;
//             for (ui v = 0; v < n; v++) if (matrix[u * n + v]) {
//                 assert(degree[v]);
//                 degree[v]--;
//                 if (degree[v] + K == best_solution_size) q.push(v);
//             }
//         }

//         if (choose_u != -1 && vis[choose_u]) return;

//         for (ui i = 0; i < n; i++) SC_rid[i] = n;

//         for (ui i = 0; i < n; i++) if (!vis[i]) {
//             SC[C_end] = i;
//             SC_rid[i] = C_end++;
//         }

//         memset(degree, 0, sizeof(ui) * n);
//         for (ui i = 0; i < C_end; i++) {
//             ui u = SC[i];
//             assert(SC_rid[u] == i);
//             for (ui j = 0; j < C_end; j++) if (matrix[u * n + SC[j]]) degree[u]++;
//         }

//         memset(level_id, 0, sizeof(ui) * n);
//         for (ui i = 0; i < C_end; i++) level_id[SC[i]] = n;

//         assert(Qv.empty());
//     }

//     void kplex_search(ui S_end, ui C_end, ui level, int choose_u)
//     {
//         // printf("S_end = %d, C_end = %d, level = %d\n", S_end, C_end, level);
//         if (S_end > best_solution_size) {
//             best_solution_size = S_end;
//             for (ui i = 0; i < best_solution_size; i++) best_solution[i] = SC[i];
// #ifndef NDEBUG
//             for (ui i = 0; i < best_solution_size; i++) assert(degree_in_S[best_solution[i]] + K >= best_solution_size);
//             printf("Find a k-plex of size: %u\n", best_solution_size);
// #endif
//         }
//         if (C_end <= best_solution_size) return;

//         // upper bound
//         ui ub = upper_bound(S_end, C_end);
//         if(ub <= best_solution_size) return;

//         ui pre_kplex_size = best_solution_size, old_C_end = C_end;
//         total++;
//         // printf("S_end = %d, C_end = %d, kplex = %d, total = %d\n", S_end, C_end, best_solution_size, total);

// #ifndef NDEBUG
//         for (ui i = 0; i < C_end; i++) assert(degree[SC[i]] + K > best_solution_size);
//         for (ui i = 0; i < S_end; i++) assert(degree_in_S[SC[i]] + K >= S_end);
//         for (ui i = S_end; i < C_end; i++) assert(degree_in_S[SC[i]] + K > S_end);
//         for (ui i = 0; i < C_end; i++) assert(level_id[SC[i]] == n);
//         for (ui i = 0; i < C_end; i++) assert(check_balance(S_end, SC[i]));
// #endif
//         // choose branching vertex
//         bool must_choose = false;
//         ui u = n;
//         if (choose_u != -1) {
//             assert(S_end == 0);
//             assert(SC_rid[choose_u] < C_end);
//             u = choose_u;
//             must_choose = true;
//         }
//         else {
//             u = choose_branch_vertex(S_end, C_end);
//         }
//         assert(u != n);
//         assert(SC[SC_rid[u]] == u && SC_rid[u] >= S_end && SC_rid[u] < C_end);
//         assert(degree[u] + K > best_solution_size);
//         assert(check_balance(S_end, u));

//         // the first branch includes u into S
//         bool pruned = move_u_to_S(S_end, C_end, level, u);
//         if (!pruned) kplex_search(S_end, C_end, level + 1, -1);
//         restore_C(S_end, C_end, old_C_end, level);
//         move_u_to_C(S_end, C_end, level);
//         assert(C_end == old_C_end);
//         assert(Qv.empty());
//         assert(u == SC[S_end]);

//         // the second branch exclude u from S
//         if (must_choose) return;
//         pruned = false;
//         if (best_solution_size > pre_kplex_size) {
//             if (C_end <= best_solution_size) return;
//             ub = upper_bound(S_end, C_end);
//             if(ub <= best_solution_size) return;
//             pruned = collect_removable_vertices(S_end, C_end, level);
//         }

//         if (!pruned) pruned = remove_u_from_C(S_end, C_end, u, level);
//         if (!pruned) {
// #ifndef NDEBUG
//             for (ui i = 0; i < C_end; i++)
//                 assert(degree[SC[i]] + K > best_solution_size);
//             for (ui i = 0; i < S_end; i++)
//                 assert(degree_in_S[SC[i]] + K >= S_end);
//             for (ui i = S_end; i < C_end; i++)
//                 assert(degree_in_S[SC[i]] + K > S_end);
//             for (ui i = 0; i < C_end; i++)
//                 assert(level_id[SC[i]] == n);
// #endif
//             kplex_search(S_end, C_end, level + 1, -1);
//         }
//         restore_C(S_end, C_end, old_C_end, level);
//     }

//     ui choose_branch_vertex(ui S_end, ui C_end)
//     {
//         return SC[S_end];
//         // ui u = n, min_degree_in_S = n;
//         // for(ui i = 0; i < C_end; i++) {
//         //     ui v = SC[i];
//         //     // if(degree[v] + K >= C_end) continue;
//         //     if(degree_in_S[v] < min_degree_in_S) {
//         //         u = v;
//         //         if(degree[v] + K >= C_end) {
//         //             min_degree_in_S = degree_in_S[v];
//         //         }
//         //     }
//         // }
//         // assert(u != n);
//         // if(min_degree_in_S == n) {
//         //     total += 100000;
//         // }
//         // if(SC_rid[u] < S_end) {
//         //     int *t_matrix = matrix+u*n;
//         //     u = n;
//         //     ui max_degree = 0;
//         //     for(ui i = S_end; i < C_end; i++) if(!t_matrix[SC[i]]) {
//         //         if(degree[SC[i]] > max_degree) {
//         //             max_degree = degree[SC[i]];
//         //             u = SC[i];
//         //         }
//         //     }
//         // }
//         // if(u == n) u = SC[S_end];
//         // assert(u != n);
//         // return u;
//     }

//     bool move_u_to_S(ui& S_end, ui& C_end, ui level, ui u)
//     {
//         swap_pos(SC_rid[u], S_end++);
//         assert(u == SC[S_end - 1]);

//         ui neighbors_n = 0, nonneighbors_n = 0;
//         for (ui i = 0; i < C_end; i++) if (SC[i] != u) {
//             if (matrix[u * n + SC[i]]) neighbors[neighbors_n++] = SC[i];
//             else nonneighbors[nonneighbors_n++] = SC[i];
//         }
//         assert(neighbors_n + nonneighbors_n == C_end - 1);
//         for (ui i = 0; i < neighbors_n; i++) degree_in_S[neighbors[i]]++;

//         assert(Qv.empty());
//         // reduce
//         if (degree_in_S[u] + K == S_end) {
//             ui i = 0;
//             while (i < nonneighbors_n && SC_rid[nonneighbors[i]] < S_end) i++;
//             for (; i < nonneighbors_n; i++) {
//                 level_id[nonneighbors[i]] = level;
//                 Qv.push(nonneighbors[i]);
//             }
//         }
//         else {
//             ui i = 0;
//             while (i < nonneighbors_n && SC_rid[nonneighbors[i]] < S_end) i++;
//             for (; i < nonneighbors_n; i++) if (degree_in_S[nonneighbors[i]] + K <= S_end) {
//                 level_id[nonneighbors[i]] = level;
//                 Qv.push(nonneighbors[i]);
//             }
//         }

//         for (ui i = 0; i < nonneighbors_n && SC_rid[nonneighbors[i]] < S_end; i++) {
//             if (degree_in_S[nonneighbors[i]] + K == S_end) {
//                 int* t_matrix = matrix + nonneighbors[i] * n;
//                 for (ui j = S_end; j < C_end; j++) if (level_id[SC[j]] == n && !t_matrix[SC[j]]) {
//                     level_id[SC[j]] = level;
//                     Qv.push(SC[j]);
//                 }
//             }
//         }

//         // check balance
//         {
//             ui i = 0;
//             while (i < neighbors_n && SC_rid[neighbors[i]] < S_end) i++;
//             for (; i < neighbors_n; i++) {
//                 ui v = neighbors[i];
//                 if (level_id[neighbors[i]] == n && !check_balance(S_end, neighbors[i])) {
//                     level_id[neighbors[i]] = level;
//                     assert(neighbors[i] != u);
//                     Qv.push(neighbors[i]);
//                 }
//             }
//         }
//         return remove_vertices(S_end, C_end, level);
//     }

//     // old_C_end->C_end
//     bool remove_vertices(ui S_end, ui& C_end, ui level)
//     {
//         while (!Qv.empty()) {
//             // remove u from C
//             ui u = Qv.front(); Qv.pop();
//             assert(level_id[u] == level);
//             assert(SC[SC_rid[u]] == u);
//             assert(SC_rid[u] >= S_end && SC_rid[u] < C_end);
//             swap_pos(SC_rid[u], --C_end);
//             assert(u == SC[C_end]);

//             int* t_matrix = matrix + u * n;
//             bool terminate = false;

//             ui neighbors_n = 0;
//             for (ui i = 0; i < C_end; i++) if (t_matrix[SC[i]]) {
//                 ui v = SC[i];
//                 neighbors[neighbors_n++] = SC[i];
//                 degree[v]--;
//                 if (degree[v] + K <= best_solution_size) {
//                     if (i < S_end) terminate = true;
//                     else if (level_id[v] == n) {
//                         level_id[v] = level;
//                         Qv.push(v);
//                     }
//                 }
//             }
//             if (terminate) return true;
//         }
//         return false;
//     }

//     // C_end->old_C_end
//     void restore_C(ui S_end, ui& C_end, ui old_C_end, ui level)
//     {
//         while (!Qv.empty()) {
//             ui u = Qv.front(); Qv.pop();
//             assert(level_id[u] == level && SC_rid[u] < C_end);
//             level_id[u] = n;
//         }
//         while (C_end < old_C_end) {
//             ui u = SC[C_end];
//             assert(level_id[u] == level && SC_rid[u] == C_end);
//             level_id[u] = n;
//             for (ui i = 0; i < C_end; i++) if (matrix[u * n + SC[i]]) degree[SC[i]]++;
//             C_end++;
//         }
//     }

//     // u: S->C
//     void move_u_to_C(ui& S_end, ui C_end, ui level)
//     {
//         assert(S_end);
//         ui u = SC[--S_end];
//         ui neighbors_n = 0;
//         int* t_matrix = matrix + u * n;
//         for (ui i = 0; i < C_end; i++) if (t_matrix[SC[i]]) degree_in_S[SC[i]]--;
//     }

//     bool collect_removable_vertices(ui S_end, ui& C_end, ui level)
//     {
//         assert(Qv.empty());
//         for (ui i = 0; i < S_end; i++) if (degree[SC[i]] + K <= best_solution_size) return true;
//         for (ui i = S_end; i < C_end; i++) if (degree[SC[i]] + K <= best_solution_size) {
//             level_id[SC[i]] = level;
//             Qv.push(SC[i]);
//         }
//         return remove_vertices(S_end, C_end, level);
//     }

//     bool remove_u_from_C(ui S_end, ui& C_end, ui u, ui level)
//     {
//         // assert(Qv.empty());
//         // if (u != SC[S_end]) return false;
//         // swap_pos(S_end, --C_end);
//         // assert(u == SC[C_end]);
//         // assert(level_id[u] == n);
//         // level_id[u] = level;

//         // ui neighbors_n = 0, nonneighbors_n = 0;
//         // for (ui i = 0; i < C_end; i++) {
//         //     assert(SC[i] != u);
//         //     if (matrix[u * n + SC[i]])
//         //         neighbors[neighbors_n++] = SC[i];
//         //     else
//         //         nonneighbors[nonneighbors_n++] = SC[i];
//         // }
//         // assert(neighbors_n + nonneighbors_n == C_end);
//         // for (ui i = 0; i < neighbors_n; i++)
//         //     degree[neighbors[i]]--;

//         // for (ui i = 0; i < neighbors_n; i++) {
//         //     ui v = neighbors[i];
//         //     if (degree[v] + K <= best_solution_size) {
//         //         if (i < S_end)
//         //             return true;
//         //         else {
//         //             assert(level_id[v] == n);
//         //             level_id[v] = level;
//         //             Qv.push(v);
//         //         }
//         //     }
//         // }
//         // return remove_vertices(S_end, C_end, level);
//         assert(Qv.empty());
//         if (u != SC[S_end]) return false;
//         assert(level_id[u] == n);
//         Qv.push(u);
//         level_id[u] = level;
//         return remove_vertices(S_end, C_end, level);
//     }

//     void swap_pos(ui i, ui j)
//     {
//         swap(SC[i], SC[j]);
//         SC_rid[SC[i]] = i;
//         SC_rid[SC[j]] = j;
//     }

//     ui upper_bound(ui S_end, ui C_end) {
//         ui ub1 = upper_bound_based_partition_1(S_end, C_end);
//         ui ub2 = upper_bound_based_partition_2(S_end, C_end);
//         ui ub3 = upper_bound_based_partition_3(S_end, C_end);
//         // assert(ub1 <= ub2);
//         // assert(ub2 == ub3);

//         return min(ub1, min(ub2, ub3));
//     }

//     ui upper_bound_based_partition_1(ui S_end, ui C_end)
//     {
//         ui ub = 0, pi0 = C_end - S_end;
//         ui* count = neighbors;
//         ui* C_missing_edges = nonneighbors;
//         memset(count, 0, sizeof(ui) * n);
//         memset(vis, 0, sizeof(ui) * n);
//         for (ui i = 0; i < S_end; i++) {
//             assert(K >= S_end - degree_in_S[SC[i]]);
//             C_missing_edges[i] = K - (S_end - degree_in_S[SC[i]]);
//             for (ui j = S_end; j < C_end; j++)
//                 if (!matrix[SC[i] * n + SC[j]])
//                     count[i]++;
//             // printf("count[%d] = %d\n", SC[i], count[i]);
//         }
//         for (ui i = 0; i < S_end; i++) {
//             if (count[i] <= C_missing_edges[i]) {
//                 vis[i] = 1;
//             }
//         }
//         while (1) {
//             // find max dise
//             int uid = -1;
//             double max_dise = 0;
//             for (ui i = 0; i < S_end; i++)
//                 if (!vis[i]) {
//                     double dise = (C_missing_edges[i] == 0) ? 0 : (1.0 * count[i] / C_missing_edges[i]);
//                     if (max_dise < dise) {
//                         max_dise = dise;
//                         uid = i;
//                     }
//                 }
//             if (uid == -1 || max_dise <= 1)
//                 break;
//             vis[uid] = 1;
//             ui u = SC[uid];
//             // printf("u = %d, C_missing_edges = %d, count = %d\n", u, C_missing_edges[uid], count[uid]);
//             // printf("u = %d, contribution = %d, max_dise = %.6lf\n", SC[uid], min(C_missing_edges[uid], count[uid]), max_dise);
//             ub = ub + min(C_missing_edges[uid], count[uid]);

//             for (ui j = S_end; j < C_end; j++) {
//                 if (!vis[j] && !matrix[u * n + SC[j]]) {
//                     vis[j] = 1;
//                     for (ui k = 0; k < S_end; k++) {
//                         if (!vis[k] && !matrix[SC[k] * n + SC[j]]) {
//                             assert(count[k] > 0);
//                             count[k]--;
//                         }
//                     }
//                     pi0--;
//                 }
//             }
//         }

//         int total = 0;
//         for (ui i = S_end; i < C_end; i++) {
//             if (!vis[i])
//                 total++;
//         }
//         // printf("C - total = %d\n", C_end - total);

//         ui* ub2_vertices = neighbors;
//         ui ub2_vertices_n = 0;
//         ui* colors = nonneighbors;
//         memset(colors, 0, sizeof(ui) * n);
//         for (ui i = S_end; i < C_end; i++) {
//             if (!vis[i])
//                 ub2_vertices[++ub2_vertices_n] = SC[i];
//         }

//         // printf("S_end = %d, pi0 = %d\n", S_end, pi0);
//         ub = S_end + pi0 + ub;
//         // printf("S_end = %u, C_end = %u, ub = %u, kplex = %u, total = %d\n", S_end, C_end, ub, best_solution_size, total);
//         assert(ub <= C_end);
//         return ub;
//     }

//     ui upper_bound_based_partition_2(ui S_end, ui C_end)
//     {
//         ui ub = 0, pi0 = C_end - S_end;
//         ui* count = neighbors;
//         ui* C_missing_edges = nonneighbors;
//         memset(count, 0, sizeof(ui) * n);
//         memset(vis, 0, sizeof(ui) * n);
//         for (ui i = 0; i < S_end; i++) {
//             C_missing_edges[i] = K - (S_end - degree_in_S[SC[i]]);
//             for (ui j = S_end; j < C_end; j++)
//                 if (!matrix[SC[i] * n + SC[j]])
//                     count[i]++;
//             // printf("count[%d] = %d\n", SC[i], count[i]);
//         }
//         ui co = 0;
//         for (ui i = 0; i < S_end; i++) {
//             // find max dise
//             int uid = i;
//             double max_dise = 0;
//             vis[uid] = 1;
//             ui u = SC[uid];
//             // printf("u = %d, C_missing_edges = %d, count = %d\n", u, C_missing_edges[uid], count[uid]);
//             // printf("u = %d, contribution = %d, max_dise = %.6lf\n", SC[uid], min(C_missing_edges[uid], count[uid]), max_dise);
//             ub = ub + min(C_missing_edges[uid], count[uid]);

//             for (ui j = S_end; j < C_end; j++) {
//                 if (!vis[j] && !matrix[u * n + SC[j]]) {
//                     vis[j] = 1;
//                     for (ui k = 0; k < S_end; k++) {
//                         if (!vis[k] && !matrix[SC[k] * n + SC[j]]) {
//                             assert(count[k] > 0);
//                             count[k]--;
//                         }
//                     }
//                     pi0--;
//                 }
//             }
//         }
//         // printf("S_end = %d, pi0 = %d\n", S_end, pi0);
//         ub = S_end + pi0 + ub;
//         // printf("S_end = %u, C_end = %u, ub = %u, kplex = %u, total = %d\n", S_end, C_end, ub, best_solution_size, total);
//         assert(ub <= C_end);
//         return ub;
//     }

//     ui upper_bound_based_partition_3(ui S_end, ui C_end)
//     {
//         ui ub = 0, pi0 = C_end - S_end;
//         ui* C_missing_edges = nonneighbors;
//         memset(vis, 0, sizeof(ui) * n);

//         for (ui i = 0; i < S_end; i++) {
//             C_missing_edges[i] = S_end - degree_in_S[SC[i]] - 1;
//         }
//         for (ui i = 0; i < S_end; i++) {
//             ui u = SC[i];
//             ui pii = 0;
//             int* t_matrix = matrix + u * n;
//             for (ui j = S_end; j < C_end; j++) {
//                 if (!vis[j] && !t_matrix[SC[j]]) {
//                     vis[j] = 1;
//                     pii++;
//                     pi0--;
//                 }
//             }
//             ub = ub + min(K - 1 - C_missing_edges[i], pii);
//         }
//         ub = S_end + pi0 + ub;
//         // printf("S_end = %u, C_end = %u, ub = %u, kplex = %u, total = %d\n", S_end, C_end, ub, best_solution_size, total);
//         assert(ub <= n);
//         return ub;
//     }

//     bool check_balance(ui S_end, ui u)
//     {
//         ui neighbors_n = 0;
//         ui* nei = vis;
//         int* t_matrix = matrix + u * n;
//         for (ui i = 0; i < S_end; i++) if (t_matrix[SC[i]]) nei[neighbors_n++] = SC[i];
//         for (ui i = 0; i < neighbors_n; i++) for (ui j = i + 1; j < neighbors_n; j++) {
//             ui v = nei[i], w = nei[j];
//             if (!matrix[v * n + w]) continue;
//             int tri = matrix[v * n + w] + t_matrix[v] + t_matrix[w];
//             assert(tri == 3 || tri == 1 || tri == -1 || tri == -3);
//             if (tri == 1 || tri == -3) return false;
//         }
//         return true;
//     }
// };

// #endif /* _SIGNED_KPLEX_ */