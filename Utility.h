#ifndef _UTILITY_H_
#define _UTILITY_H_

#define NDEBUG
// #define NPRINT
#include <bits/stdc++.h>
using namespace std;

typedef unsigned int ui;

using ui = unsigned int;
using ept = unsigned long;
using pii = pair<ui, ui>;

static long long dfs_cnt = 0;
static long long dfs_cnt_1 = 0;
static long long dfs_cnt_2 = 0;
static long long dfs_cnt_after_prune = 0;

struct Edge
{
    int a, b, c;
    Edge(int a, int b, int c) : a(a), b(b), c(c) {}
};

static string integer_to_string(long long number)
{
    std::vector<ui> sequence;
    if (number == 0)
        sequence.push_back(0);
    while (number > 0) {
        sequence.push_back(number % 1000);
        number /= 1000;
    }

    char buf[5];
    std::string res;
    for (unsigned int i = sequence.size(); i > 0; i--) {
        if (i == sequence.size())
            sprintf(buf, "%u", sequence[i - 1]);
        else
            sprintf(buf, ",%03u", sequence[i - 1]);
        res += std::string(buf);
    }
    return res;
}

template <typename T>
int find(T *st, T *ed, T x)
{
    return lower_bound(st, ed, x) - st;
}

#endif /* _UTILITY_H_ */