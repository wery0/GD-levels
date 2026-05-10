#pragma GCC optimize("O3")
#include "bits/stdc++.h"
using namespace std;

const int dx[] = { -1, 0, 1, 0, -1, 1, 1, -1};
const int dy[] = {0, 1, 0, -1, 1, 1, -1, -1};

const int inf = 1e9;

template<typename I> vector<pair<int, int>> get_segs_of_eq_elements(I first, I last, function<bool(const typename iterator_traits<I>::value_type&, const typename iterator_traits<I>::value_type&)> cmp = [](const auto& l, const auto& r) {return l == r;}) {using T = typename iterator_traits<I>::value_type; vector<pair<int, int>> ans; if (first == last) return ans; int l = 0, r = 1; T prev = *first; for (auto cit = next(first); cit != last; ++cit, ++r) {if (!cmp(*cit, prev)) {ans.emplace_back(l, r - 1); l = r;} prev = *cit;} ans.emplace_back(l, r - 1); return ans;}


/*
Input for level 7:
14 21
.*...................
.*...................
.*...................
.*...................
*****................
...2****.............
.4P.*................
.2.3***4...E.........
********.............
*3*.*2.4.............
*.*.*.*..............
*..2.3*.2.*..........
*******...*..........
..........*..........

Input for level 8:
14 21
.2*....*...2*....*...
..*....*....*2...*...
..*....*....*....*...
..*.2..*...2*....*...
*********************
*2...*..***..*....*..
*.P..*..3***.*....*..
*...2*...***3*....*..
*2...*....****....*..
*********************
...*2...*....*....*..
..2*....*...2*...2*..
...*....*....*....*.E
...*....*....*....*..
*/
int main() {
    int a, b; cin >> a >> b;
    vector<string> m(a);
    vector num(a, vector<int>(b, -1));
    int cd = 0;
    int sx = -1, sy = -1;
    int ex = -1, ey = -1;
    for (int q = 0; q < a; ++q) {
        cin >> m[q];
        for (int w = 0; w < b; ++w) {
            if (isdigit(m[q][w])) {
                num[q][w] = cd++;
            } else if (m[q][w] == 'P') {
                sx = q, sy = w;
                m[q][w] = '.';
            } else if (m[q][w] == 'E') {
                ex = q, ey = w;
                m[q][w] = '.';
            } else {
                assert(m[q][w] == '.' || m[q][w] == '*');
            }
        }
    }
    auto can_jump = [&](int mask, int x, int y) -> bool {
        if (m[x][y] == '.') return false;
        if (m[x][y] == '*') return false;
        assert(isdigit(m[x][y]));
        return ~mask >> num[x][y] & 1;
    };
    auto calc = [&](int mask, int cx, int cy, int d) -> array<int, 3> {
        int nx = cx + dx[d];
        int ny = cy + dy[d];
        int nm = mask;
        if (nx < 0 || ny < 0 || nx >= a || ny >= b) return {mask, cx, cy};
        while (can_jump(nm, nx, ny)) {
            nm |= 1 << num[nx][ny];
            int j = m[nx][ny] - 48;
            nx += dx[d] * j;
            ny += dy[d] * j;
            nx = clamp(nx, 0, a - 1); ny = clamp(ny, 0, b - 1);
        }
        if (m[nx][ny] == '*') return {mask, cx, cy};
        return {nm, nx, ny};
    };
    vector dst(1 << cd, vector(a, vector<int>(b, inf)));
    vector pr(1 << cd, vector(a, vector<array<int, 4>>(b, { -1, -1, -1, -1})));
    deque<array<int, 3>> dq;
    dq.push_back({0, sx, sy});
    dst[0][sx][sy] = 0;
    for (; dq.size();) {
        auto [mask, x, y] = dq.front(); dq.pop_front();
        for (int d = 0; d < 4; ++d) {
            auto [fm, fx, fy] = calc(mask, x, y, d);
            if (dst[mask][x][y] + 1 < dst[fm][fx][fy]) {
                dst[fm][fx][fy] = dst[mask][x][y] + 1;
                pr[fm][fx][fy] = {mask, x, y, "URDL"[d]};
                dq.push_back({fm, fx, fy});
            }
        }
    }
    assert(pr[(1 << cd) - 1][ex][ey][3] != -1);
    string o;
    for (int mask = (1 << cd) - 1, x = ex, y = ey; ;) {
        auto [pm, px, py, dir] = pr[mask][x][y];
        if (dir == -1) break;
        o += dir;
        mask = pm, x = px, y = py;
    }
    reverse(o.begin(), o.end());
    cout << o << endl;
    for (auto [l, r] : get_segs_of_eq_elements(o.begin(), o.end())) {
        cout << r - l + 1 << " " << o[l] << endl;
    } 
}
