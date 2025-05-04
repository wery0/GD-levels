#include "bits/stdc++.h"
using namespace std;

string D = "URDL";
const int dx[] = {-1, 0, 1, 0};
const int dy[] = {0, 1, 0, -1};

/*
  0 - empty
  1 - stone
  2 - start
  3 - finish
  4 - bomb
  5 - stand

  Input for 2nd level:
  6 6
  300000
  000000
  000001
  100000
  000000
  000002
*/
int main() {
    int a, b; cin >> a >> b;
    vector<string> m(a);
    for (auto& s : m) cin >> s, assert(s.size() == b);
    for (auto& i : m) for (auto& j : i) j -= 48;
    map<vector<string>, string> mp;
    vector<string> need = m;
    for (auto& i : need) for (auto& j : i) j = j == 2 || j == 3 ? 0 : j;
    deque<vector<string>> dq = {m};
    mp[m] = "";
    while (dq.size())  {
        auto cur = dq[0]; dq.pop_front();
        for (int d = 0; d < 4; ++d) {
            auto nw = cur;
            int lox = 0;
            vector stand(a, vector<int>(b, 0));
            for (;;) {
                int fl = 0;
                for (int x = 0; x < a; ++x) {
                    for (int y = 0; y < b; ++y) {
                        if (nw[x][y] != 2) continue;
                        if (stand[x][y]) continue;
                        int nx = x + dx[d];
                        int ny = y + dy[d];
                        if (nx < 0 || ny < 0 || nx >= a || ny >= b) continue;
                        if (nw[nx][ny] == 4) {lox = 1; continue;}
                        if (nw[nx][ny] == 1) continue;
                        if (nw[nx][ny] == 2) continue;
                        if (m[nx][ny] == 5) nw[nx][ny] = 5;
                        fl = 1;
                        if (nw[nx][ny] == 0) {
                            swap(nw[x][y], nw[nx][ny]);
                            x = nx, y = ny;
                        } else if (nw[nx][ny] == 3) {
                            nw[x][y] = 0;
                            nw[nx][ny] = 0;
                        } else {
                            assert(nw[nx][ny] == 5);
                            nw[x][y] = 0;
                            nw[nx][ny] = 2;
                            stand[nx][ny] = 1;
                        }
                    }
                }
                if (!fl) break;
            }
            for (int q = 0; q < a; ++q) {
                for (int w = 0; w < b; ++w) {
                    if (m[q][w] == 5) nw[q][w] = nw[q][w] == 0 ? 5 : nw[q][w];
                }
            }
            if (!lox && !mp.count(nw)) mp[nw] = mp[cur] + D[d], dq.push_back(nw);
        }
    }
    if (mp.count(need)) cout << mp[need];
    else cout << "NO SOLUTION";
}
