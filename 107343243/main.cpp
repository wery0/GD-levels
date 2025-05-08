#pragma GCC optimize("Ofast")
#include "bits/stdc++.h"
using namespace std;

map<char, int> mp = {
    {'2', 2},
    {'3', 3},
    {'4', 4},
    {'5', 5},
    {'6', 6},
    {'7', 7},
    {'8', 8},
    {'9', 9},
    {'D', 10},
    {'A', 1},
    {'J', 10},
    {'K', 10},
    {'Q', 10},
};
int a, b;
vector<int> pwa;
vector<string> board;
string ranks = "KQJD98765432A";
int kek[128], lol[128];
vector<int> dp;
vector<uint64_t> pr;
vector<vector<int>> dp2;
vector<vector<uint64_t>> pr2;

vector<int> get_state(int mask) {
    vector<int> state(b);
    for (int u = mask, i = 0; i < b; ++i, u /= (a + 1)) {
        state[i] = u % (a + 1);
    }
    return state;
}

int get_mask(vector<int> state) {
    int u = 0;
    for (int i = 0; i < b; ++i) {
        u += pwa[i] * state[i];
    }
    return u;
}

int eval_move(const string& s) {
    assert(s.size());
    int o = s[0] == 'J' ? 2 : 0;
    {
        for (int l = 0; l + 2 < s.size(); ++l) {
            int u = 0, br = -1;
            for (int r = l; r < s.size(); ++r) {
                if (u & (1 << kek[s[r]])) break;
                u |= 1 << kek[s[r]];
                if (r - l + 1 >= 3 && __lg(u) - __builtin_ctz(u) == r - l) {
                    o += r - l + 1;
                    br = r;
                }
            }
            if (br != -1) {
                //o += br - l + 1;
                l = br;
            }
        }
    }
    int d = 0;
    for (int pr = 0, sum = 0; auto c : s) {
        sum += lol[c];
        if (sum == 15) o += 2;
        if (sum == 31) o += 2;
        if (c == pr) {
            ++d;
            assert(d <= 4);
        } else {
            o += d == 2 ? 2 : d == 3 ? 8 : d == 4 ? 20 : 0;
            d = 1;
        }
        pr = c;
    }
    o += d == 2 ? 2 : d == 3 ? 8 : d == 4 ? 20 : 0;
    return o;
}

int get_new_mask(int mask, uint64_t move) {
    int o = mask, sz = move & 255; move >>= 8;
    assert(sz <= 13);
    for (int j = 0; j < sz; ++j) {
        o -= pwa[move & 3];
        move >>= 2;
    }
    return o;
}

string calc_smove_by_mask_and_move(int mask, uint64_t move) {
    int sz = move & 255; move >>= 8;
    assert(sz <= 13);
    string res;
    res.reserve(sz);
    for (int j = 0; j < sz; ++j) {
        int i = move & 3; move >>= 2;
        int d = mask / pwa[i] % (a + 1);
        assert(d > 0);
        res += board[--d][i];
        mask -= pwa[i];
    }
    return res;
}

int calc_clicks_of_move(uint64_t move) {
    int pos = -1, sz = move & 255; move >>= 8;
    int res = 0;
    for (int j = 0; j <= sz; ++j) {
        int i = j == sz ? -1 : move & 3; move >>= 2;
        res += abs(pos - i);
        res += 1;
        pos = i;
    }
    return res;
}

vector<uint64_t> gen_all_possible_moves(int mask) {
    static vector<uint64_t> res;
    static unordered_set<uint64_t> cache;
    res.clear();
    cache.clear();
    function<void(int, int, uint64_t, uint64_t)> go = [&](int mask, int sum, uint64_t move, uint64_t cards) {
        if (cache.count(cards * pwa[b] + mask)) return;
        cache.insert(cards * pwa[b] + mask);
        int fl = 0, dep = move & 255;
        for (uint64_t i = 0; i < b; ++i) {
            int d = mask / pwa[i] % (a + 1);
            if (d == 0 || sum + lol[board[d - 1][i]] > 31) continue;
            fl = 1;
            char card = board[d - 1][i];
            go(mask - pwa[i], sum + lol[card], (move + 1) | (i << (8 + dep * 2)), cards * 13 + kek[card]);
        }
        if (!fl) {
            res.push_back(move);
        }
    };
    go(mask, 0, 0, 0);
    return res;
}

int calc(int mask) {
    if (mask == 0) return 0;
    if (pr[mask]) return dp[mask];
    auto moves = gen_all_possible_moves(mask);
    // cout << mask << ": " << moves.size() << endl;
    for (int i = 0; i < moves.size(); ++i) {
        int new_mask = get_new_mask(mask, moves[i]);
        int val = eval_move(calc_smove_by_mask_and_move(mask, moves[i]));
        auto nw = calc(new_mask) + val;
        if (nw > dp[mask]) {
            dp[mask] = nw;
            pr[mask] = moves[i];
        }
    }
    return dp[mask];
}

//dp2[mask][score] = min number of moves to get score >= 61 and clear the board if now board state if mask and score is score
uniform_int_distribution<int> gen(0, 100000);
int calc2(int mask, int score) {
    // const int G = 1e6+6;
    // static int nm[G], ev[G], clk[G];
    if (mask == 0) return score >= 61 ? 0 : -1;
    if (pr2[mask][score]) return dp2[mask][score];
    auto moves = gen_all_possible_moves(mask);
    int ms = moves.size();
    vector<int> nm(ms), ev(ms), clk(ms);
    for (int i = 0; i < ms; ++i) {
        nm[i] = get_new_mask(mask, moves[i]);
        assert(nm[i] < mask);
        ev[i] = eval_move(calc_smove_by_mask_and_move(mask, moves[i]));
        clk[i] = calc_clicks_of_move(moves[i]);
    }
    fill(pr2[mask].begin(), pr2[mask].end(), -1);
    const int h = dp2[0].size();
    for (int i = 0; i < ms; ++i) {
        for (int j = 0; j < h; ++j) {
            assert(moves[i] & 255);
            int new_mask = nm[i];
            int val = ev[i];
            int cost = clk[i];
            int tyt = calc2(new_mask, min(61, j + val));
            if (tyt == -1) continue;
            if (dp2[mask][j] == -1 || tyt + cost < dp2[mask][j]) {
                dp2[mask][j] = tyt + cost;
                pr2[mask][j] = moves[i];
            }
        }
    }
    return dp2[mask][score];
}

vector<int> mask_to_move(uint64_t mask) {
    int sz = mask & 255; mask >>= 8;
    assert(sz <= 13);
    vector<int> res(sz);
    for (int i = 0; i < sz; ++i) {
        res[i] = mask & 3; mask >>= 2;
    }
    return res;
}

string move_to_clicks(uint64_t move) {
    int sz = move & 255; move >>= 8;
    int pos = -1;
    string res;
    for (int j = 0; j <= sz; ++j) {
        int i = j == sz ? -1 : move & 3; move >>= 2;
        int u = abs(i - pos);
        char c = i > pos ? 'R' : 'L';
        for (int i = 0; i < u; ++i) res += c;
        res += '+';
        pos = i;
    }
    return res;
}

void restore_answer(int mask, int start_score, int mode) {
    auto state = get_state(mask);
    cout << "STARTING BOARD: " << endl;
    for (int q = 0; q < a; ++q) {
        for (int w = 0; w < b; ++w) {
            if (state[w] > q) cout << board[q][w];
            else cout << "*";
        }
        cout << endl;
    }
    cout << endl;
    cout << "STARTING SCORE: " << start_score << endl;
    string clicks;
    int score = start_score;
    int clicks_cnt = 0;
    for (; mask;) {
        auto best_move = mode == 1 ? pr[mask] : pr2[mask][min(61, score)];
        auto best_smove = calc_smove_by_mask_and_move(mask, best_move);
        int val = eval_move(best_smove);
        clicks_cnt += calc_clicks_of_move(best_move);
        clicks += move_to_clicks(best_move) + "\n";
        score += val;
        cout << "     MOVE: "; for (int i : mask_to_move(best_move)) cout << i << " "; cout << endl;
        cout << "    SMOVE: " << best_smove << endl;
        cout << " MOVE VAL: " << val << endl;
        cout << "CUR SCORE: " << score << endl;
        cout << endl;
        int new_mask = get_new_mask(mask, best_move);
        mask = new_mask;
    }
    cout << "       SCORE: " << score << endl;
    cout << "CLICKS_COUNT: " << clicks_cnt << endl;
    cout << "      CLICKS:\n" << clicks << endl;
}

/*
  D stands for 10
  MODE == 1 finds highest possible score and a path to it
  MODE == 2 find shortest possible path to score >= 61

  Example of input:
  13 4
  634A
  D59D
  A2KJ
  735Q
  4726
  JJ6Q
  52KK
  A8K4
  78AD
  8853
  QDQ2
  9J34
  6799
*/
int main() {
    for (int i = 0; char c : ranks) kek[c] = i++;
    for (auto [ch, vl] : mp) lol[ch] = vl;
    assert(eval_move("J588") == 8);
    assert(eval_move("JQAQ") == 4);
    assert(eval_move("DDD") == 8);
    assert(eval_move("JJJ") == 10);
    assert(eval_move("6666") == 20);
    assert(eval_move("5555") == 22);
    assert(eval_move("555547") == 24);
    assert(eval_move("5554444") == 32);
    assert(eval_move("33332222AAAA") == 60);
    assert(eval_move("6547333") == 24);
    assert(eval_move("3245A67") == 27);

    cin >> a >> b;
    board.resize(a);
    for (auto& s : board) cin >> s; 
    {
        map<char, int> dd;
        for (auto s : board) for (auto c : s) ++dd[c];
        int fl = 1;
        for (auto [x, c] : dd) {
            if (c != 4) {
                cout << x << ": " << c << '\n';
                fl = 0;
            }
        }
        assert(fl && "Wrong board");
    }
    pwa.resize(b + 1, 1);
    for (int q = 1; q <= b; ++q) pwa[q] = pwa[q - 1] * (a + 1);
    int A = pwa[b];
    dp.resize(A, -1);
    pr.resize(A);
    dp2.resize(A, vector<int>(62, -1));
    pr2.resize(A, vector<uint64_t>(62));
    dp[0] = 0;
    dp2[0][0] = 0;
    vector<int> start_state(b, a);
    vector<int> used = {0, 0, 0, 0};
    assert(used.size() == b);
    for (int i = 0; i < b; ++i) start_state[i] -= used[i];
    int need = get_mask(start_state);
    
    int MODE = 1;
    assert(MODE == 1 || MODE == 2);
    MODE == 1 ? calc(need) : calc2(need, 0);
    restore_answer(need, 0, MODE);
}
