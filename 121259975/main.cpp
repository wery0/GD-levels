#pragma GCC optimize("Ofast")
// #pragma GCC target("avx,avx2,fma")

#include "bits/stdc++.h"

//#define NDEBUG
#define F first
#define S second
#define vec vector
#define pb push_back
#define pll pair<ll, ll>
#define all(m) m.begin(), m.end()
#define rall(m) m.rbegin(), m.rend()
#define uid uniform_int_distribution
#define timeStamp() std::chrono::steady_clock::now()
#define unify(m) sort(all(m)), m.erase(unique(all(m)), m.end());
#define duration_micro(a) chrono::duration_cast<chrono::microseconds>(a).count()
#define duration_milli(a) chrono::duration_cast<chrono::milliseconds>(a).count()
#define fast cin.tie(0), cout.tie(0), cin.sync_with_stdio(0), cout.sync_with_stdio(0);
using namespace std;
using str = string;
using ll = long long;
using ld = long double;
mt19937 rnd(timeStamp().time_since_epoch().count());
mt19937_64 rndll(timeStamp().time_since_epoch().count());
template<typename T, typename U> bool chmin(T& a, const U& b) {return (T)b < a ? a = b, 1 : 0;}
template<typename T, typename U> bool chmax(T& a, const U& b) {return (T)b > a ? a = b, 1 : 0;}
struct custom_hash {static uint64_t xs(uint64_t x) {x += 0x9e3779b97f4a7c15; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9; x = (x ^ (x >> 27)) * 0x94d049bb133111eb; return x ^ (x >> 31);} template<typename T> size_t operator()(T x) const {static const uint64_t C = timeStamp().time_since_epoch().count(); return xs(hash<T> {}(x) + C);}};
template<typename K> using uset = unordered_set<K, custom_hash>;
template<typename K, typename V> using umap = unordered_map<K, V, custom_hash>;
template<typename T1, typename T2> ostream& operator<<(ostream& out, const pair<T1, T2>& x) {return out << x.F << ' ' << x.S;}
template<typename T1, typename T2> istream& operator>>(istream& in, pair<T1, T2>& x) {return in >> x.F >> x.S;}
template<typename T, size_t N> istream& operator>>(istream& in, array<T, N>& a) {for (auto &x : a) in >> x; return in;}
template<typename T, size_t N> ostream& operator<<(ostream& out, const array<T, N>& a) {for (size_t i = 0; i < a.size(); ++i) {out << a[i]; if (i + 1 < a.size()) out << ' ';} return out;}
template<typename T> istream& operator>>(istream& in, vector<T>& a) {for (auto& x : a) in >> x; return in;}
template<typename T> ostream& operator<<(ostream& out, const vector<T>& a) {for (size_t i = 0; i < a.size(); ++i) {out << a[i]; if (i + 1 < a.size()) out << ' ';} return out;}

array<int, 5> mask_to_arr(int mask) {
    assert(0 <= mask && mask < 7776);
    array<int, 5> m;
    for (int q = 0; q < 5; ++q, mask /= 6) {
        m[q] = mask % 6 + 1;
    }
    return m;
}

int arr_to_mask(array<int, 5> m) {
    int mask = 0;
    for (int q = 4; q >= 0; --q) {
        assert(1 <= m[q] && m[q] <= 6);
        mask = mask * 6 + (m[q] - 1);
    }
    return mask;
}

vector<array<int, 5>> arrs;
int mask_to_num(int mask) {
    assert(0 <= mask && mask < 7776);
    static vector<int> res(7776, -1);
    if (res[mask] == -1) {
        auto arr = mask_to_arr(mask);
        sort(arr.begin(), arr.end());
        auto it = lower_bound(arrs.begin(), arrs.end(), arr);
        assert(*it == arr);
        res[mask] = it - arrs.begin();
    }
    return res[mask];
}

int num_to_mask(int num) {
    return arr_to_mask(arrs[num]);
}

int calc_score(int mask, int i, int num) {
    assert(0 <= i && i < 13);
    assert(0 <= num && num < 252);
    static vector res(13, vector<int>(252, -1));
    auto calc_mxseq = [](array<int, 5> arr) {
        int mxc = 1, c = 1;
        for (int q = 1; q < 5; ++q) {
            if (arr[q] == arr[q - 1] + 1) ++c, mxc = max(mxc, c);
            else if (arr[q] > arr[q - 1] + 1) c = 1;
        }
        return mxc;
    };
    if (i == 11 || res[i][num] == -1) {
        auto arr = arrs[num];
        if (i < 6) res[i][num] = count(arr.begin(), arr.end(), i + 1) * (i + 1);
        else if (i == 6) {
            res[i][num] = arr[0] == arr[2] || arr[1] == arr[3] || arr[2] == arr[4] ? 6 + arr[2] * 3 : 0;
        } else if (i == 7) {
            res[i][num] = arr[0] == arr[3] || arr[1] == arr[4] ? 7 + arr[2] * 3 : 0;
        }
        else if (i == 8) {
            res[i][num] = arr[0] == arr[1] && arr[1] != arr[2] && arr[2] == arr[4] ||
                          arr[0] == arr[2] && arr[2] != arr[3] && arr[3] == arr[4] ? 25 : 0;
        } else if (i == 9) {
            res[i][num] = calc_mxseq(arr) >= 4 ? 30 : 0;
        } else if (i == 10) {
            res[i][num] = calc_mxseq(arr) >= 5 ? 40 : 0;
        } else if (i == 11) {
            int kek = 0;
            for (int j = 0; j < 6; ++j) {
                if (mask >> j & 1) kek += (j + 1) * 3;
            }
            res[i][num] = arr[0] == arr[4] ? 50 + 0 * max(0, 63 - kek) : 0;
        } else if (i == 12) {
            res[i][num] = accumulate(arr.begin(), arr.end(), 0);
        }
    }
    return res[i][num];
}

string bit_to_name(int i) {
    assert(0 <= i && i < 13);
    static vector<string> m = {"Ones", "Twos", "Threes", "Fours", "Fives", "Sixes", "Three of a kind", "Four of a kind", "Full house", "Sm straight", "Lg straight", "Yahtzee", "Chance"};
    return m[i];
}

string describe_move(array<int, 5> m, int b) {
    if (b < 13) {
        return "Take " + bit_to_name(b);
    }
    b -= 100;
    assert(0 <= b && b < 32);
    multiset<int> who_to_lock;
    auto sm = m;
    sort(sm.begin(), sm.end());
    for (int i = 0; i < 5; ++i) {
        if (b >> i & 1) who_to_lock.insert(sm[i]);
    }
    string res = "Lock ";
    for (int i = 0; i < 5; ++i) {
        if (who_to_lock.count(m[i])) {
            who_to_lock.erase(who_to_lock.find(m[i]));
            res += "1";
        } else res += "0";
    }
    return res + " and roll";
}

/*
 Takes about 35 seconds to compute expected values of all positions.
 Enter dices without spaces, i.e. like 63123.
 Lock 01001 means locking second and fifth dices.
 Scores for Three/Four of a kind are approximate since idk what exact formula is used for them
*/
int main() {
    fast;
    auto st = timeStamp();
    function<void(vector<int>)> gen_arrs = [&](vector<int> m) {
        if (m.size() == 5) {
            array<int, 5> u;
            copy(m.begin(), m.end(), u.begin());
            arrs.push_back(u);
            return;
        }
        for (int i = m.size() ? m.back() : 1; i <= 6; ++i) {
            m.push_back(i);
            gen_arrs(m);
            m.pop_back();
        }
    };
    gen_arrs({});
    vector rl(252, vector(32, vector<pair<double, int>>()));
    for (int n = 0; n < 252; ++n) {
        for (int msk = 0; msk < 32; ++msk) {
            auto arr = arrs[n];
            map<int, int> cnt;
            function<void(int)> g = [&](int i) {
                if (i == 5) {
                    ++cnt[mask_to_num(arr_to_mask(arr))];
                    return;
                }
                int L = 1, R = 6;
                if (msk >> i & 1) L = R = arr[i];
                for (int c = L; c <= R; ++c) {
                    arr[i] = c;
                    g(i + 1);
                }
            };
            g(0);
            double tot = powl(6, 5 - __builtin_popcount(msk));
            for (auto [nwnum, c] : cnt) {
                rl[n][msk].emplace_back(c / tot, nwnum);
            }
        }
    }
    cout << "Precalc done in " << duration_milli(timeStamp() - st) << " ms" << endl;
    int moves = 13;
    vector dp(1 << moves, vector(252, array<double, 4> { -1, -1, -1, -1}));
    vector bm(1 << moves, vector(252, array<int, 4> { -1, -1, -1, -1}));
    function<double(int, int, int)> f = [&](int mask, int num, int rol) -> double {
        if (dp[mask][num][rol] >= 0) return dp[mask][num][rol];
        if (mask == dp.size() - 1) return 0;
        double best_take = -1;
        int best_take_bit = -1;
        for (int i = 0; i < moves; ++i) {
            if (mask >> i & 1) continue;
            double tyt = calc_score(mask, i, num) + f(mask | (1 << i), 0, 3);
            if (tyt > best_take) best_take = tyt, best_take_bit = i;
        }
        if (rol == 3) best_take = 0;
        bm[mask][num][rol] = best_take_bit;
        if (rol == 0) return dp[mask][num][rol] = best_take;
        double best_roll = 0;
        int best_roll_mask = 0;
        for (int i = 0; i < (rol == 3 ? 1 : 32); ++i) {
            double tyt = 0;
            for (auto [prob, nwnum] : rl[num][i]) {
                tyt += prob * f(mask, nwnum, rol - 1);
            }
            if (tyt > best_roll) best_roll = tyt, best_roll_mask = i;
        }
        if (best_roll > best_take) bm[mask][num][rol] = 100 + best_roll_mask;
        return dp[mask][num][rol] = max(best_take, best_roll);
    };
    f(0, 0, 3);
    cout << "Calculated every state in " << duration_milli(timeStamp() - st) / 1000.0 << " s" << endl;
    cout << "Expected score: " << dp[0][0][3] << endl;
    int cm = 0, cr = 3;
    array<int, 5> cu{1, 1, 1, 1, 1};
    int score = 0;
    for (; cm < dp.size() - 1;) {
        cout << endl;
        int cn = mask_to_num(arr_to_mask(cu));
        cout << "         Cur state: " << cm << " {" << cu[0] << ", " << cu[1] << ", " << cu[2] << ", " << cu[3] << ", " << cu[4] << "} " << cr << endl;
        cout << "         Cur score: " << score << endl;
        cout << "Cur expected score: " << score + dp[cm][cn][cr] << endl;
        int b = bm[cm][cn][cr];
        cout << "         Best move: " << describe_move(cu, b) << endl;
        if (b < moves) {
            score += calc_score(cm, b, cn);
            cm |= 1 << b;
            cu = {1, 1, 1, 1, 1};
            cr = 3;
        } else {
            --cr;
            assert(cr >= 0);
            cout << "Enter current dices: " << endl;
            string s; cin>>s;
            assert(s.size() == 5);
            for (int i = 0; i < 5; ++i) {
                assert('1' <= s[i] && s[i] <= '6');
                cu[i] = s[i] - '0';
            }
        }
    }
}
