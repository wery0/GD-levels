#pragma GCC optimize("Ofast")
#include "bits/stdc++.h"
#define timeStamp() std::chrono::steady_clock::now()
#define duration_milli(a) chrono::duration_cast<chrono::milliseconds>(a).count()
using namespace std;
mt19937 rnd(timeStamp().time_since_epoch().count());

struct action {
    char c = 0;

    action() = default;
    action(char c): c(c) {
        // assert(find(possible_actions.begin(), possible_actions.end(), *this) != possible_actions.end());
    }
    bool operator==(const action& rhs) const {return c == rhs.c;}

    void inverse() {c = c == 'A' ? 'D' : c == 'D' ? 'A' : c == 'L' ? 'R' : c == 'R' ? 'L' : c == '+' ? '-' : '+';}
    action get_inverse() const {action res = *this; res.inverse(); return res;}
    bool is_doable_from(int from) {
        assert(from == 1 || from == 2);
        return c == '+' ? from == 1 : c == '-' ? from == 2 : true;
    }

    string to_string() const {return string(1, c);}
    friend ostream& operator<<(ostream& out, const action& rhs) {return out << rhs.to_string();}
};

vector<action> possible_actions;

const uint64_t M = 120632132875995;
uint64_t cache[16] = {0};
struct state {
    uint64_t mask = 0;
    uint8_t w = 0;

    struct hsh {
        size_t operator()(const state& s) const {
            static size_t n[8] = {rnd(), rnd(), rnd(), rnd(), rnd(), rnd(), rnd(), rnd()};
            return s.mask ^ s.w;
            // return custom_hash().xs(s.mask) ^ n[s.w];
        }
    };

    state() = default;
    state(uint64_t mask, uint8_t w): mask(mask), w(w) {}
    state(vector<string> s, uint8_t w): w(w) {
        assert(s.size() == 4);
        for (int i = 0; i < s.size(); ++i) {
            assert(s[i].size() == 4);
            for (int j = 0; j < s[i].size(); ++j) {
                assert(isdigit(s[i][j]));
                uint64_t d = s[i][j] - '0';
                assert(0 <= d && d <= 3);
                mask |= d << ((i * 4 + j) * 3);
            }
        }
        assert((mask & M) == mask);
    }

    bool operator==(const state& rhs) const {return mask == rhs.mask && w == rhs.w;}
    bool operator<(const state& rhs) const {return mask < rhs.mask || mask == rhs.mask && w < rhs.w;}

    bool can_apply_action(const action& a) const {
        return true;
    }

    void apply_action(const action& a) {
        assert(can_apply_action(a));
        if (!cache[w]) {
            int x = w / 4, y = w % 4;
            for (int i = x < 2 ? 0 : x - 1; i < (x < 3 ? x + 2 : x + 1); ++i) {
                for (int j = y < 2 ? 0 : y - 1; j < (y < 3 ? y + 2 : y + 1); ++j) {
                    cache[w] += 1ull << ((i * 4 + j) * 3);
                }
            }
        }
        if (a.c == 'A') w += (w & 3) == 0 ? 3 : -1;
        else if (a.c == 'D') w += (w & 3) == 3 ? -3 : 1;
        else if (a.c == 'L') w += (w & 12) == 0 ? 12 : -4;
        else if (a.c == 'R') w += (w & 12) == 12 ? -12 : 4;
        else if (a.c == '+') {
            mask += cache[w];
            mask &= M;
        } else if (a.c == '-') {
            mask += cache[w] * 3;
            mask &= M;
        } else {
            cout << "UNKNOWN ACTION: " << a << endl;
            assert(0);
        }
        assert(0 <= w && w < 16);
    }

    void undo_action(const action& a) {
        apply_action(a.get_inverse());
    }

    vector<pair<action, state>> get_adjacent_states(int from) const {
        vector<pair<action, state>> res;
        auto process = [&](action a) {
            if (a.is_doable_from(from) && can_apply_action(a)) {
                state ns = *this;
                ns.apply_action(a);
                res.emplace_back(a, ns);
            }
        };
        for (action a : possible_actions) {
            process(a);
        }
        return res;
    }

    int calc_distance_from(const state& rhs) const {
        int res = w < rhs.w ? rhs.w - w : w - rhs.w;
        res = 0;
        for (int i = 0; i < 16; ++i) {
            int x = mask >> (i * 3) & 3;
            int y = rhs.mask >> (i * 3) & 3;
            int d = abs(x - y);
            res += min(d, 4 - d);
        }
        return res;
    }

    string to_string() const {
        string res;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                res += '0' + ((mask >> ((i * 4 + j) * 3)) & 3);
            }
            res += '\n';
        }
        res += std::to_string(int(w) / 4) + " " + std::to_string(int(w) % 4);
        return res;
    }

    friend ostream& operator<<(ostream& out, const state& rhs) {return out << rhs.to_string();}
};

/*
  If A is closer to 0, then path is searching slower, but it is shorter. A = 0 will find the shortest possible path.
  If A is closer to 1, then path is searching faster, but it is longer.
*/
const double A = 0.7;
uniform_real_distribution<double> gen(0, 1);


/*
  Example of input:
  3333
  0121
  3333
  0223
*/
int main() {
    {
        string s1 = "ADLR+-";
        for (int i = 0; i < s1.size(); ++i) {
            possible_actions.push_back(action(s1[i]));
        }
    }
    vector<string> input(4);
    for (auto& s : input) cin >> s;

    vector<state> start_states = {state(input, 0)};
    vector<state> end_states;
    for (int w = 0; w < 16; ++w) {
        end_states.emplace_back(0, w);
    }
    assert(start_states.size() && end_states.size());
    cout << "START STATES:\n"; for (auto s : start_states) cout << s << '\n'; cout << '\n';
    cout << "END STATES:\n"; for (auto s : end_states) cout << s << '\n'; cout << endl;

    unordered_map<state, pair<uint8_t, action>, state::hsh> us(1 << 25);
    deque<state> dq(start_states.begin(), start_states.end());
    dq.insert(dq.end(), end_states.begin(), end_states.end());
    bool found = false;
    state ans_start_state, ans_end_state;
    vector<action> ans_actions;
    for (auto s : start_states) us[s].first = 1;
    for (auto s : end_states) {
        if (us[s].first == 1) found = true, ans_start_state = s, ans_end_state = s, us[s].first = 3;
        else us[s].first = 2;
    }
    uint64_t visited_states = 0, looked_again_states = 0;
    const uint64_t STEP = 5e6;
    auto start_time = timeStamp();
    for (uint64_t d = 0, eso = dq.size(), thr = 0; dq.size() && !found;) {
        if (eso == 0) eso = dq.size(), ++d; --eso;
        if (us.size() > thr) {cout << "CACHE SIZE: " << thr / 1e6 << "M, time " << duration_milli(timeStamp() - start_time) / 1000.0 << "s" << endl; thr += STEP;}
        state cur_state = dq.front(); dq.pop_front();
        auto [from, last_action] = us[cur_state];
        auto states_to_explore = cur_state.get_adjacent_states(from);
        if (d > 1 && gen(rnd) < A) {
            auto cmp = [&](const auto& l, const auto& r) {
                const auto& m = from == 1 ? end_states : start_states;
                auto dl = l.second.calc_distance_from(m[0]);
                auto dr = r.second.calc_distance_from(m[0]);
                return dl < dr;
            };
            nth_element(states_to_explore.begin(), states_to_explore.begin() + 1, states_to_explore.end(), cmp);
            states_to_explore.resize(1);
        }
        for (auto [action, next_state] : states_to_explore) {
            auto& [nm, npa] = us[next_state];
            if (nm != from) {
                if (nm == 0) {
                    nm = from, npa = action;
                    dq.emplace_back(next_state);
                    visited_states += 1;
                    // cout << cur_state << " -(" << int(from) << ' ' << action << ")> " << next_state << endl;
                } else {
                    assert(nm == 3 - from);
                    cout << "FOUND!" << endl;
                    found = true;
                    ans_start_state = from == 1 ? cur_state : next_state;
                    while (find(start_states.begin(), start_states.end(), ans_start_state) == start_states.end()) {
                        assert(us.count(ans_start_state) && "Check your undoing!");
                        ans_actions.push_back(us[ans_start_state].second);
                        ans_start_state.undo_action(us[ans_start_state].second);
                    }
                    reverse(ans_actions.begin(), ans_actions.end());
                    ans_actions.push_back(from == 1 ? action : action.get_inverse());
                    ans_end_state = from == 1 ? next_state : cur_state;
                    while (find(end_states.begin(), end_states.end(), ans_end_state) == end_states.end()) {
                        assert(us.count(ans_end_state) && "Check your undoing!");
                        ans_actions.push_back(us[ans_end_state].second.get_inverse());
                        ans_end_state.undo_action(us[ans_end_state].second);
                    }
                    break;
                }
            } else {
                ++looked_again_states;
            }
        }
    }
    cout << "VISITED STATES: " << visited_states << '\n';
    cout << "CACHE SIZE: " << us.size() << '\n';
    cout << "LOOKED AGAIN STATES: " << looked_again_states << '\n';
    double t = duration_milli(timeStamp() - start_time) / 1000.0;
    cout << "TIME: " << t << "s" << '\n';
    cout << "SPEED: " << int((visited_states + looked_again_states) / t) / 10000 * 10000 / 1e6 << "M states / second\n\n";
    if (found) {
        cout << "FOUND!\n";
        cout << "START STATE: " << ans_start_state << '\n';
        cout << "END STATE: " << ans_end_state << '\n';
        cout << "ACTIONS (" << ans_actions.size() << "): "; for (action a : ans_actions) cout << a; cout << '\n';
    } else {
        cout << "NOT FOUND!";
    }
    cout << endl;
}
