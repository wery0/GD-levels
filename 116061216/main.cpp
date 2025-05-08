#pragma GCC optimize("Ofast")
#include "bits/stdc++.h"
#define timeStamp() std::chrono::steady_clock::now()
#define duration_milli(a) chrono::duration_cast<chrono::milliseconds>(a).count()
using namespace std;
mt19937 rnd(timeStamp().time_since_epoch().count());

#include <ext/pb_ds/assoc_container.hpp>
using namespace __gnu_pbds;

struct action {
    char c;

    action() = default;
    action(char c): c(c) {assert(tolower(c) == 'l' || tolower(c) == 'r');}
    bool operator==(const action& rhs) const {return c == rhs.c;}

    void inverse() {c = c == 'l' ? 'r' : c == 'r' ? 'l' : c == 'L' ? 'R' : 'L';}
    action get_inverse() const {action res = *this; res.inverse(); return res;}

    string to_string() const {return string(1, c);}
    friend ostream& operator<<(ostream& out, const action& rhs) {return out << rhs.to_string();}
};

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
    state(string s) {
        assert(s.size() == 18);
        for (int i = 0; char c : s) {
            assert(c == '[' || c == ']' || isdigit(c) || 'A' <= c && c <= 'F');
            if (c == '[') w = i;
            else if (c != ']') mask = mask << 4 | (isdigit(c) ? c - 48 : 10 + c - 'A');
            ++i;
        }
    }

    bool operator==(const state& rhs) const {return mask == rhs.mask && w == rhs.w;}
    bool operator<(const state& rhs) const {return mask < rhs.mask || mask == rhs.mask && w < rhs.w;}

    bool can_apply_action(const action& a) const {
        return !(w == 0 && a.c == 'l') && !(w == 8 && a.c == 'r');
    }

    void apply_action(const action& a) {
        assert(can_apply_action(a));
        if (a.c == 'l') --w;
        else if (a.c == 'r') ++w;
        else if (a.c == 'L') {
            int i = 8 - w;
            uint64_t s1 = mask & (((uint64_t(1) << 28) - 1) << (i * 4));
            uint64_t s2 = mask & (((uint64_t(1) << 4) - 1) << (28 + i * 4));
            mask = mask ^ (s1 | s2);
            s1 <<= 4;
            s2 >>= 28;
            mask |= s1 | s2;
        } else if (a.c == 'R') {
            int i = 8 - w;
            uint64_t s1 = mask & (((uint64_t(1) << 4) - 1) << (i * 4));
            uint64_t s2 = mask & (((uint64_t(1) << 28) - 1) << (4 + i * 4));
            mask = mask ^ (s1 | s2);
            s1 <<= 28;
            s2 >>= 4;
            mask |= s1 | s2;
        } else assert(0);
    }

    void undo_action(const action& a) {
        apply_action(a.get_inverse());
    }

    vector<pair<action, state>> get_adjacent_states() const {
        vector<pair<action, state>> res;
        res.reserve(4);
        for (action a : {'l', 'r', 'L', 'R'}) {
            if (can_apply_action(a)) {
                state ns = *this;
                ns.apply_action(a);
                res.emplace_back(a, ns);
            }
        }
        return res;
    }

    vector<pair<action, state>> get_adjacent_states_except_parent(action last_action) const {
        last_action.inverse();
        vector<pair<action, state>> res;
        res.reserve(4);
        for (action a : {'l', 'r', 'L', 'R'}) {
            if (a == last_action) continue;
            if (can_apply_action(a)) {
                state ns = *this;
                ns.apply_action(a);
                res.emplace_back(a, ns);
            }
        }
        return res;
    }

    int calc_distance_from(const state& rhs) const {
        int res = w < rhs.w ? rhs.w - w : w - rhs.w;
        for (int i = 0; i < 16; ++i) {
            int x = mask >> (i * 4) & 15;
            int y = rhs.mask >> (i * 4) & 15;
            res += abs(x - y);
        }
        return res;
    }

    string to_string() const {
        string res;
        for (int i = 15; i >= 0; --i) {
            int u = mask >> (i * 4) & 15;
            res += u < 10 ? 48 + u : 'A' + u - 10;
        }
        res.insert(res.begin() + w + 8, ']');
        res.insert(res.begin() + w, '[');
        return res;
    }

    friend ostream& operator<<(ostream& out, const state& rhs) {return out << rhs.to_string();}
};

//O(nlogn)
int64_t calc_inversions(uint64_t m) {
    uint64_t o = 0;
    uint16_t kek = 0;
    for (int i = 0; i < 16; ++i) {
        int x = m & 15; m >>= 4;
        o += __builtin_popcount(kek >> x);
        kek |= 1 << x;
    }
    return o;
}

/*
  If A is closer to 0, then path is searching slower, but it is shorter. A = 0 will find the shortest possible path.
  If A is closer to 1, then path is searching faster, but it is longer.
*/
const double A = 0;
uniform_real_distribution<double> gen(0, 1);

//State format is like [2AED3675]10C4B98F
/*
  Some states, actually encountered in the level by me (A = 0 results):
  [E460C12D]3B79A5F8 -> 3.1G RAM, 43.3M states, 13s, LLLrLrrLLrrLLLrLLllRlRRrLrrrrR
  [4D073586]CF21A9EB -> 33.4G RAM, 465M states, 263s, RrLLrrRrrrrRrLLLLlLllRllRlllLrrRlLrL
  [E95DC476]8AF2B013 -> 61.4G RAM, 982M states, 918s, RRrRrrLrRrrrLrLLlLlRrLLllRRlRRllLLlLlL
  [643A0FBE]8572C91D -> 61.8G RAM, 988M states, 710s, rrrrRrRrrLrRllLrLLllLlLLlLlLLlLLrLllLL
  [2BF650DC]9EA34718 -> 90.3G RAM, 1.36B states, 1355s, rrLrrrLrrRrLLllRRlLllllLLlLLrLLLrLrRrR
*/
int main() {
    cout << "STARTED" << endl;
    assert(0 <= A && A <= 1);
    string input; cin >> input;
    assert(input.size() == 18 && "WRONG LENGTH");
    bool is_correct = true;
    uint64_t mask = 0;
    for (int po = -1, pc = -1, i = 0; auto c : input) {
        if (c == '[') po = i;
        if (c == ']') pc = i;
        if (po != -1 && pc != -1) assert(pc == po + 9);
        int u = c == '[' ? 16 : c == ']' ? 17 : isdigit(c) ? c - 48 : c + 10 - 'A';
        assert(0 <= u && u < 18 && "WRONG INPUT, unknown character!");
        mask |= 1 << u;
        ++i;
    }
    assert(mask == 262143 && "WRONG INPUT, some character from [0-9][A-F] is missing!");

    vector<state> start_states = {input};
    cout << "START STATES:\n"; for (auto s : start_states) cout << s << '\n'; cout << '\n';

    vector<state> end_states;
    for (int w = 0; w <= 8; ++w) {
        end_states.emplace_back(81985529216486895, w);
    }
    cout << "END STATES:\n"; for (auto s : end_states) cout << s << '\n'; cout << '\n';

    unordered_map<state, pair<uint8_t, action>, state::hsh> us(1 << 30);
    // gp_hash_table<state, pair<uint8_t, action>, state::hsh, equal_to<state>, direct_mask_range_hashing<>, linear_probe_fn<>, hash_standard_resize_policy<
    // hash_exponential_size_policy<>, hash_load_check_resize_trigger<>, true>> us; us.resize(1 << 30);

    deque<state> dq(start_states.begin(), start_states.end());
    dq.insert(dq.end(), end_states.begin(), end_states.end());

    bool found = false;
    state ans_start_state;
    vector<action> ans_actions;
    state ans_end_state;
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
        if (us.size() > thr) {cout << "CACHE SIZE: " << thr << ", time " << duration_milli(timeStamp() - start_time) / 1000.0 << "s" << endl; thr += STEP;}
        state cur_state = dq.front(); dq.pop_front();
        auto [from, last_action] = us[cur_state];
        auto states_to_explore = d ? cur_state.get_adjacent_states_except_parent(last_action) : cur_state.get_adjacent_states();
        if (d > 1 && gen(rnd) < A) {
            auto cmp = [&](const auto& l, const auto& r) {
                return l.second.calc_distance_from(end_states[3]) < r.second.calc_distance_from(end_states[3]);
                // return calc_inversions(l.S.mask) < calc_inversions(r.S.mask);
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
                    // cout << cur_state << " -(" << from << ' ' << action << ")> " << next_state << endl;
                } else {
                    assert(us[next_state].first == 3 - from);
                    found = true;
                    ans_start_state = from == 1 ? cur_state : next_state;
                    while (find(start_states.begin(), start_states.end(), ans_start_state) == start_states.end()) {
                        ans_actions.push_back(us[ans_start_state].second);
                        ans_start_state.undo_action(us[ans_start_state].second);
                    }
                    reverse(ans_actions.begin(), ans_actions.end());
                    ans_actions.push_back(from == 1 ? action : action.get_inverse());
                    ans_end_state = from == 1 ? next_state : cur_state;
                    while (find(end_states.begin(), end_states.end(), ans_end_state) == end_states.end()) {
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
