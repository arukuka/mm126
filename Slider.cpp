// C++11
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <memory>
#include <chrono>
#include <cstring>

struct Point {
  int r, c;
};

constexpr double TLE = 9.0;
constexpr int MAX_N = 10;
constexpr int MAX_C = 9;
constexpr std::array<std::array<int, 2>, 4> OFS{{
  {1, 0},
  {0, 1},
  {-1, 0},
  {0, -1}
}};
constexpr char COMMANDS[5] = "RULD";

struct Timer {
  std::chrono::high_resolution_clock::time_point start;
  Timer() : start(std::chrono::high_resolution_clock::now()) {}
  bool TLE() const {
    const auto now = std::chrono::high_resolution_clock::now();
    const auto secs = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count() / 1e3;
    return ::TLE < secs;
  }
};

Timer timer;

int N, C, H;

bool is_out(const int r, const int c) {
  return r < 0 || N <= r
         || c < 0 || N <= c;
}

bool is_out(const Point& p) {
  return is_out(p.r, p.c);
}

struct Field {
  int Z;
  int grid[MAX_N][MAX_N];
  int score;
  std::shared_ptr<Field> parent;
  std::string prev_command;

  int* operator [](const int index) {
    return grid[index];
  }
  const int* operator [](const int index) const {
    return grid[index];
  }

  int& operator()(const int r, const int c) {
    return grid[r][c];
  }
  const int& operator()(const int r, const int c) const {
    return grid[r][c];
  }
  int& at(int r, int c) {
    return grid[r][c];
  }
  const int& at(const int r, const int c) const {
    return grid[r][c];
  }
};

std::shared_ptr<Field> field;

struct Item {
  int color;
  int r, c;
  int score;
};

std::string solve() {
  auto main_queue_compare = [](const Item& l, const Item& r) {
    return l.score < r.score;
  };
  std::priority_queue<Item, std::vector<Item>, decltype(main_queue_compare)> main_queue{main_queue_compare};
  for (int r = 0; r < N; ++r) {
    for (int c = 0; c < N; ++c) {
      if (field->at(r, c)> 0) {
        Item item{field->at(r, c), r, c};
        struct Node {
          Point p;
          int score;
        };
        const auto compare = [](const Node& l, const Node& r) {
          return l.score > r.score;
        };
        std::priority_queue<Node, std::vector<Node>, decltype(compare)> queue{compare};
        Node init{{r, c}, 0};
        queue.push(init);
        bool done[MAX_N][MAX_N];
        int best = N * N;
        while (!queue.empty()) {
          Node node = queue.top();
          queue.pop();
          if (done[node.p.r][node.p.c]) {
            continue;
          }
          done[node.p.r][node.p.c] = true;
          if (field->at(node.p.r, node.p.c) == -1) {
            best = node.score;
            break;
          }
          for (const auto& D : OFS) {
            const Point np{node.p.r + D[1], node.p.c + D[0]};
            if (is_out(np)) {
              continue;
            }
            Node next{np, node.score + 1};
          }
        }
        item.score = -(best * (MAX_C + 1) + field->at(r, c));
        main_queue.push(item);
      }
    }
  }
  std::shared_ptr<Field> best = field;

  for (;;) {
    if (main_queue.empty()) {
      break;
    }
    if (timer.TLE()) {
      break;
    }

    std::vector<Item> negatives;
    while (!main_queue.empty()) {
      Item item = main_queue.top();
      main_queue.pop();

      struct Node {
        Point p;
        int prev;
      };
      std::queue<Node> queue;
      int memo[MAX_N][MAX_N];
      std::memset(memo, -1, sizeof(memo));
      queue.push({{item.r, item.c}, OFS.size()});
      Point hole{-1, -1};
      while(!queue.empty()) {
        Node node = queue.front();
        queue.pop();
        if (memo[node.p.r][node.p.c] != -1) {
          continue;
        }
        memo[node.p.r][node.p.c] = node.prev;
        if (best->at(node.p.r, node.p.c) == -1) {
          hole = node.p;
          break;
        }
        for (int index = 0; index < OFS.size(); ++index) {
          Point np = {node.p.r + OFS[index][1], node.p.c + OFS[index][0]};
          if (is_out(np)) {
            continue;
          }
          if (best->at(np.r, np.c) != 0) {
            continue;
          }
          queue.push({np, index});
        }
      }
      if (hole.r == -1) {
        negatives.push_back(item);
        continue;
      }
      std::stringstream ss;
      Point ite = hole;
      while (memo[ite.r][ite.c] != OFS.size()) {
        const int index = memo[ite.r][ite.c];
        ss << COMMANDS[index];
        ite.r -= OFS[index][1];
        ite.c -= OFS[index][0];
      }
      std::shared_ptr<Field> next = std::make_shared<Field>(*best);
      next->parent = best;
      next->at(item.r, item.c) = 0;
      next->prev_command = ss.str();
      next->Z -= next->prev_command.size();
      next->score += std::max(0, next->Z * (item.color - 1));
      best = next;
    }

    for (const auto& item : negatives) {
      main_queue.push(item);
    }
  }

  std::cerr << "Score: " << best->score << std::endl;
  std::stringstream ss;
  std::shared_ptr<Field> ite = best;
  while (ite->parent) {
    ss << ite->prev_command;
    ite = ite->parent;
  }
  std::string ans = ss.str();
  std::reverse(ans.begin(), ans.end());

  return ans;
}

int main() {
  field = std::make_shared<Field>();
  std::cin >> N;
  std::cin >> C;
  std::cin >> H;

  // read grid
  for (int k = 0; k < N * N; k++) {
    int r = k / N;
    int c = k % N;
    std::cin >> field->at(r, c);
  }

  field->Z = N * N;

  const auto ans = solve();

  std::cout << ans << std::endl;

  std::cout.flush();
}
