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
constexpr int MAX_N = 30;
constexpr int MAX_C = 9;
constexpr std::array<std::array<int, 2>, 4> OFS{{
  {1, 0},
  {0, 1},
  {-1, 0},
  {0, -1}
}};
constexpr char DIR_COMMANDS[5] = "RDLU";

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

struct Command {
  Point point;
  char type;
  char dir;
};

struct Field {
  int Z;
  int grid[MAX_N][MAX_N];
  int score;
  std::shared_ptr<Field> parent;
  std::vector<Command> prev_command;

  Field(const std::shared_ptr<Field> parent)
      : Z(parent->Z), score(parent->score)
      , parent(parent), prev_command() {
    std::memcpy(grid, parent->grid, sizeof(grid));
  }
  Field(const Field& parent)
      : Z(parent.Z), score(parent.score)
      , prev_command() {
    std::memcpy(grid, parent.grid, sizeof(grid));
  }
  Field() {}

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

std::ostream& operator<<(std::ostream& os, const Point& p)
{
  os << "(" << p.r << ", " << p.c << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, const Item& item)
{
  os << "{point: " << Point{item.r, item.c} << ", color: " << item.color << ", score" << item.score << "}";
  return os;
}

std::ostream& operator<<(std::ostream& os, const Command& cmd)
{
  os << cmd.point.r << " " << cmd.point.c << " " << cmd.type << " " << cmd.dir;
  return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
  os << "[";
  for (size_t i = 0; i < vec.size(); ++i) {
    if (i) os << ", ";
    os << vec[i];
  }
  os << "]";
  return os;
}

std::vector<Command> solve() {
  auto main_queue_compare = [](const Item& l, const Item& r) {
    return l.score < r.score;
  };
  std::priority_queue<Item, std::vector<Item>, decltype(main_queue_compare)> main_queue{main_queue_compare};
  for (int r = 0; r < N; ++r) {
    for (int c = 0; c < N; ++c) {
      if (field->at(r, c) > 1) {
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
        bool done[MAX_N][MAX_N] = {0};
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
            if (field->at(np.r, np.c) == 1) {
              continue;
            }
            Node next{np, node.score + 1};
            queue.push(next);
          }
        }
        item.score = -(best * (MAX_C + 1) + MAX_C - field->at(r, c));
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
        Point point;
        Command command;
      };
      std::queue<Node> queue;
      Command memo[MAX_N][MAX_N];
      std::memset(memo, -1, sizeof(memo));
      queue.push({{item.r, item.c}, {-1, -1, '\0', '\0'}});
      Point hole{-1, -1};
      while(!queue.empty()) {
        Node node = queue.front();
        queue.pop();
        if (memo[node.point.r][node.point.c].type != -1) {
          continue;
        }
        memo[node.point.r][node.point.c] = node.command;
        if (best->at(node.point.r, node.point.c) == -1) {
          hole = node.point;
          break;
        }
        for (int index = 0; index < OFS.size(); ++index) {
          Point np = {node.point.r + OFS[index][1], node.point.c + OFS[index][0]};
          if (is_out(np)) {
            continue;
          }
          if (best->at(np.r, np.c) > 0) {
            continue;
          }
          queue.push({np, {node.point, 'M', DIR_COMMANDS[index]}});
          for (;;) {
            np.r += OFS[index][1];
            np.c += OFS[index][0];
            if (is_out(np) || best->at(np.r, np.c) > 0) {
              np.r -= OFS[index][1];
              np.c -= OFS[index][0];
              break;
            }
            if (best->at(np.r, np.c) == -1) {
              break;
            }
          }
          queue.push({np, {node.point, 'S', DIR_COMMANDS[index]}});
        }
      }
      if (hole.r == -1) {
        negatives.push_back(item);
        continue;
      }
      std::vector<Command> command;
      Point ite = hole;
      while (memo[ite.r][ite.c].type != '\0') {
        const Command cmd = memo[ite.r][ite.c];
        ite = cmd.point;
        command.push_back(cmd);
      }
      std::shared_ptr<Field> next = std::make_shared<Field>(best);
      next->parent = best;
      next->at(item.r, item.c) = 0;
      next->prev_command = command;
      next->Z -= next->prev_command.size();
      next->score += std::max(0, (next->Z + 1) * (item.color - 1));
      best = next;
    }

    for (const auto& item : negatives) {
      main_queue.push(item);
    }
  }

  std::cerr << "Score: " << best->score << std::endl;
  std::vector<Command> ans;
  std::shared_ptr<Field> ite = best;
  while (ite->parent) {
    ans.insert(ans.end(), ite->prev_command.begin(), ite->prev_command.end());
    ite = ite->parent;
  }
  std::reverse(ans.begin(), ans.end());

  if (ans.size() > N * N) {
    ans.resize(N * N);
  }

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

  std::cout << ans.size() << std::endl;
  for (const auto& c : ans) {
    std::cout << c << std::endl;
  }

  std::cout.flush();
}
