// C++11
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <memory>
#include <chrono>
#include <cstring>
#include <random>
#include <unordered_map>
#include <cassert>

template<typename T>
void debug_print(const std::string file, const int line, const std::string func, std::string name, const T& t)
{
  std::cerr << file << ":" << line << " (" << func << ") " << name << ": " << t << std::endl;
}

#define DBG(x) debug_print(__FILE__, __LINE__, __func__, #x, x)

struct Point {
  std::int8_t r, c;

  bool operator==(const Point& p) const {
    return r == p.r && c == p.c;
  }
};

constexpr double TLE = 9.0;
constexpr int MAX_N = 30;
constexpr int MAX_C = 9;
constexpr std::array<std::array<std::int8_t, 2>, 4> OFS{{
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
std::mt19937 random_engine(std::random_device{}());

int N, C, H;

bool is_out(const int r, const int c) {
  return r < 0 || N <= r
         || c < 0 || N <= c;
}

bool is_out(const Point& p) {
  return is_out(p.r, p.c);
}

int manhattan_distance(int x1, int y1, int x2, int y2) {
  return std::abs(x2 - x1) + std::abs(y2 - y1);
}

struct Command {
  Point point;
  char type;
  char dir;
};

template <int... args>
struct _Field;

template<>
struct _Field<0> {
  using this_type = _Field<0>;
  using cell_type = std::int8_t;
  static constexpr cell_type HOLE = -1;

  int Z;
  cell_type grid[MAX_N][MAX_N];
  int score;
  std::shared_ptr<this_type> parent;
  std::vector<Command> prev_command;
  int prev_target_color;

  _Field<0>(const std::shared_ptr<this_type> parent)
      : Z(parent->Z), score(parent->score)
      , parent(parent), prev_command() {
    std::memcpy(grid, parent->grid, sizeof(grid));
  }
  _Field<0>(const this_type& parent)
      : Z(parent.Z), score(parent.score)
      , prev_command() {
    std::memcpy(grid, parent.grid, sizeof(grid));
  }
  _Field<0>() {}

private:
  cell_type& at(int r, int c) {
#ifndef NDEBUG
    assert(!is_out(r, c));
#endif
    return grid[r][c];
  }
  const cell_type& at(const int r, const int c) const {
#ifndef NDEBUG
    assert(!is_out(r, c));
#endif
    return grid[r][c];
  }
public:

  bool is_block(const int r, const int c) const {
    return grid[r][c] > 0;
  }

  bool is_hole(const int r, const int c) const {
    return grid[r][c] == HOLE;
  }

  int get_block_color(const int r, const int c) const {
    assert(is_block(r, c));
    return grid[r][c];
  }

  void set(const int r, const int c, const cell_type v) {
    grid[r][c] = v;
  }

  void clear(const int r, const int c) {
    grid[r][c] = 0;
  }

  static std::shared_ptr<this_type> make_child(std::shared_ptr<this_type> parent, std::vector<Command>& command) {
    const auto point = command.back().point;
    const auto color = parent->at(point.r, point.c);
    std::shared_ptr<this_type> next = std::make_shared<this_type>(parent);
    next->parent = parent;
    next->at(point.r, point.c) = 0;
    next->prev_command = command;
    next->Z -= next->prev_command.size();
    next->score += std::max(0, (next->Z + 1) * (color - 1));
    next->prev_target_color = color;
    return next;
  }
};

template<int GRID_SIZE>
struct _Field<1, GRID_SIZE> {
  using this_type = _Field<1, GRID_SIZE>;
  using cell_type = std::uint8_t;
  using cell_vec_type = std::uint64_t;
  static constexpr int CELL_BITS = 4;
  static constexpr int CELL_VEC_BITS = sizeof(cell_vec_type) * 8;
  static constexpr cell_type MASK = (1 << CELL_BITS) - 1;
  static constexpr cell_type HOLE = (1 << CELL_BITS) - 1;
  static_assert((static_cast<cell_type>(-1) & MASK) == HOLE, "");
  static constexpr int CELL_VEC_NUMS
      = (GRID_SIZE * GRID_SIZE + (GRID_SIZE - 1)) // maximum access
        * CELL_BITS / CELL_VEC_BITS
        + 1;

  int Z;
  cell_vec_type grid[CELL_VEC_NUMS];
  int score;
  std::shared_ptr<this_type> parent;
  std::vector<Command> prev_command;
  cell_type prev_target_color;

  _Field<1, GRID_SIZE>(const std::shared_ptr<this_type> parent)
      : Z(parent->Z), score(parent->score)
      , parent(parent), prev_command() {
    std::memcpy(grid, parent->grid, sizeof(grid));
  }
  _Field<1, GRID_SIZE>(const this_type& parent)
      : Z(parent.Z), score(parent.score)
      , prev_command() {
    std::memcpy(grid, parent.grid, sizeof(grid));
  }
  _Field<1, GRID_SIZE>() {}

private:
  struct CellAccessInfo {
    std::uint16_t index;
    std::uint16_t shift;
    CellAccessInfo(const int r, const int c) {
      const int bits = (r * GRID_SIZE + c) * CELL_BITS;
      this->index = bits / CELL_VEC_BITS;
      this->shift = bits % CELL_VEC_BITS;
      assert(this->index < CELL_VEC_NUMS);
    }
  };
  const cell_type at(const int r, const int c) const {
#ifndef NDEBUG
    assert(!is_out(r, c));
#endif
    const CellAccessInfo info{r, c};
    return static_cast<cell_type>((grid[info.index] >> info.shift) & MASK);
  }
public:

  bool is_block(const int r, const int c) const {
    return !is_hole(r, c) && at(r, c) > 0;
  }

  bool is_hole(const int r, const int c) const {
    return at(r, c) == HOLE;
  }

  cell_type get_block_color(const int r, const int c) const {
    assert(!is_hole(r, c));
    return at(r, c);
  }

  void set(const int r, const int c, const cell_type v) {
    clear(r, c);
    const CellAccessInfo info{r, c};
    grid[info.index] |= static_cast<cell_vec_type>(v & MASK) << info.shift;
  }

  void clear(const int r, const int c) {
    const CellAccessInfo info{r, c};
    grid[info.index] &= ~(static_cast<cell_vec_type>(MASK) << info.shift);
  }

  static std::shared_ptr<this_type> make_child(std::shared_ptr<this_type> parent, std::vector<Command>& command) {
    const auto point = command.back().point;
    const auto color = parent->at(point.r, point.c);
    std::shared_ptr<this_type> next = std::make_shared<this_type>(parent);
    next->parent = parent;
    next->clear(point.r, point.c);
    next->prev_command = command;
    next->Z -= next->prev_command.size();
    next->score += std::max(0, (next->Z + 1) * (color - 1));
    next->prev_target_color = color;
    return next;
  }
};

using Field = _Field<0>;

std::shared_ptr<Field> field;

std::ostream& operator<<(std::ostream& os, const Point& p)
{
  os << "(" << static_cast<int>(p.r) << ", " << static_cast<int>(p.c) << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, const Command& cmd)
{
  os << static_cast<int>(cmd.point.r) << " " << static_cast<int>(cmd.point.c) << " " << cmd.type << " " << cmd.dir;
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

std::shared_ptr<Field> bfs(const std::shared_ptr<Field> src, const Point& target_point)
{
  struct Node {
    Point point;
    Command command;
  };
  std::queue<Node> queue;
  Command memo[MAX_N][MAX_N];
  std::memset(memo, -1, sizeof(memo));
  queue.push({target_point, {-1, -1, '\0', '\0'}});
  Point hole{-1, -1};
  while(!queue.empty()) {
    Node node = queue.front();
    queue.pop();
    if (memo[node.point.r][node.point.c].type != -1) {
      continue;
    }
    memo[node.point.r][node.point.c] = node.command;
    if (src->is_hole(node.point.r, node.point.c)) {
      hole = node.point;
      break;
    }
    for (int index = 0; index < OFS.size(); ++index) {
      Point np = {static_cast<int8_t>(node.point.r + OFS[index][1]), static_cast<int8_t>(node.point.c + OFS[index][0])};
      if (is_out(np)) {
        continue;
      }
      if (src->is_block(np.r, np.c)) {
        continue;
      }
      queue.push({np, {node.point, 'M', DIR_COMMANDS[index]}});
      for (;;) {
        np.r += OFS[index][1];
        np.c += OFS[index][0];
        if (is_out(np) || src->is_block(np.r, np.c)) {
          np.r -= OFS[index][1];
          np.c -= OFS[index][0];
          break;
        }
        if (src->is_hole(np.r, np.c)) {
          break;
        }
      }
      queue.push({np, {node.point, 'S', DIR_COMMANDS[index]}});
    }
  }
  if (hole.r == -1) {
    return src;
  }
  std::vector<Command> command;
  Point ite = hole;
  while (memo[ite.r][ite.c].type != '\0') {
    const Command cmd = memo[ite.r][ite.c];
    ite = cmd.point;
    command.push_back(cmd);
  }
  std::shared_ptr<Field> next = Field::make_child(src, command);

  return next;
}


std::shared_ptr<Field> pseudo_dijkstra(const std::shared_ptr<Field> src, int filter_color = 0)
{
  struct Item {
    Point point;
    int color;
    std::size_t index;

    bool operator==(const Item& item) const {
      return index == item.index && point == item.point;
    }
    bool operator!=(const Item& item) const {
      return !(this->operator==(item));
    }
  };
  std::vector<Item> items;
  for (std::int8_t r = 0; r < N; ++r) {
    for (std::int8_t c = 0; c < N; ++c) {
      if (!src->is_block(r, c)) {
        continue;
      }
      if (src->get_block_color(r, c) <= filter_color) {
        continue;
      }
      items.emplace_back(Item{Point{r, c}, src->get_block_color(r, c), items.size()});
    }
  }
  if (items.size() == 0) {
    return src;
  }

  struct Node {
    Item item;
    Command command;
    int score;
  };
  struct ItemHasher {
    std::size_t operator()(const Item& item) const {
      return item.index * N * N + (item.point.r * N + item.point.c);
    }
  };
  auto node_comparator = [](const Node& lhs, const Node& rhs) {
    return lhs.score < rhs.score;
  };
  std::priority_queue<Node, std::vector<Node>, decltype(node_comparator)> queue{node_comparator};
  std::unordered_map<Item, Command, ItemHasher> memo;
  for (const auto& item : items) {
    queue.emplace(Node{item, Command{-1, -1, '\0', '\0'}, src->Z * item.color});
  }
  Item hole{{-1, -1}, -1, items.size()};
  while(!queue.empty()) {
    Node node = queue.top();
    queue.pop();
    if (memo.count(node.item)) {
      continue;
    }
    memo[node.item] = node.command;
    if (src->is_hole(node.item.point.r, node.item.point.c)) {
      hole = node.item;
      break;
    }

    for (int index = 0; index < OFS.size(); ++index) {
      Point np = {static_cast<std::int8_t>(node.item.point.r + OFS[index][1]), static_cast<std::int8_t>(node.item.point.c + OFS[index][0])};
      if (is_out(np)) {
        continue;
      }
      if (src->is_block(np.r, np.c)) {
        continue;
      }
      queue.emplace(Node{
        Item{np, node.item.color, node.item.index},
        Command{node.item.point, 'M', DIR_COMMANDS[index]},
        node.score - node.item.color
      });
      for (;;) {
        np.r += OFS[index][1];
        np.c += OFS[index][0];
        if (is_out(np) || src->is_block(np.r, np.c)) {
          np.r -= OFS[index][1];
          np.c -= OFS[index][0];
          break;
        }
        if (src->is_hole(np.r, np.c)) {
          break;
        }
      }
      queue.emplace(Node{
        Item{np, node.item.color, node.item.index},
        Command{node.item.point, 'S', DIR_COMMANDS[index]},
        node.score - node.item.color
      });
    }
  }
  if (hole.color == -1) {
    return src;
  }

  const Item item = items[hole.index];
  std::vector<Command> command;
  Item ite = hole;
  while (memo[ite].type != '\0') {
    const Command cmd = memo[ite];
    ite.point = cmd.point;
    command.push_back(cmd);
  }
  std::shared_ptr<Field> next = Field::make_child(src, command);

  return next;
}

std::shared_ptr<Field> solve_greedy(const std::shared_ptr<Field> src) {
  struct Item {
    int color;
    Point point;
    int score;
  };

  std::vector<Point> holes;
  for (std::int8_t r = 0; r < N; ++r) {
    for (std::int8_t c = 0; c < N; ++c) {
      if (src->is_hole(r, c)) {
        holes.push_back({r, c});
      }
    }
  }

  auto main_queue_compare = [](const Item& l, const Item& r) {
    return l.score < r.score;
  };
  std::priority_queue<Item, std::vector<Item>, decltype(main_queue_compare)> main_queue{main_queue_compare};
  for (std::int8_t r = 0; r < N; ++r) {
    for (std::int8_t c = 0; c < N; ++c) {
      if (src->is_block(r, c)) {
        Item item{src->get_block_color(r, c), {r, c}};
        int dist = std::numeric_limits<int>::max();
        for (const auto& hole : holes) {
          dist = std::min(dist, manhattan_distance(c, r, hole.c, hole.r));
        }
        item.score = -(dist * (MAX_C + 1) + MAX_C - src->get_block_color(r, c));
        main_queue.push(item);
      }
    }
  }

  std::shared_ptr<Field> best = src;

  while (!main_queue.empty()) {
    Item item = main_queue.top();
    main_queue.pop();

    best = bfs(best, item.point);
  }
  return best;
}

bool check(std::shared_ptr<Field> src, std::vector<Command>& commands, bool deep = false)
{
  Point point = commands.back().point;
  if (!src->is_block(point.r, point.c)) {
    return false;
  }
  if (!deep) {
    return true;
  }
  int rev_type[256];
  for (int i = 0; DIR_COMMANDS[i] != '\0'; ++i) {
    rev_type[DIR_COMMANDS[i]] = i;
  }
  for (auto command_ite = commands.rbegin(); command_ite != commands.rend(); command_ite++) {
    const Command& command = *command_ite;
    int dx = OFS[rev_type[command.dir]][0];
    int dy = OFS[rev_type[command.dir]][1];
    point.r += dy;
    point.c += dx;
    if (is_out(point) || src->is_block(point.r, point.c)) {
      return false;
    }
    if (command.type == 'S') {
        for (;;) {
          point.r += dy;
          point.c += dx;
          if (is_out(point) || src->is_block(point.r, point.c)) {
            point.r -= dy;
            point.c -= dx;
            break;
          }
          if (src->is_hole(point.r, point.c)) {
            break;
          }
        }
    }
  }
  return src->is_hole(point.r, point.c);
}

std::shared_ptr<Field> use_history(
    const std::shared_ptr<Field> src,
    const std::vector<std::shared_ptr<Field>>& histories) {
  std::shared_ptr<Field> applied = src;
  for (const std::shared_ptr<Field> history : histories) {
    if (check(applied, history->prev_command)) {
      if (!check(applied, history->prev_command, true)) {
        applied = bfs(applied, history->prev_command.back().point);
      } else {
        applied = Field::make_child(applied, history->prev_command);
      }
    }
  }
  return applied;
}

std::shared_ptr<Field> neighbor_insert(const std::shared_ptr<Field> src)
{
  struct ScoredField {
    std::shared_ptr<Field> field;
    int score;
  };
  std::vector<ScoredField> nodes;
  {
    std::vector<ScoredField> all_state;
    std::shared_ptr<Field> ite = src;
    while (ite) {
      ScoredField node;
      node.field = ite;
      node.score = std::max(0, ite->Z + static_cast<int>(ite->prev_command.size()) - 1) * (MAX_C + 1 - ite->prev_target_color);
      all_state.push_back(node);
      ite = ite->parent;
    }
    if (all_state.size() <= 2) {
      return src;
    }
    std::reverse(all_state.begin(), all_state.end());
    nodes.insert(nodes.end(), all_state.begin() + 1, all_state.end() - 1);
  }
  std::vector<int> score_roulette;
  int score_acc = 0;
  for (const auto& node : nodes) {
    score_acc += node.score;
    score_roulette.emplace_back(score_acc);
  }
  const auto target_ite = std::lower_bound(
    score_roulette.begin(),
    score_roulette.end(),
    std::uniform_int_distribution<>{0, score_acc - 1}(random_engine)
  );
  if (target_ite == score_roulette.end()) {
    DBG("OMG!");
    return src;
  }
  const std::size_t target_index = std::distance(score_roulette.begin(), target_ite);
  const std::shared_ptr<Field> target_child = nodes[target_index].field;
  const std::shared_ptr<Field> target = target_child->parent;
  const int target_child_score = target_child->score - target->score;

  std::shared_ptr<Field> next = pseudo_dijkstra(target, target_child->prev_target_color);

  std::vector<std::shared_ptr<Field>> histories;
  for (int i = target_index; i < nodes.size(); ++i) {
    histories.emplace_back(nodes[i].field);
  }
  next = use_history(next, histories);

  next = solve_greedy(next);
  return next;
}

std::vector<Command> solve() {
  std::shared_ptr<Field> best = solve_greedy(field);

  int iteration;
  for (int iteration = 0;; ++iteration) {
    if (timer.TLE()) {
      DBG(iteration);
      break;
    }

    std::shared_ptr<Field> next = neighbor_insert(best);

    if (next->score > best->score) {
      best = next;
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

void cal_max_score() {
  std::vector<int> blocks;
  for (int r = 0; r < N; ++r) {
    for (int c = 0; c < N; ++c) {
      if (!field->is_block(r, c)) {
        continue;
      }
      blocks.push_back(field->get_block_color(r, c));
    }
  }
  std::sort(blocks.begin(), blocks.end(), std::greater<int>{});
  int score = 0;
  for (int i = 0; i < blocks.size(); ++i) {
    score += (field->Z - i) * (blocks[i] - 1);
  }
  std::cerr << score << std::endl;
}

int main(const int argc, const char** argv) {
  field = std::make_shared<Field>();
  std::cin >> N;
  std::cin >> C;
  std::cin >> H;

  // read grid
  for (int k = 0; k < N * N; k++) {
    int r = k / N;
    int c = k % N;
    int v;
    std::cin >> v;
    field->set(r, c, v);
  }

  field->Z = N * N;

  if (argc >= 2 && std::string(argv[1]) == "cal_max_score") {
    cal_max_score();
    return 0;
  }

  const auto ans = solve();

  std::cout << ans.size() << std::endl;
  for (const auto& c : ans) {
    std::cout << c << std::endl;
  }

  std::cout.flush();
}
