// C++11
#include <algorithm>
#include <cmath>
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

// https://github.com/kmyk/atcoder-heuristic-contest-002/blob/074128b66b88dc9083f2bd4650c9b86a71b428cb/main.cpp#L9-L30
class xor_shift_128 {
public:
    typedef uint32_t result_type;
    xor_shift_128(uint32_t seed = 42) {
        set_seed(seed);
    }
    void set_seed(uint32_t seed) {
        a = seed = 1812433253u * (seed ^ (seed >> 30));
        b = seed = 1812433253u * (seed ^ (seed >> 30)) + 1;
        c = seed = 1812433253u * (seed ^ (seed >> 30)) + 2;
        d = seed = 1812433253u * (seed ^ (seed >> 30)) + 3;
    }
    uint32_t operator() () {
        uint32_t t = (a ^ (a << 11));
        a = b; b = c; c = d;
        return d = (d ^ (d >> 19)) ^ (t ^ (t >> 8));
    }
    static constexpr uint32_t max() { return std::numeric_limits<result_type>::max(); }
    static constexpr uint32_t min() { return std::numeric_limits<result_type>::min(); }
private:
    uint32_t a, b, c, d;
};

struct Point {
  std::int16_t r, c;

  bool operator==(const Point& p) const {
    return r == p.r && c == p.c;
  }
  bool operator!=(const Point& p) const {
    return !this->operator==(p);
  }
  std::int16_t val() const {
    return (static_cast<std::int16_t>(r) << 5) | c;
  }
  bool operator<(const Point& p) const {
    return val() < p.val();
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
  double secs = 0;
  Timer() : start(std::chrono::high_resolution_clock::now()) {}
  bool TLE() {
    const auto now = std::chrono::high_resolution_clock::now();
    const auto secs = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count() / 1e3;
    this->secs = secs;
    return ::TLE < secs;
  }
};

Timer timer;
xor_shift_128 random_engine(20210505);

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
  using cell_type = std::int16_t;
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
  static std::shared_ptr<this_type> make_child(std::shared_ptr<this_type> parent, std::vector<Command>& command, const Point& next_point) {
    const auto point = command.back().point;
    const auto color = parent->at(point.r, point.c);
    std::shared_ptr<this_type> next = std::make_shared<this_type>(parent);
    next->parent = parent;
    next->at(point.r, point.c) = 0;
    next->prev_command = command;
    next->Z -= next->prev_command.size();
    next->prev_target_color = color;
    if (next_point.r != -1) {
      next->set(next_point.r, next_point.c, color);
    } else {
      next->score += std::max(0, (next->Z + 1) * (color - 1));
    }
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

std::shared_ptr<Field> bfs(const std::shared_ptr<Field> src, const Point& from, const Point& to)
{
  if (from == to) {
    return src;
  }
  struct Node {
    Point point;
    Command command;
  };
  std::queue<Node> queue;
  Command memo[MAX_N][MAX_N];
  std::memset(memo, -1, sizeof(memo));
  queue.push({from, {-1, -1, '\0', '\0'}});
  while(!queue.empty()) {
    Node node = queue.front();
    queue.pop();
    if (memo[node.point.r][node.point.c].type != -1) {
      continue;
    }
    memo[node.point.r][node.point.c] = node.command;
    if (src->is_hole(node.point.r, node.point.c)) {
      continue;
    }
    if (node.point == to) {
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
  if (memo[to.r][to.c].type == -1) {
    return src;
  }
  std::vector<Command> command;
  Point ite = to;
  while (memo[ite.r][ite.c].type != '\0') {
    const Command cmd = memo[ite.r][ite.c];
    ite = cmd.point;
    command.push_back(cmd);
  }
  std::shared_ptr<Field> next = Field::make_child(src, command, to);

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

std::int8_t normalize(std::int8_t x)
{
  if (x == 0) {
    return 0;
  }
  return x / std::abs(x);
}

bool in_range(const Point& from, const Point& to, const Point& point) {
  Point min{std::min(from.r, to.r), std::min(from.c, to.c)};
  Point max{std::max(from.r, to.r), std::max(from.c, to.c)};
  return min.r <= point.r && point.r <= max.r
         && min.c <= point.c && point.c <= max.c;
}

std::pair<Point, int> search_receiver(const std::shared_ptr<Field> src, const Point& receiver, const std::set<Point>& filter) {
  if (src->is_block(receiver.r, receiver.c)) {
    return {receiver, 0};
  }

  struct Info {
    int count;
    Point from;
  };
  struct Node {
    Point point;
    Info info;
  };

  std::vector<Node> targets;
  for (const auto& d1 : OFS) {
    Point ite = receiver;
    ite.r += d1[1];
    ite.c += d1[0];
    if (is_out(ite.r, ite.c)) {
      continue;
    }
    if (src->is_block(ite.r, ite.c)) {
      targets.emplace_back(Node{ite, Info{0, ite}});
      continue;
    }
    for (;;) {
      if (is_out(ite.r, ite.c)) {
        break;
      }
      if (src->is_block(ite.r, ite.c)) {
        targets.emplace_back(Node{ite, Info{0, ite}});
        break;
      }
      for (const auto& d2 : OFS) {
        if (d2 == d1 || (d2[0] * -1 == d1[0] && d2[1] * -1 == d1[1])) {
          continue;
        }
        Point ite2{ite};
        ite2.r += d2[1];
        ite2.c += d2[0];
        if (is_out(ite2.r, ite2.c)) {
          continue;
        }
        if (src->is_block(ite2.r, ite2.c)) {
          targets.emplace_back(Node{ite2, Info{0, ite2}});
          continue;
        }
        for (;;) {
          if (is_out(ite2.r, ite2.c)) {
            break;
          }
          if (src->is_block(ite2.r, ite2.c)) {
            targets.emplace_back(Node{ite2, Info{0, ite2}});
            break;
          }
          ite2.r += d2[1];
          ite2.c += d2[0];
        }
      }
      ite.r += d1[1];
      ite.c += d1[0];
    }
  }

  std::queue<Node> queue;
  for (const auto& target : targets) {
    if (filter.count(target.point)) {
      continue;
    }
    queue.push(target);
  }

  Info memo[MAX_N][MAX_N];
  std::memset(memo, -1, sizeof(memo));
  while(!queue.empty()) {
    Node node = queue.front();
    queue.pop();
    if (src->is_hole(node.point.r, node.point.c)) {
      continue;
    }
    if (memo[node.point.r][node.point.c].count != -1) {
      continue;
    }
    memo[node.point.r][node.point.c] = node.info;
    if (node.point == receiver) {
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
      queue.push({np, Info{node.info.count + 1, node.info.from}});
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
      queue.push({np, Info{node.info.count + 1, node.info.from}});
    }
  }
  Info info = memo[receiver.r][receiver.c];
  return {info.from, info.count};
}

std::pair<Point, int> search_receiver(const std::shared_ptr<Field> src, const Point& receiver, const Point& filter) {
  std::set<Point> set;
  set.insert(filter);
  return search_receiver(src, receiver, set);
}

struct Directions {
  int Z;
  int score;
  Point start;
  struct Info {
    Point dir;
    Point stop_point;
    std::pair<Point, int> target_receiver;
  };
  std::vector<Info> infos;

  bool operator<(const Directions& d) const {
    if (Z != d.Z) {
      return Z < d.Z;
    }
    return score < d.score;
  }
};

constexpr std::int8_t NOT_USED = 32;
constexpr std::int8_t ALREADY_DONE = -1;

int cal_move_times(const Directions& dirs, const std::vector<Directions::Info>::const_iterator end) {
  int sum = 0;
  Point point = dirs.start;
  for (auto ite = dirs.infos.cbegin(); ite != end; ++ite) {
    if (ite->target_receiver.first.r != NOT_USED) {
      ++sum;
    } else {
      sum += manhattan_distance(
          point.c, point.r,
          ite->stop_point.c, ite->stop_point.r
      );
    }
    point = ite->stop_point;
  }
  return sum;
}


Directions make_move(
    const std::shared_ptr<Field> src,
    const Point& to,
    const Point& now,
    const Point& move,
    const bool use_receiver,
    const Directions& directions)
{
  Directions next{directions};

  const Point stop_point{
    static_cast<std::int8_t>(now.r + move.r),
    static_cast<std::int8_t>(now.c + move.c)
  };
  const Point normalized_diff{
    normalize(move.r),
    normalize(move.c)
  };
  const Point receiver{
    static_cast<std::int8_t>(now.r - normalized_diff.r),
    static_cast<std::int8_t>(now.c - normalized_diff.c)
  };

  std::set<Point> filter;
  filter.insert(to);
  for (const auto& info : directions.infos) {
    filter.insert(info.target_receiver.first);
  }

  std::pair<Point, int> target_receiver{{NOT_USED, -1}, -1};
  bool done = false;
  if (is_out(receiver.r, receiver.c)
      || src->is_hole(now.r, now.c)) {
    done = true;
    target_receiver.first.r = ALREADY_DONE;
  }
  if (use_receiver && !done) {
    target_receiver = search_receiver(src, receiver, filter);
  }

  next.infos.emplace_back(Directions::Info{normalized_diff, stop_point, target_receiver});

  return next;
}

void calculate_score(
    const std::shared_ptr<Field> src,
    Directions& dirs)
{
  std::map<Point, std::pair<int, int>> dic;
  std::vector<bool> dones(dirs.infos.size());

  for (std::size_t i = 0; i < dirs.infos.size(); ++i) {
    const auto& info = dirs.infos[i];
    if (info.target_receiver.first.r == ALREADY_DONE) {
      dones[i] = true;
    } else if (info.target_receiver.first.r != NOT_USED) {
      dic[info.target_receiver.first] = {i, info.target_receiver.second};
    }
  }

  Point now = dirs.start;
  std::size_t depth = 0;
  for (auto ite_info = dirs.infos.cbegin(); ite_info != dirs.infos.cend(); ++ite_info) {
    if (!dones[depth]
        && ite_info->target_receiver.first.r != NOT_USED
        && !in_range(ite_info->stop_point, now, ite_info->target_receiver.first)) {
      dirs.Z -= ite_info->target_receiver.second;
      dones[depth] = true;
    }
    Point ite{now.r, now.c};
    auto run = [&](){
      if (dic.count(ite) && !dones[dic[ite].first]) {
        dirs.Z -= dic[ite].second;
        dones[dic[ite].first] = true;
      } else if (src->is_block(ite.r, ite.c)) {
        const int dist = manhattan_distance(
          now.c, now.r,
          ite.c, ite.r
        );
        if (dones[depth]) {
          dirs.Z -= std::min(1, dist);
        } else {
          dirs.Z -= dist;
        }
        dirs.Z -= cal_move_times(dirs, ite_info);
        dirs.score += (src->get_block_color(ite.r, ite.c) - 1) * std::max(0, dirs.Z + 1);
      }
    };
    while (ite != ite_info->stop_point) {
      run();
      ite.r += ite_info->dir.r;
      ite.c += ite_info->dir.c;
    }
    assert(ite == ite_info->stop_point);
    if (depth + 1 >= dirs.infos.size()) {
      run();
    }
    now = ite_info->stop_point;
    ++depth;
  }
}

std::shared_ptr<Field> apply_directions(
    const std::shared_ptr<Field> src,
    Directions& dirs)
{
  std::map<Point, std::pair<int, int>> dic;
  std::vector<bool> dones(dirs.infos.size());

  for (std::size_t i = 0; i < dirs.infos.size(); ++i) {
    const auto& info = dirs.infos[i];
    if (info.target_receiver.first.r == ALREADY_DONE) {
      dones[i] = true;
    } else if (info.target_receiver.first.r != NOT_USED) {
      dic[info.target_receiver.first] = {i, info.target_receiver.second};
    }
  }

  std::shared_ptr<Field> next = src;
  Point now = dirs.start;
  std::size_t depth = 0;
  for (auto ite_info = dirs.infos.cbegin(); ite_info != dirs.infos.cend(); ++ite_info) {
    const Point receiver{
      static_cast<std::int8_t>(now.r - ite_info->dir.r),
      static_cast<std::int8_t>(now.c - ite_info->dir.c)
    };
    if (!dones[depth]
        && ite_info->target_receiver.first.r != NOT_USED
        && !in_range(ite_info->stop_point, now, ite_info->target_receiver.first)) {
      next = bfs(next, ite_info->target_receiver.first, receiver);
      dones[depth] = true;
    }
    Point ite{now.r, now.c};
    auto run = [&](){
      if (dic.count(ite) && !dones[dic[ite].first]) {
        next = bfs(next, ite, receiver);
        dones[dic[ite].first] = true;
      } else if (src->is_block(ite.r, ite.c)) {
        next = bfs(next, ite);
      }
    };
    while (ite != ite_info->stop_point) {
      run();
      ite.r += ite_info->dir.r;
      ite.c += ite_info->dir.c;
    }
    assert(ite == ite_info->stop_point);
    if (depth + 1 >= dirs.infos.size()) {
      run();
    }
    now = ite_info->stop_point;
    ++depth;
  }

  return next;
}

Directions search_directions(const std::shared_ptr<Field> src, const Point now, const Point goal, const Directions& directions, const int length, const Point& prev_dir) {
  if (length <= 0) {
    Directions ans{directions};
    calculate_score(src, ans);
    return ans;
  }
  const std::int8_t dr = goal.r - now.r;
  const std::int8_t dc = goal.c - now.c;
  Directions ans{std::numeric_limits<int>::min(), 0};
  if (length == 1) {
    if (dr == 0 || dc == 0) {
      Point move{dr, dc};
      for (int i = 0; i < 2; ++i) {
        auto next = make_move(src, goal, now, move, i > 0, directions);
        auto dir = next.infos.back().dir;
        if (prev_dir == dir || (prev_dir.r * -1 == dir.r && prev_dir.c * -1 == dir.c)) {
          break;
        }
        auto res = search_directions(
            src, next.infos.back().stop_point,
            goal, next, length - 1, dir);
        if (ans < res) {
          ans = res;
        }
      }
    } else {
      return {std::numeric_limits<int>::min()};
    }
  } else if (length == 2) {
    if (dr == 0 || dc == 0) {
      return {std::numeric_limits<int>::min()};
    }
    const std::array<Point, 2> moves{{
      {0, dc},
      {dr, 0}
    }};
    for (int i = 0; i < 2; ++i) {
      const auto move = moves[i];
      for (int j = 0; j < 2; ++j) {
        auto next = make_move(src, goal, now, move, j > 0, directions);
        auto dir = next.infos.back().dir;
        if (prev_dir == dir || (prev_dir.r * -1 == dir.r && prev_dir.c * -1 == dir.c)) {
          break;
        }
        auto res = search_directions(
            src, next.infos.back().stop_point,
            goal, next, length - 1, dir);
        if (ans < res) {
          ans = res;
        }
      }
    }
  } else {
    for (int i = 0; i < 2; ++i) {
      for (const auto& ofs : OFS) {
        Point move{ofs[1], ofs[0]};
        for (;;) {
          if (is_out(now.r + move.r, now.c + move.c)) {
            break;
          }
          auto next = make_move(src, goal, now, move, i > 0, directions);
          auto dir = next.infos.back().dir;
          if (prev_dir == dir || (prev_dir.r * -1 == dir.r && prev_dir.c * -1 == dir.c)) {
            break;
          }
          auto res = search_directions(
              src, next.infos.back().stop_point,
              goal, next, length - 1, dir);
          if (ans < res) {
            ans = res;
          }

          move.r += ofs[1];
          move.c += ofs[0];
        }
      }
    }
  }

  return ans;
}

std::shared_ptr<Field> solve_greedy_ver2(const std::shared_ptr<Field> src) {
  std::vector<Point> holes;
  for (std::int8_t r = 0; r < N; ++r) {
    for (std::int8_t c = 0; c < N; ++c) {
      if (src->is_hole(r, c)) {
        holes.push_back({r, c});
      }
    }
  }

  struct LightItem {
    int color;
    Point point;
    int score;
  };
  std::vector<LightItem> raw_targets;
  for (std::int8_t r = 0; r < N; ++r) {
    for (std::int8_t c = 0; c < N; ++c) {
      if (src->is_block(r, c)) {
        int dist = std::numeric_limits<int>::max();
        for (const auto& hole : holes) {
          dist = std::min(dist, manhattan_distance(c, r, hole.c, hole.r));
        }
        int color = src->get_block_color(r, c);
        constexpr int MUL = MAX_N * 2;
        raw_targets.emplace_back(LightItem{color, {r, c}, color * MUL + MUL - 1 - dist});
      }
    }
  }
  if (raw_targets.size() == 0) {
    return src;
  }
  std::sort(raw_targets.begin(), raw_targets.end(), [](const LightItem& l, const LightItem& r) {
    return l.score > r.score;
  });

  Directions best_dirs{std::numeric_limits<int>::min()};
  int max = raw_targets[0].color;
  if (max == 1) {
    return src;
  }
  std::int64_t best_score = std::numeric_limits<std::int64_t>::min();
  for (const auto& target : raw_targets) {
    if (target.color < max) {
      break;
    }
    for (const auto& hole : holes) {
      const Directions atom{src->Z, 0, hole};
      constexpr Point prev_dir{0, 0};
      for (int len = 1; len <= 2; ++len) {
        const auto trial = search_directions(src, hole, target.point, atom, len, prev_dir);
        if (trial.infos.size() == 0) {
          continue;
        }
        std::int64_t next_score = trial.Z;
        next_score *= 8 * NOT_USED * NOT_USED;
        next_score += trial.score;
        {
          int sub = 0;
          const auto dir = trial.infos.back().dir;
          if (len == 1) {
            for (const auto& ofs : OFS) {
              if (dir.r * -1 == ofs[1] && dir.c * -1 == ofs[0]) {
                continue;
              }
              Point ite{target.point};
              ite.r += ofs[1];
              ite.c += ofs[0];
              int Z = trial.Z - 1;
              if (dir.r == ofs[1] && dir.c == ofs[0]) {
                ++Z;
              }
              while (!is_out(ite)) {
                if (src->is_block(ite.r, ite.c)) {
                  sub = std::max(
                      sub,
                      Z * (src->get_block_color(ite.r, ite.c) - 1)
                  );
                  --Z;
                }
                ite.r += ofs[1];
                ite.c += ofs[0];
              }
            }
          } else {
            Point ite{target.point};
            ite.r += dir.r;
            ite.c += dir.c;
            int Z = trial.Z - 1;
            while (!is_out(ite)) {
              if (src->is_block(ite.r, ite.c)) {
                sub = std::max(
                    sub,
                    (Z - 1) * (src->get_block_color(ite.r, ite.c) - 1)
                );
                Z -= 2;
              }
              ite.r += dir.r;
              ite.c += dir.c;
            }
          }
          next_score *= 8 * NOT_USED * NOT_USED;
          next_score += sub;
        }
        next_score *= 256;
        next_score += std::uniform_int_distribution<>{0, 255}(random_engine);
        if (best_score < next_score) {
          best_dirs = trial;
          best_score = next_score;
        }
      }
    }
  }

  if (best_dirs.Z == std::numeric_limits<int>::min()) {
    return src;
  }

  auto next = apply_directions(src, best_dirs);
  assert(best_dirs.Z <= next->Z);

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
  std::shared_ptr<Field> best = field;
  for (;;) {
    if (best->Z <= 0 || timer.TLE()) {
      break;
    }
    auto next = solve_greedy_ver2(best);
    if (best == next) {
      break;
    }
    best = next;
  }

  using field_type = std::shared_ptr<Field>;
  auto field_compare = [](const field_type l, const field_type r) {
    return l->score < r->score;
  };
  using queue_type = std::priority_queue<field_type, std::vector<field_type>, decltype(field_compare)>;
  std::vector<queue_type> targets;
  for (int i = 0; i <= field->Z; ++i) {
    targets.emplace_back(queue_type{field_compare});
  }
  targets[field->Z].push(field);

  for (int iteration = 0;; ++iteration) {
    if (timer.TLE()) {
      DBG(iteration);
      break;
    }

    for (int i = field->Z; i > 0; --i) {
      if (timer.TLE()) {
        break;
      }

      auto& queue = targets[i];
      if (!queue.size()) {
        continue;
      }
      auto node = queue.top();
      queue.pop();

      if (node->score > best->score) {
        best = node;
      }

      field_type next_1 = solve_greedy_ver2(node);
      if (node == next_1) {
        continue;
      }
      targets[std::max(0, next_1->Z)].push(next_1);


      int color = next_1->prev_target_color;
      field_type next_2 = pseudo_dijkstra(node, color - 1);
      if (node != next_2) {
        targets[std::max(0, next_2->Z)].push(next_2);
      }
    }
    {
      auto& queue = targets[0];
      if (queue.size()) {
        auto last = queue.top();
        queue.pop();
        if (last->score > best->score) {
          best = last;
        }
      }
      targets[field->Z].push(field);
    }
  }
  DBG(timer.TLE());
  DBG(timer.secs);

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
