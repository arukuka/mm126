// C++11
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

constexpr int MAX_N = 10;

int C, H;

struct Field {
  int N;
  int grid[MAX_N][MAX_N];

  int* operator [](const int index) {
    return grid[index];
  }

  int& operator()(const int r, const int c) {
    return grid[r][c];
  }

  bool is_out(const int r, const int c) {
    return r < 0 || N <= r
           || c < 0 || N <= c;
  }
};

Field field;

std::string solve() {
  std::string ans;
  return ans;
}

int main() {
  std::cin >> field.N;
  std::cin >> C;
  std::cin >> H;

  // read grid
  for (int k = 0; k < field.N * field.N; k++) {
    int r = k / field.N;
    int c = k % field.N;
    std::cin >> field(r, c);
  }

  const auto ans = solve();

  std::cout << ans << std::endl;

  std::cout.flush();
}
