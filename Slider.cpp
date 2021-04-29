// C++11
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <set>
#include <string>

using namespace std;

int main() 
{
  int N;
  int C;
  int H;

  cin >> N;
  cin >> C;
  cin >> H;
  vector<int> grid(N*N);
  //read grid
  int moves = 0;
  for (int k=0; k<N*N; k++) 
  {
    cin >> grid[k];
    if (grid[k]>0) moves++;
  }
  cout << moves << endl;

  for (int r=0; r<N; r++)
    for (int c=0; c<N; c++)
      if (grid[r*N+c]>0)
        cout << r << " " << c << " S L" << endl;

  cout.flush();
}