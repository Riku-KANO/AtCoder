#include<iostream>
#include<string>
using namespace std;

int main(void) {
  string s; cin >> s;
  int remain;
  if(s == "Monday") remain=5;
  else if(s == "Tuesday") remain=4;
  else if(s == "Wednesday") remain=3;
  else if(s == "Thursday") remain=2;
  else if(s == "Friday") remain=1;
  cout << remain << endl;
  return 0;
}