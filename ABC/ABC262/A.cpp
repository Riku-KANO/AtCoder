#include<iostream>
using namespace std;

int main() {
    int y; cin >> y;
    while(1) {
        if(y % 4 == 2) {
            cout << y << endl;
            return 0;
        }
        y++;
    }
    return 0;
}