#include<iostream>
#include<cmath>
#include<limits>
#include<iomanip>
using namespace std;

int main() {
    cout << fixed << setprecision(10);
    int a, b, d; cin >> a >> b >> d;
    double r = sqrt(a*a+b*b);
    double theta = atan2(b, a);
    double phi = theta + d / 180. * acos(-1);
    cout << r * cos(phi) << " " << r * sin(phi) << endl;
    return 0;
}