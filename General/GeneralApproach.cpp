#include <iostream>
#include <chrono>
using namespace std;
double f(double x) {
    return x * x;
}

double trapezoidal_integration(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));

    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += f(x);
    }

    return sum * h;
}

int main() {
    double a = 0.0;
    double b = 10.0;
    int n = 900000000;

    auto start = chrono::high_resolution_clock::now();
    double result = trapezoidal_integration(a, b, n);
    auto end = chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Result: " << result << endl;
    std::cout << "Time: " << elapsed.count() << " seconds" << endl;

    return 0;
}