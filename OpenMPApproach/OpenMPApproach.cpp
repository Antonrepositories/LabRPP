#include <iostream>
#include <omp.h>
#include <chrono>
using namespace std;

double f(double x) {
    return x * x; 
}

double trapezoidal_integration(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));

#pragma omp parallel for reduction(+:sum)
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
    for (int num_threads : {2, 4, 8}) {
        omp_set_num_threads(num_threads);

        auto start = chrono::high_resolution_clock::now();
        double result = trapezoidal_integration(a, b, n);
        auto end = chrono::high_resolution_clock::now();

        chrono::duration<double> elapsed = end - start;
        cout << "Number of threads: " << num_threads << endl;
        cout << "Result: " << result << endl;
        cout << "Time: " << elapsed.count() << " seconds" << endl;
        cout << "------------------------------" << endl;
    }
    return 0;
}