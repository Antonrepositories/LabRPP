#include <iostream>
#include <mpi.h>
#include <chrono>
#include<string>
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

void run_mpi_test() {
    MPI_Init(NULL, NULL);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double a, b;
    int n;

    if (rank == 0) {
        a = 0.0;
        b = 10.0;
        n = 900000000;
    }

    // Розсилка даних всім процесам
    MPI_Bcast(&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double local_a = a + rank * (b - a) / size;
    double local_b = a + (rank + 1) * (b - a) / size;
    int local_n = n / size;

    auto start = chrono::high_resolution_clock::now();
    double local_result = trapezoidal_integration(local_a, local_b, local_n);
    double global_result = 0.0;

    MPI_Reduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    auto end = chrono::high_resolution_clock::now();

    if (rank == 0) {
        chrono::duration<double> elapsed = end - start;
        cout << "Processes: " << size << endl;
        cout << "Result: " << global_result << endl;
        cout << "Time: " << elapsed.count() << " seconds" << endl;
        cout << "-------------------------" << endl;
    }

    MPI_Finalize();
}


int main() {
    run_mpi_test();
    return 0;
}
//mpiexec -n 8 MPIApproach.exe in sourse file
//Debug the app after changes otherwise wont work