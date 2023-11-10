#include <iostream>
#include <random> // rand(), srand() functions
#include <time.h> // clock(), time() functions
using namespace std;

#define N 2       // Matrix dimensions
#define EPSILON 1e-15
#define MIN -9
#define MAX 9

void create_randox_matrix(double A[N][N + 1]) {
    srand(time(NULL)); // set seed for random as current date
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N+1; j++) {
            A[i][j] = MIN + (double)(rand()) / ((double)(RAND_MAX / (MAX - MIN))); // inputs random double from MIN to MAX into matrix
        }
    }
}

void copy_matrix(double A[N][N + 1], double copy[N][N + 1]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N + 1; j++) {
            copy[i][j] = A[i][j];
        }
    }
}

void check_answers(double original_matrix[N][N + 1], double X[N]) {
    double sum;
    for (int i = 0; i < N; i++) {
        sum = 0;
        for (int j = 0; j < N; j++) {
            sum += original_matrix[i][j]*X[j];
        }
        cout << sum << " =?= " << original_matrix[i][N] << endl;
    }
}

void print_matrix(double A[N][N + 1]) {
    cout.precision(2);
    for (int i = 0; i < N; i++){
        cout << fixed << endl;
        for (int j = 0; j <= N; j++){
            cout << A[i][j] << "\t";
        }
    }
    cout << endl;
}

void print_solution_vector(double X[N]) {
    cout << "The solution for the system:" << endl;
    for (int i = 0; i < N; i++) {
        cout.precision(2);
        cout << fixed << "x" << i + 1 << " = " << X[i] << "\t";
        if (i % 15 == 0) {
            cout << endl;
        }
    }
    cout << endl;
}

void swap_row(double A[N][N + 1], int i, int j) {
    for (int k = 0; k <= N; k++) {
        double temp = A[i][k];
        A[i][k] = A[j][k];
        A[j][k] = temp;
    }
}

int forward_elimination(double A[N][N + 1]) {
    for (int k = 0; k < N; k++) {
        int index_of_max_element_in_row = k;
        int value_of_max_element_in_row = A[index_of_max_element_in_row][k];
        for (int i = k + 1; i < N; i++)
            if (abs(A[i][k]) > value_of_max_element_in_row)
                value_of_max_element_in_row = A[i][k], index_of_max_element_in_row = i;
        if (abs(A[k][index_of_max_element_in_row])<EPSILON) {
            return k;
            }
        if (index_of_max_element_in_row != k) {
            swap_row(A, k, index_of_max_element_in_row);
        }
        for (int i = k + 1; i < N; i++) {
            double c = A[i][k] / A[k][k];
            for (int j = k + 1; j <= N; j++)
                A[i][j] -= c * A[k][j];
            A[i][k] = 0;
        }
    }
    return -1;
}

void back_substitution(double A[N][N + 1], double X[N]) {
    for (int i = N - 1; i >= 0; i--) {
        X[i] = A[i][N];
        for (int j = i + 1; j < N; j++) {
            X[i] -= A[i][j] * X[j];
        }
        X[i] = X[i] / A[i][i];
    }
}

void Gauss(double A[N][N + 1], double X[N]) {
    int ld_flag = forward_elimination(A);
    if (ld_flag != -1) {
        cout << "At least two rows of the matrix are linearly dependent";
        return;
    }
    back_substitution(A, X);
}

int main() {
    /*double A[N][N + 1] = { {2, 1, 1, 2},
                          {1, -1, 0, -2},
                          {3, -1, 2, 2} };
                          */
    //double A[N][N + 1]; // fails with stack overflow (16kb size) error for huge matrices!
    auto A = new double[N][N+1];
    auto A_copy = new double[N][N + 1];
    auto X = new double[N];
    create_randox_matrix(A);
    copy_matrix(A, A_copy);
    cout << "Solving matrix:" << endl;
    print_matrix(A_copy);
    clock_t time_start = clock();
    Gauss(A, X);
    clock_t time_finish = clock();
    print_solution_vector(X);
    //cout << "Got solution for the matrix with Gauss method in " << (double)(time_finish - time_start) / CLOCKS_PER_SEC << " seconds" << endl;
    cout << "Checking answers:" << endl;
    check_answers(A_copy, X);
    return 0;
}