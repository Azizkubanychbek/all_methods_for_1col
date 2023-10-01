#include <iostream>
#include <cmath>
#include <vector>

using namespace std;


// Уравнение: cos(4 * x) - x + 0.5
double equation(double x) {
    return cos(4 * x) - x + 0.5;
}

// Производная уравнения: -4 * sin(4 * x) - 1
double derivative(double x) {
    return -4 * sin(4 * x) - 1;
}

// Метод Ньютона
double newtonMethod(double x0, double epsilon, int maxIterations) {
    int iteration = 0;
    double x = x0;

    while (iteration < maxIterations) {
        double fx = equation(x);
        double dfx = derivative(x);

        if (abs(dfx) < epsilon) {
            cerr << "Division by zero (Newton's method)" << endl;
            return NAN; // Защита от деления на ноль
        }

        x = x - fx / dfx;

        if (abs(fx) < epsilon) {
            cout << "Newton's method converged after " << iteration << " iterations" << endl;
            return x; // Найдено приближенное решение
        }

        cout << "Iteration " << iteration << ": x = " << x << endl;
        iteration++;
    }

    cerr << "Newton's method did not converge" << endl;
    return NAN;
}

// Метод итераций (простая итерация)
double iterationMethod(double x0, double epsilon, int maxIterations) {
    int iteration = 0;
    double x = x0;

    while (iteration < maxIterations) {
        x = equation(x);

        cout << "Iteration " << iteration << ": x = " << x << endl;

        if (abs(equation(x)) < epsilon) {
            cout << "Iteration method converged after " << iteration << " iterations" << endl;
            return x; // Найдено приближенное решение
        }

        iteration++;
    }

    cerr << "Iteration method did not converge" << endl;
    return NAN;
}

// Метод деления пополам
double bisectionMethod(double a, double b, double epsilon, int maxIterations) {
    int iteration = 0;

    while (iteration < maxIterations) {
        double c = (a + b) / 2;
        double fa = equation(a);
        double fb = equation(b);
        double fc = equation(c);

        if (abs(fc) < epsilon) {
            cout << "Bisection method converged after " << iteration << " iterations" << endl;
            return c; // Найдено приближенное решение
        }

        if (fa * fc < 0) {
            b = c;
        }
        else {
            a = c;
        }

        cout << "Iteration " << iteration << ": x = " << c << endl;
        iteration++;
    }

    cerr << "Bisection method did not converge" << endl;
    return NAN;
}
// Метод хорд
double chordMethod(double a, double b, double epsilon, int maxIterations) {
    int iteration = 0;
    double x0 = a;

    while (iteration < maxIterations) {
        double fa = equation(a);
        double fb = equation(b);
        double x1 = (a * fb - b * fa) / (fb - fa);

        if (abs(equation(x1)) < epsilon) {
            cout << "Chord method converged after " << iteration << " iterations" << endl;
            return x1; // Найдено приближенное решение
        }

        a = x0;
        b = x1;

        cout << "Chord Iteration " << iteration << ": x1 = " << x1 << endl;
        iteration++;
    }

    cerr << "Chord method did not converge" << endl;
    return NAN;
}

// Метод секущих
double secantMethod(double x0, double x1, double epsilon, int maxIterations) {
    int iteration = 0;

    while (iteration < maxIterations) {
        double f0 = equation(x0);
        double f1 = equation(x1);

        if (abs(f1 - f0) < epsilon) {
            cout << "Secant method converged after " << iteration << " iterations" << endl;
            return x1; // Найдено приближенное решение
        }

        double x2 = x1 - (f1 * (x1 - x0)) / (f1 - f0);

        x0 = x1;
        x1 = x2;

        cout << "Secant Iteration " << iteration << ": x2 = " << x2 << endl;
        iteration++;
    }

    cerr << "Secant method did not converge" << endl;
    return NAN;
}
//метод Гауса
void printMatrix(vector<vector<double>>& matrix) {
    int n = matrix.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n + 1; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

// Вычисление определителя матрицы
double determinant(vector<vector<double>>& matrix) {
    int n = matrix.size();
    if (n == 1) {
        return matrix[0][0];
    }
    else {
        double det = 0;
        for (int i = 0; i < n; i++) {
            vector<vector<double>> submatrix(n - 1, vector<double>(n - 1));
            for (int j = 1; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    if (k < i) {
                        submatrix[j - 1][k] = matrix[j][k];
                    }
                    else if (k > i) {
                        submatrix[j - 1][k - 1] = matrix[j][k];
                    }
                }
            }
            det += (i % 2 == 0 ? 1 : -1) * matrix[0][i] * determinant(submatrix);
        }
        return det;
    }
}

// Метод Гаусса для решения системы линейных уравнений и вычисления определителя
pair<vector<double>, double> gauss(vector<vector<double>>& matrix) {
    int n = matrix.size();
    vector<double> result(n);

    // Прямой ход
    for (int i = 0; i < n; i++) {
        cout << "Итерация " << i + 1 << ":" << endl;
        printMatrix(matrix);

        // Нормализация текущей строки
        double maxEl = abs(matrix[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(matrix[k][i]) > maxEl) {
                maxEl = abs(matrix[k][i]);
                maxRow = k;
            }
        }

        for (int k = i; k < n + 1; k++) {
            swap(matrix[maxRow][k], matrix[i][k]);
        }

        // Обнуление нижележащих строк
        for (int k = i + 1; k < n; k++) {
            double factor = -matrix[k][i] / matrix[i][i];
            for (int j = i; j < n + 1; j++) {
                if (i == j) {
                    matrix[k][j] = 0;
                }
                else {
                    matrix[k][j] += factor * matrix[i][j];
                }
            }
        }
    }

    // Вычисление определителя
    double det = 1.0;
    for (int i = 0; i < n; i++) {
        det *= matrix[i][i];
    }

    // Обратный ход
    for (int i = n - 1; i >= 0; i--) {
        result[i] = matrix[i][n] / matrix[i][i];
        for (int k = i - 1; k >= 0; k--) {
            matrix[k][n] -= matrix[k][i] * result[i];
        }
    }

    return make_pair(result, det);
}

// Метод Крамера для решения системы линейных уравнений
vector<double> cramer(vector<vector<double>>& matrix) {
    int n = matrix.size();
    vector<double> result(n);
    double detA = determinant(matrix);

    for (int i = 0; i < n; i++) {
        vector<vector<double>> tempMatrix = matrix;
        for (int j = 0; j < n; j++) {
            tempMatrix[j][i] = matrix[j][n];
        }
        double detAi = determinant(tempMatrix);
        result[i] = detAi / detA;
    }

    return result;
}

int main() {
    double x1 = 1.4;
    double x2 = 1.7;
    double epsilon = 1e-3;
    int maxIterations = 1000;

    double resultNewton = newtonMethod(x1, epsilon, maxIterations);
    double resultIteration = iterationMethod(x1, epsilon, maxIterations);
    double resultBisection = bisectionMethod(x1, x2, epsilon, maxIterations);
    double resultChord = chordMethod(0, 1, epsilon, maxIterations);
    double resultSecant = secantMethod(x1, x2, epsilon, maxIterations);

    cout << "Newton's method result: " << resultNewton << endl;
    cout << "Iteration method result: " << resultIteration << endl;
    cout << "Bisection method result: " << resultBisection << endl;
    cout << "Chord method result: " << resultChord << endl;
    cout << "Secant method result: " << resultSecant << endl;

    vector<vector<double>> matrix = {
         {-7, 3, 7, -37},
         {5, 8, -1,-59},
         {-7, -4, 4, 59}
    };


    cout << "Метод Гаусса:" << endl;
    pair<vector<double>, double> gaussResult = gauss(matrix);
    cout << "Решение:" << endl;
    for (int i = 0; i < gaussResult.first.size(); i++) {
        cout << "x[" << i << "] = " << gaussResult.first[i] << endl;
    }
    cout << "Определитель A: " << gaussResult.second << endl;

    cout << endl;

    cout << "Метод Крамера:" << endl;
    vector<double> cramerResult = cramer(matrix);
    cout << "Решение:" << endl;
    for (int i = 0; i < cramerResult.size(); i++) {
        cout << "x[" << i << "] = " << cramerResult[i] << endl;
    }


    return 0;
}
