#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

const int N = 10;
const double EPS = 1e-4;

double f(double x) {
    return x * x * x - 3 * x * x - 17 * x + 22 + sin(x);
}

/*vector<double> chebyshevNodes(double a, double b, int n) {
    vector<double> nodes(n);
    for (int i = 0; i < n; ++i) {
        nodes[i] = (a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * M_PI / (2 * n));
        cout << "Chebyshev node " << i + 1 << ": x = " << nodes[i] << endl;
    }
    return nodes;
}*/

vector<double> equidistantNodes(double a, double b, int n) {
    vector<double> nodes(n);
    for (int i = 0; i < n; ++i) {
        nodes[i] = a + i * (b - a) / (n - 1);
        cout << "Equidistant node " << i + 1 << ": x = " << nodes[i] << endl;
    }
    return nodes;
}

vector<vector<double>> dividedDifferences(const vector<double>& nodes) {
    int n = nodes.size();
    vector<vector<double>> dd(n, vector<double>(n));
    for (int i = 0; i < n; ++i) {
        dd[i][0] = f(nodes[i]);
    }
    for (int j = 1; j < n; ++j) {
        for (int i = 0; i < n - j; ++i) {
            dd[i][j] = (dd[i + 1][j - 1] - dd[i][j - 1]) / (nodes[i + j] - nodes[i]);
        }
    }
    return dd;
}

double newtonPolynomial(double x, const vector<double>& nodes, const vector<vector<double>>& dd) {
    double result = dd[0][0];
    double term = 1;
    for (int i = 1; i < nodes.size(); ++i) {
        term *= (x - nodes[i - 1]);
        result += term * dd[0][i];
    }
    return result;
}

vector<double> naturalCubicSplineCoefficients(const vector<double>& nodes, const vector<double>& values) {
    int n = nodes.size();
    vector<double> h(n - 1), alpha(n - 1), l(n), mu(n), z(n), c(n), d(n - 1);

    for (int i = 0; i < n - 1; ++i) {
        h[i] = nodes[i + 1] - nodes[i];
        alpha[i] = 3 * ((values[i + 1] - values[i]) / h[i] - (values[i] - values[i - 1]) / h[i - 1]);
    }

    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for (int i = 1; i < n - 1; ++i) {
        l[i] = 2 * (nodes[i + 1] - nodes[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n - 1] = 1;
    z[n - 1] = 0;
    c[n - 1] = 0;

    for (int j = n - 2; j >= 0; --j) {
        c[j]

                = z[j] - mu[j] * c[j + 1];
        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }

    vector<double> coefficients(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        coefficients[i] = (values[i + 1] - values[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3;
    }

    return coefficients;
}

double cubicSpline(double x, const vector<double>& nodes, const vector<double>& coefficients) {
    int i = 0;
    while (i < nodes.size() - 1 && x > nodes[i + 1]) {
        ++i;
    }
    double dx = x - nodes[i];
    double h = nodes[i + 1] - nodes[i];
    double a = f(nodes[i]);
    double b = (f(nodes[i + 1]) - f(nodes[i])) / h - h * (2 * coefficients[i] + coefficients[i + 1]) / 3;
    double c = coefficients[i];
    double d = (coefficients[i + 1] - coefficients[i]) / (3 * h);
    return a + b * dx + c * dx * dx + d * dx * dx * dx;
}

double findSmallestPositiveRoot(const vector<double>& nodes, const vector<vector<double>>& dd) {
    double a = 0.0;
    double b = 2.0;
    double c;
    while (b - a > EPS) {
        c = (a + b) / 2.0;
        if (newtonPolynomial(c, nodes, dd) * newtonPolynomial(a, nodes, dd) < 0) {
            b = c;
        } else {
            a = c;
        }
    }
    return c;
}

int main() {
    vector<string> terms = { "x", "f(x)", "Newton", "Cubic Spline" };
    double a = 0.0;
    double b = 2.0;

    vector<double> nodes = equidistantNodes(a, b, N);
    vector<vector<double>> dd = dividedDifferences(nodes);
    vector<double> coefficients = naturalCubicSplineCoefficients(nodes, dd[0]);

    double smallest_positive_root = findSmallestPositiveRoot(nodes, dd);

    cout << "\nSmallest positive root: " << smallest_positive_root << endl << endl;
    cout.width(20); cout << left << terms[0] << setw(20) << terms[1] << setw(20) << terms[2] << setw(20) << terms[3] << endl;

    for (int i = 0; i < nodes.size(); ++i) {
        double x = nodes[i];
        double fx = f(x);
        double newton = newtonPolynomial(x, nodes, dd);
        double spline = cubicSpline(x, nodes, coefficients);
        cout.width(20); cout << left << x << setw(20) << fx << setw(20) << newton << setw(20) << spline << endl;
    }

    return 0;
}
