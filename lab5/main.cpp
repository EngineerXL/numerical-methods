#include <bits/stdc++.h>

using namespace std;

using pii = pair<int, int>;

const double EPS = 1e-6;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);

    int n, K;
    cin >> n >> K;
    double l, T, h, tau;
    cin >> l >> T >> h >> tau;
    double a, b, c;
    cin >> a >> b >> c;
    vector< vector<double> > u(n, vector<double>(K));
    for (int i = 0; i < n; ++i) {
        cin >> u[i][0];
    }
    for (int k = 0; k < K; ++k) {
        cin >> u[0][k];
    }
    for (int k = 0; k < K; ++k) {
        cin >> u.back()[k];
    }

    double h2 = h * h;
    for (int k = 0; k < K - 1; ++k) {
        for (int i = 1; i < n - 1; ++i) {
            double uik = u[i][k];
            double dudx = (u[i + 1][k] - u[i - 1][k]) / (2 * h);
            double d2udx2 = (u[i - 1][k] - 2 * uik + u[i + 1][k]) / h2;
            u[i][k + 1] = uik + tau * (a * d2udx2 + b * dudx + c * uik);
        }
    }
    cout.precision(12);
    cout << fixed;
    for (int k = 0; k < K; ++k) {
        for (int i = 0; i < n; ++i) {
            if (i) {
                cout << ' ';
            }
            cout << u[i][k];
        }
        cout << '\n';
    }
}
