#ifndef SRC_HPP
#define SRC_HPP

#include "fraction.hpp"
#include <vector>

class resistive_network {
private:
    int n, m;
    std::vector<int> u, v;
    std::vector<fraction> r, g; // resistance and conductance
    std::vector<std::vector<fraction>> L; // Laplacian

    static bool is_zero(const fraction &x) { return x == fraction(0); }

    // Solve L_red * x = b_red with node `ground` fixed at 0 voltage.
    // Returns full vector U of size n with U[ground]=0.
    std::vector<fraction> solve_with_ground(const std::vector<std::vector<fraction>> &Lap,
                                            const std::vector<fraction> &b,
                                            int ground) {
        int N = (int)Lap.size();
        std::vector<int> idx_map; idx_map.reserve(N-1);
        for (int i = 0; i < N; ++i) if (i != ground) idx_map.push_back(i);
        int K = N - 1;
        std::vector<std::vector<fraction>> A(K, std::vector<fraction>(K));
        std::vector<fraction> bb(K);
        for (int i = 0; i < K; ++i) {
            int ii = idx_map[i];
            bb[i] = b[ii];
            for (int j = 0; j < K; ++j) {
                int jj = idx_map[j];
                A[i][j] = Lap[ii][jj];
            }
        }
        // Gaussian elimination
        for (int col = 0; col < K; ++col) {
            int pivot = col;
            while (pivot < K && is_zero(A[pivot][col])) pivot++;
            if (pivot == K) throw resistive_network_error();
            if (pivot != col) {
                std::swap(A[pivot], A[col]);
                std::swap(bb[pivot], bb[col]);
            }
            fraction diag = A[col][col];
            for (int row = col + 1; row < K; ++row) {
                if (is_zero(A[row][col])) continue;
                fraction factor = A[row][col] / diag;
                for (int c = col; c < K; ++c) A[row][c] = A[row][c] - factor * A[col][c];
                bb[row] = bb[row] - factor * bb[col];
            }
        }
        std::vector<fraction> x(K);
        for (int i = K - 1; i >= 0; --i) {
            fraction sum = bb[i];
            for (int j = i + 1; j < K; ++j) sum = sum - A[i][j] * x[j];
            x[i] = sum / A[i][i];
        }
        std::vector<fraction> U(N, fraction(0));
        for (int i = 0; i < K; ++i) U[idx_map[i]] = x[i];
        U[ground] = fraction(0);
        return U;
    }

public:
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[]) {
        n = interface_size_; m = connection_size_;
        u.resize(m); v.resize(m); r.resize(m); g.resize(m);
        L.assign(n, std::vector<fraction>(n, fraction(0)));
        for (int k = 0; k < m; ++k) {
            int a = from[k] - 1; int b = to[k] - 1;
            u[k] = a; v[k] = b; r[k] = resistance[k];
            g[k] = fraction(1) / r[k];
            L[a][a] = L[a][a] + g[k];
            L[b][b] = L[b][b] + g[k];
            L[a][b] = L[a][b] - g[k];
            L[b][a] = L[b][a] - g[k];
        }
    }

    ~resistive_network() = default;

    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        int a = interface_id1 - 1, b = interface_id2 - 1;
        std::vector<fraction> I(n, fraction(0));
        I[a] = I[a] + fraction(1);
        I[b] = I[b] - fraction(1);
        int ground = n - 1;
        std::vector<fraction> U = solve_with_ground(L, I, ground);
        return U[a] - U[b];
    }

    fraction get_voltage(int id, fraction current[]) {
        std::vector<fraction> I(n, fraction(0));
        for (int i = 0; i < n; ++i) I[i] = current[i];
        int ground = n - 1;
        std::vector<fraction> U = solve_with_ground(L, I, ground);
        return U[id - 1];
    }

    fraction get_power(fraction voltage[]) {
        // Sum over edges: (u_i - u_j)^2 / r_ij = g_ij * (u_i - u_j)^2
        fraction P(0);
        for (int k = 0; k < m; ++k) {
            fraction du = voltage[u[k]] - voltage[v[k]];
            P = P + g[k] * du * du;
        }
        return P;
    }
};

#endif // SRC_HPP
