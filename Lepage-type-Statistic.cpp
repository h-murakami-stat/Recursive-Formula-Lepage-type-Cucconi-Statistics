#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <map>
#include <tuple>
#include <algorithm>
#include <numeric>

using namespace Rcpp;

// Function of k(x)
inline double kk(double x) {
    return x * (x * x - 1) / 12.0;
}

// Lepage-type Statistic
inline double M_Lepage(int m1, int m2, double W, double MO) {
    int N = m1 + m2;
    double e_w = m1 * (m1 + m2 + 1) / 2.0;
    double v_w = m1 * m2 * (m1 + m2 + 1) / 12.0;

    double e_mo = m1 * (N * N - 1) / 12.0;
    double v_mo = m1 * m2 * (N + 1) * (N * N - 4) / 180.0;

    double A = (W - e_w) / sqrt(v_w);
    double B = (MO - e_mo) / sqrt(v_mo);

    return A * A + B * B;
}

// Probability of Wilcoxon and Mood statistics
std::map<std::tuple<int, int, int, double>, double> memo;

double MLep_rec(int m1, int m2, int W, double MO) {
    if (m1 < 0 || m2 < 0 || W < 0 || MO < 0) {
        return 0;
    }

    auto key = std::make_tuple(m1, m2, W, MO);
    if (memo.find(key) != memo.end()) {
        return memo[key];
    }

    if (m1 == 0 && W == 0 && MO == 0) {
        return memo[key] = 1;
    }

    if (m2 == 0 && W == m1 * (m1 + 1) / 2 && MO == kk(m1)) {
        return memo[key] = 1;
    }

    if (m1 == 1 && m2 == 1 && (W == 1 || W == 2) && MO == 0.25) {
        return memo[key] = 0.5;
    }

    if (m1 == 1 && m2 >= 2) {
        int N1 = 1 + m2;
        double pp = 0;
        for (int i = 1; i <= N1; ++i) {
            if (W == i && MO == std::pow(i - (N1 + 1) / 2.0, 2)) {
                pp += 1.0 / N1;
            }
        }
        return memo[key] = pp;
    }

    if (m1 >= 2 && m2 == 1) {
        int N2 = m1 + 1;
        double pp2 = 0;
        for (int i = 1; i <= N2; ++i) {
            if (W == (m1 + 1) * (m1 + 2) / 2 - i && MO == kk(N2) - std::pow(i - (N2 + 1) / 2.0, 2)) {
                pp2 += 1.0 / N2;
            }
        }
        return memo[key] = pp2;
    }

    if (m1 >= 2 && m2 >= 2) {
        int N3 = m1 + m2;
        double result = m1 * (m1 - 1) / (static_cast<double>(N3) * (N3 - 1)) *
                        MLep_rec(m1 - 2, m2, W - 2 * m1 - m2 + 1, MO - 2 * std::pow((N3 - 1) / 2.0, 2)) +
                        m1 * m2 / (static_cast<double>(N3) * (N3 - 1)) *
                        MLep_rec(m1 - 1, m2 - 1, W - m1, MO - std::pow((N3 - 1) / 2.0, 2)) +
                        m2 * (m2 - 1) / (static_cast<double>(N3) * (N3 - 1)) *
                        MLep_rec(m1, m2 - 2, W - m1, MO) +
                        m1 * m2 / (static_cast<double>(N3) * (N3 - 1)) *
                        MLep_rec(m1 - 1, m2 - 1, W - 2 * m1 - m2 + 1, MO - std::pow((N3 - 1) / 2.0, 2));
        return memo[key] = result;
    }

    return memo[key] = 0;
}

// [[Rcpp::export]]
Rcpp::List Lepage_Type_r(int m1, int m2, double W, double MO) {
    double TT = M_Lepage(m1, m2, W, MO);
    double P = MLep_rec(m1, m2, W, MO);
    return Rcpp::List::create(
        Rcpp::Named("TT") = TT,
        Rcpp::Named("P") = P
    );
}

// Function of Sum of the first n Sequences When N is Even
// [[Rcpp::export]]
int sum_sequence1(int n) {
    int sum = 0;
    int value;
    for (int i = 1; i <= n; ++i) {
        if (i % 2 == 0) {
            value = i-1;
        }else{
            value = i;
        }
        sum += value * value;
    }
    return sum;
}



// Function of Sum of the first n Sequences When N is Odd
// [[Rcpp::export]]
int sum_sequence2(int n) {
    int sum = 0;
    int value;
    for (int i = 1; i <= n; ++i) {
        if (i % 2 == 0) {
            value = i / 2;
        } else {
            value = (i - 1) / 2;
        }
        sum += value * value;
    }
    return sum;
}



