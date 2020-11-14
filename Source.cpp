#include <iostream>
#include <vector>
#include <map>
#include <iterator>
#include <cmath>
double _function(double& c, double& d, double&& x) {
    return (c * x + d);
}
double _function(double& c, double& d, const double& x) {
    return (c * x + d);
}
double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
double min_square(std::map<double, double>& _map, double& w1, double& w0) { //MNK
    double sum = 0;
    for (auto it : _map) {
        sum += (_function(w1, w0, it.first) - it.second) * (_function(w1, w0, it.first) - it.second);
    }
    return sum;
}
double find_min_dikh(std::map<double, double>& _map, double& a, double& b) { //a=cmin, b = cmax
    const double eps = 0.01;
    const double delta = 0.01;
    double xL, xR;
    double min_L, min_R;
    double w0 = 0;

    while ((b - a) > eps) {
        xL = 0.5 * (b + a) - delta;
        xR = 0.5 * (b + a) + delta;
        if (min_square(_map, xL, w0) < min_square(_map, xR, w0)) {
            b = xR;
        }
        else {
            a = xL;
        }
    } 
    if (min_square(_map, a, w0) < min_square(_map, b, w0)) {
        return a;
    } else {
        return b;
    }
}
double find_min_gold(std::map<double, double>& _map, double& a, double& b, double& w1) { //a=dmin, b = dmax
    const double t = 0.5 * (1 + std::sqrt(5));
    const double eps = 0.01;
    double xL = a + (1 - 1 / t) * b;
    double xR = a + b / t;
    while ((b - a) > eps) {
        if (min_square(_map, w1, xL) < min_square(_map, w1, xR)) {
            b = xR;
            xR = a + b - xL;
        }
        else {
            a = xL;
            xL = a + b - xR;
        }
        if (xL > xR) {
            std::swap(xL, xR);
        }
    } 
    if (min_square(_map, w1, a) < min_square(_map, w1, b)) {
        return a;
    } else {
        return b;
    }
}
int main() {
    double c = 1;
    double d = 0;
    double a = -2;
    double b = 2;
    int N = 24;
    int A = 2;
    double x_step = (b - a) / N;
    std::map<double, double> x_t; //for y =cx + d
    std::map<double, double> x_t_er; //for y with errors
    double error;
    for (auto i = 0; i < N + 1; ++i) {
        x_t.insert(std::make_pair((a + i * x_step), _function(c, d, (a + i * x_step))));
        error = A * fRand(-0.5, 0.5);
        x_t_er.insert(std::make_pair(((a + i * x_step) * error),
            _function(c, d, (a + i * x_step) * error)));
    }
    double cmin; // values c, d in f(x) = cx + d
    double cmax;
    double dmin;
    double dmax;
    std::map<double, double>::iterator it_x_t_begin = x_t.begin();
    std::map<double, double>::reverse_iterator it_x_t_end = x_t.rbegin();
    do {
        cmin = it_x_t_begin->second / it_x_t_begin->first;
    } while (it_x_t_begin->second == 0 && it_x_t_begin->first == 0);
    cmax = cmin;
    for (auto it : x_t) {
        if (it.first != 0 && cmin > it.second) {
            cmin = it.second / it.first;
        }
    }
    for (auto it : x_t) {
        if (it.first != 0 && cmax > it.second) {
            cmax = it.second / it.first;
        }
    }
    dmin = it_x_t_begin->first;
    dmax = it_x_t_end->first;
    c = find_min_dikh(x_t, cmin, cmax);
    d = find_min_gold(x_t, dmin, dmax, c);
    std::cout << "without errors: c=" << c << " d=" << d;
    std::map<double, double>::iterator it_x_t_er_begin = x_t_er.begin();
    std::map<double, double>::reverse_iterator it_x_t_er_end = x_t_er.rbegin();
    do {
    cmin = it_x_t_er_begin->second / it_x_t_er_begin->first;
    } while (it_x_t_er_begin->second == 0 && it_x_t_er_begin->first == 0);
    cmax = cmin;
    for (auto it : x_t_er) {
        if (it.first != 0 && cmin > it.second) {
            cmin = it.second / it.first;
        }
    }
    for (auto it : x_t_er) {
        if (it.first != 0 && cmax > it.second) {
            cmax = it.second / it.first;
        }
    }
    dmin = it_x_t_er_begin->first;
    dmax = it_x_t_er_end->second;
    c = find_min_dikh(x_t_er, cmin, cmax);
    d = find_min_gold(x_t_er, dmin, dmax, c);
    std::cout << "\nwith errors: c=" << c << " d=" << d;
    return 0;
}