/**
 * @file ShapeFunction.h
 * @author k.ueda
 * @date July, 2024
*/

#ifndef SPLINE_H
#define SPLINE_H

class Spline
{
    public:
        void setPoints(const std::vector<double>& x, const std::vector<double>& y)
        {
            this->x = x;
            this->y = y;
            n = x.size();

            std::vector<double> h(n - 1), alpha(n - 1);
            for (int i = 0; i < n - 1; ++i) {
                h[i] = x[i + 1] - x[i];
                alpha[i] = (i == 0) ? 0 : (3 / h[i]) * (y[i + 1] - y[i]) - (3 / h[i - 1]) * (y[i] - y[i - 1]);
            }

            std::vector<double> l(n), mu(n), z(n);
            l[0] = 1.0;
            mu[0] = 0.0;
            z[0] = 0.0;

            for(int i = 1; i < n - 1; ++i) {
                l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
                mu[i] = h[i] / l[i];
                z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
            }

            l[n - 1] = 1.0;
            z[n - 1] = 0.0;
            c.resize(n);
            b.resize(n - 1);
            d.resize(n - 1);
            c[n - 1] = 0.0;

            for(int j = n - 2; j >= 0; --j) {
                c[j] = z[j] - mu[j] * c[j + 1];
                b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
                d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
            }
        }

        double operator()(double x_val) const
        {
            if(x_val < x.front() || x_val > x.back()) {
                throw std::out_of_range("x_val is out of range");
            }

            int i = 0;
            int j = n - 1;
            while(i < j){
                int k = (i + j) / 2;
                if(x_val <= x[k]){
                    j = k;
                }else{
                    i = k + 1;
                }
            }
            --i;

            double dx = x_val - x[i];
            return y[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
    }

    private:
        int n;
        std::vector<double> x, y, b, c, d;
};

#endif