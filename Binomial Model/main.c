#include <math.h>
#include <stdio.h>
#include <stdlib.h>


double option_price(double S0, double sigma, double K, double t, double r, double q, int N, int call, int american){
    double dt = t / N;
    double u = exp(sigma * sqrt(dt));
    double d = 1.0 / u;
    double p = (exp((r - q) * dt) - d) / (u - d);
    double discount = exp(-r * dt);

    double* option_values = malloc((N + 1) * sizeof(double));
    double* next_values = malloc((N + 1) * sizeof(double));
    double* stock = malloc((N + 1) * sizeof(double));

    double S = S0 * pow(d, N);
    double factor = u / d;

    for (int i = 0; i <= N; i++){
        stock[i] = S;
        option_values[i] = call ? fmax(S - K, 0.0) : fmax(K - S, 0.0);
        S *= factor;
    }

    for (int i = N - 1; i >= 0; i--){
        for (int j = 0; j <= i; j++){
            double value = discount * (p * option_values[j + 1] + (1 - p) * option_values[j]);
            stock[j] /= d;
            next_values[j] = american ? fmax((call ? fmax(stock[j] - K, 0) : fmax(K - stock[j], 0)), value) : value;
        }

        double* aux = option_values;
        option_values = next_values;
        next_values = aux;
    }

    double price = option_values[0];
    free(option_values);
    free(next_values);
    free(stock);

    return price;
}


int main(void){
    double S0 = 100.0;  // Spot price
    double sigma = 0.25;  // Volatility
    double K = 100.0;  // Strike price
    double t = 2.0;  // Time to maturity (years)
    double r = 0.04;  // Risk-free rate
    double q = 0.01;
    int N = 1000;  // Time steps

    double euro_put = option_price(S0, sigma, K, t, r, q, N, 0, 0);
    double euro_call = option_price(S0, sigma, K, t, r, q, N, 1, 0);
    double american_put = option_price(S0, sigma, K, t, r, q, N, 0, 1);
    double american_call = option_price(S0, sigma, K, t, r, q, N, 1, 1);

    printf("Price of european put: %.6f\n", euro_put);
    printf("Price of european call: %.6f\n", euro_call);
    printf("Price of american put: %.6f\n", american_put);
    printf("Price of american call: %.6f\n", american_call);

    return 0;
}