import numpy as np
import math
from scipy import special


def stirling_first_kind(n_max):
    # n corresponds to i, k corresponds to j below in compute_gregory_weights()
    # Create a (n_max+1) x (n_max+1) array initialized with zeros
    s = np.zeros((n_max + 1, n_max + 1), dtype=int)

    # Base case
    s[0, 0] = 1

    for n in range(1, n_max + 1):
        for k in range(1, n + 1):
            s[n, k] = s[n - 1, k - 1] - (n - 1) * s[n - 1, k]

    return s


# gregory_order corresponds to n0


def compute_coefficients(gregory_order):
    # Compute coefficients c_i.

    # The i-th entry of bernoulli_numbers corresponds to the i-th Bernoulli number B_i
    bernoulli_numbers = [1, -1/2, 1/6, 0, -1/30,
                         0, 1/42, 0, -1/30, 0, 5/66, 0, -691/2730]

    stirling_numbers = stirling_first_kind(12)
    # print(stirling_numbers)

    coefficients = []
    for i in range(1, gregory_order):
        # print("i: ", i)
        coefficients.append(0)
        for j in range(1, int(np.floor((i+1)/2))+1):
            # print("j: ", j)

            product_term = 1
            for l in range(2*j, i+1):
                # print("l: ", l)
                product_term *= l

            coefficients[i-1] += np.power(-1, i-2*j+1) * (bernoulli_numbers[2*j]
                                                          * stirling_numbers[i, 2*j-1]/(math.factorial(2*j)*product_term))

    return np.array(coefficients)


def compute_gregory_weights(gregory_order, n):
    # Compute Gregory weights based on coefficients c_i.
    coefficients = compute_coefficients(gregory_order)
    # Set weights of composite trapezoidal rule. Corrective terms will be added below.
    gregory_weights = np.ones(n+1)
    gregory_weights[0] = 0.5
    gregory_weights[-1] = 0.5

    # print(gregory_weights*12)

    # Add corrective terms.
    for i in range(1, gregory_order):
        # print("")
        # print("i: ", i)
        # print("coeff: ", coefficients[i-1]*24)

        for k in range(0, i+1):
            # print("k: ", k)
            # print((-1)**(i-k))
            # print(np.power(-1, k))
            # print(special.comb(i, k))

            # Forward difference terms for correction at beginning of integration interval.
            # print("added term: ",
            #       np.power(-1, i-k) * special.comb(i, k))
            gregory_weights[k] += coefficients[i-1] * \
                np.power(-1, i-k) * special.comb(i, k)

            # Backwards difference terms for correction at end of integration interval.
            # print("substracted term: ",
            #       np.power(-1, k) * special.comb(i, k))
            gregory_weights[-1-k] -= coefficients[i-1] * \
                np.power(-1, k) * special.comb(i, k)

    return gregory_weights


def main():
    gregory_order = 3
    coefficients = compute_coefficients(gregory_order)
    n = 8
    print(compute_gregory_weights(gregory_order, n)*24)


if __name__ == '__main__':
    main()
