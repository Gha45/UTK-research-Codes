import math

def compute_with_uncertainty(d, delta_d, A, delta_A, B):

    pi = math.pi
    f = 100 * (d *pi * B) / A

    # Error propagation
    term1 = (pi * B / A) * delta_d
    term2 = (pi * B * d / A**2) * delta_A
    print(term1, term2)
    delta_f = 100 * math.sqrt(term1**2 + term2**2)
    print(delta_f / 100)

    return f, delta_f

# Example usage
d = 10.1
delta_d = 0.2
A = 41.9
delta_A = 2.1


B = 0.10115


result, uncertainty = compute_with_uncertainty(d, delta_d, A, delta_A, B)
print(f"f = {result:.4f} ± {uncertainty:.4f}")
