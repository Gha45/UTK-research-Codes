import numpy as np

# Bragg's law for d-spacing calculation
def bragg_d(two_theta_deg, wavelength):
    theta_rad = np.radians(two_theta_deg / 2)
    return wavelength / (2 * np.sin(theta_rad))

# Generate allowed FCC reflections and their (h²+k²+l²) sums
def generate_fcc_reflections(max_index):
    allowed_sums = set()
    hkl_dict = {}
    for h in range(max_index+1):
        for k in range(max_index+1):
            for l in range(max_index+1):
                if (h == k == l == 0):
                    continue
                if (h % 2 == k % 2 == l % 2):
                    hkl_sum = h**2 + k**2 + l**2
                    allowed_sums.add(hkl_sum)
                    if hkl_sum in hkl_dict:
                        hkl_dict[hkl_sum].append((h,k,l))
                    else:
                        hkl_dict[hkl_sum] = [(h,k,l)]
    return sorted(list(allowed_sums)), hkl_dict

# Experimental data (replace with your measured angles)
two_theta_deg = [9.39	,	10.821	,	15.31	,	18.005]
wavelength = 0.4132
# Calculate experimental d-spacings
d_obs = np.array([bragg_d(angle, wavelength) for angle in two_theta_deg])
# Generate FCC reflections
hkl_sums, hkl_dict = generate_fcc_reflections(max_index=5)
# Matching experimental d-spacings to FCC reflections
matched_hkls = []
a_values = []
for d in d_obs:
    best_match = None
    smallest_diff = float('inf')
    for hkl_sum in hkl_sums:
        a_calc = d * np.sqrt(hkl_sum)
        diff = abs(a_calc - np.mean(a_values)) if a_values else 0
        if diff < smallest_diff:
            smallest_diff = diff
            best_match = hkl_sum
            best_a = a_calc
    matched_hkls.append((best_match, hkl_dict[best_match]))
    a_values.append(best_a)
# Average lattice parameter
avg_a = np.mean(a_values)
std_a = np.std(a_values)
# Results
print("Matched FCC reflections and lattice parameters:")
for i, (hkl_sum, hkls) in enumerate(matched_hkls):
    print(f"d_obs: {d_obs[i]:.4f} Å, (h²+k²+l²): {hkl_sum}, possible hkls: {hkls}, a = {a_values[i]:.4f} Å")
print(f"\nEstimated average lattice parameter: a = {avg_a:.4f} ± {std_a:.4f} Å")