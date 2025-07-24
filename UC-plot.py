import matplotlib.pyplot as plt


pressure = [1.7 ,
10.18,
18.81 ,
24.84,
30.29 ,
34.34 ,]


a_values = ([ 4.5257,
4.4702,
4.4304,
4.4174,
4.3963,
4.38])


a_errors =([ 0.0031,
0.0011,
0.013,
0.0027,
0.0011,
0.0054])

plt.figure(figsize=(8,6))
plt.errorbar(pressure, (a_values), yerr=(a_errors), fmt='o', ecolor='r', capsize=5, linestyle='-', marker='s', color='blue', markersize=7)
plt.xlabel('Pressure (GPa)', fontsize=12)
plt.ylabel('Lattice Parameter (Å)', fontsize=12)
plt.title('Lattice Parameter vs Pressure for ML8', fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.show()
