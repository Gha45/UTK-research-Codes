import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.interpolate import interp1d
import scipy
file_path = '/Users/georgeadamson/Desktop/UTK python programs for reaserch/946mevAu_(YbTbGdErDy)2Ti2O7.csv'
df = pd.read_csv(file_path,skiprows=1, header=None)

x_value = 12.5  # Replace with the x-coordinate where you want the vertical line
interpolated_y = np.interp(x_value, df[0], df[1])
new_point = pd.DataFrame({0: [x_value], 1: [interpolated_y]})
df1 = pd.concat([df, new_point]).sort_values(by=0).reset_index(drop=True)

###
filtered_df1 = df1[df1[0] < x_value].copy()

# Interpolating function for the entire dataset
f = interp1d(filtered_df1[0], filtered_df1[1], kind='linear')
# Integrate the area under the curve up to fill_to_x
x_values = filtered_df1[0]
y_values = filtered_df1[1]

area_under_curve = simps(y_values, x_values)
# Calculate the length of the curve up to fill_to_x
length_of_curve = max(x_values)
# Calculate the desired averaged value for the energy loss
#desired_value = area_under_curve / length_of_curve
desired_value = area_under_curve/length_of_curve

# error calculated with standard error fo the mean

std_deviation = scipy.stats.sem(y_values)

###

plt.figure(figsize=(10, 10))
plt.tight_layout()
plt.subplots_adjust(bottom=0.15)
# plotting the curve
plt.plot(df.iloc[:,0], df.iloc[:,1],color='maroon',linewidth=2.5)
plt.axvline(x=x_value, color='black', linestyle='--', label='Vertical Line',linewidth=2.5)
plt.fill_between(df1[0], df1[1], where=(df1[0] <= x_value ), color='skyblue', alpha=0.4)
###
# Adding a double-sided arrow and text "thickness"
plt.annotate('', xy=(x_value , max(y_values)/2), xytext=(0, max(y_values)/2),
             arrowprops=dict(arrowstyle='<->', color='black', lw=2.5))
plt.text(x_value / 2, ((max(y_values)/2)+2)*1.015, 'Sample', verticalalignment='bottom', horizontalalignment='center',
         color='black', fontsize=25)
plt.text(x_value / 2, (max(y_values)/2)*1.015, 'Thickness', verticalalignment='bottom', horizontalalignment='center',
         color='black', fontsize=25)
#plt.text(x_value / 2, (max(y_values)/2)*1.015, 'Sample Thickness', verticalalignment='bottom', horizontalalignment='center',
#         color='black', fontsize=16)
###
plt.text(15,5,"$(YbTbGdErDy)_2Ti_2O_7$", fontsize=30, color='black')
plt.text(15,1,f" {round(desired_value,2)} ± {round(std_deviation,2)} Kev/nm",ha='left', fontsize=30, color='black')

plt.text(0.5, (max(y_values))*1.125, '946 Mev Au',
         color='maroon', fontsize=30)

plt.xlabel(r"Depth ($\mu m$)", fontsize=30)
plt.ylabel("Energy loss (keV/nm)", fontsize=30)
plt.xlim(left=0)
plt.ylim(bottom=0,top = (max(y_values)+5))
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
#plt.grid(False)
#plt.savefig("(YbTbGdErDy)2Ti2O7_energyloss.svg")
plt.show()
plt.savefig('squareloss.svg')
#####
# lets do ionic speed loss
ions_table= {'au':183504.318,'gold':183504.318,'xe':122305.95,'xenon':122305.95, 'ar':37213.1853,'argon':37213.1853,'u':221742.1467,'uranium':221742.1467}

ion =input('What Ion was chosen:').strip().lower()
ion_charge = float(input('What is the charge of that ion:').strip())
chosen_ion = float(ions_table.get(ion))

if chosen_ion is not None:
    print(f"The mass of {ion} is {chosen_ion} MeV.")
else:
    print("Ion not found. Please check your input.")

MCI = (chosen_ion- (ion_charge * 0.511))*1000

# equation based off relativistic kinetic energy ot calculate beta
def ion_speed(MCI,T):
    beta = (1-((1/((T/MCI)+1)))**2)**(0.5)
    return beta
# Calculate ion speed for each row
df[5] = df.iloc[:,4].apply(lambda x: ion_speed(MCI, x))

interpolated_y2 = np.interp(x_value, df[0], df[5])
new_point2 = pd.DataFrame({0: [x_value], 5: [interpolated_y2]})
df2 = pd.concat([df, new_point2]).sort_values(by=0).reset_index(drop=True)
###

filtered_df2 = df2[df2[0] < x_value].copy()

# Interpolating function for the entire dataset
f2 = interp1d(filtered_df2[0], filtered_df2[5], kind='linear')
# Integrate the area under the curve up to fill_to_x
x_values2 = filtered_df2[0]
y_values2 = filtered_df2[5]
area_under_curve2 = simps(y_values2, x_values2)
# Calculate the length of the curve up to fill_to_x
length_of_curve2 =max(x_values2)
# Calculate the desired averaged value for ion speed
#beta = area_under_curve2 / length_of_curve2
beta = area_under_curve2/length_of_curve2
# error calculated with standard error fo the mean

error_beta = scipy.stats.sem(y_values2)
###
plt.figure(figsize=(10, 6))
plt.tight_layout()
plt.subplots_adjust(bottom=0.15)
plt.plot(df.iloc[:,0], df.iloc[:,5],color='maroon',linewidth=2.5)
plt.axvline(x=x_value, color='black', linestyle='--', label='Vertical Line',linewidth=2.5)
plt.fill_between(df2[0], df2[5], where=(df2[0] <= x_value ), color='skyblue', alpha=0.4)
###
# Adding a double-sided arrow and text "thickness"
plt.annotate('', xy=(x_value , max(y_values2)/2), xytext=(0, max(y_values2)/2),
             arrowprops=dict(arrowstyle='<->', color='black', lw=2.5))
plt.text(x_value / 2, (max(y_values2)/2)*1.015, 'Sample Thickness', verticalalignment='bottom', horizontalalignment='center',
         color='black', fontsize=16)
###
plt.text(1.2,0.017,"$(YbTbGdErDy)_2Ti_2O_7$", fontsize=18, color='black')
plt.text(1,0.01,f" Inital:{round(filtered_df2[5].iloc[0],4)} ",ha='left', fontsize=18, color='black')
plt.text(1,0.003,f" Exiting:{round(filtered_df2[5].iloc[-1],4)} ± {round(error_beta,4)} ",ha='left', fontsize=18, color='black')

plt.text(1, (max(y_values2))*1.015, '946 Mev Au',
         color='maroon', fontsize=20)

plt.xlabel(r"Depth ($\mu m$)", fontsize=20)
plt.ylabel("Ion speed(beta)", fontsize=20)
plt.xlim(left=0)
plt.ylim(bottom=0,top = (max(y_values2)))
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
#plt.grid(False)
#plt.savefig("(YbTbGdErDy)2Ti2O7_energyloss.svg")
plt.show()

### energy density
df[6]= (df[1])/(np.pi*df[5])

### Interpolation
interpolated_y3 = np.interp(x_value, df[0], df[6])
new_point3 = pd.DataFrame({0: [x_value], 6: [interpolated_y3]})
df3 = pd.concat([df, new_point3]).sort_values(by=0).reset_index(drop=True)
###

filtered_df3 = df3[df3[0] <= x_value].copy()

# Interpolating function for the entire dataset
f3 = interp1d(filtered_df3[0], filtered_df3[5], kind='linear')
# Integrate the area under the curve up to fill_to_x
x_values3 = filtered_df3[0]
y_values3 = filtered_df3[6]

area_under_curve3 = simps(y_values3, x_values3)
# Calculate the length of the curve up to fill_to_x
length_of_curve3 = max(x_values3)
# Calculate the desired averaged value of energy density
desired_value3 = area_under_curve3 / length_of_curve3
# did error propegation of the error
std_deviation3 =np.sqrt(((std_deviation/(np.pi*beta))**2)+(((-desired_value*error_beta)/(np.pi*(beta)**2))**2))
print(desired_value)
print(std_deviation)
print(beta)
print(error_beta)
print(desired_value3)
print(std_deviation3)
### graphing
plt.figure(figsize=(10, 6))
plt.tight_layout()
plt.subplots_adjust(bottom=0.15)
plt.plot(df.iloc[:,0], df.iloc[:,6],color='maroon',linewidth=2.5)
plt.axvline(x=x_value, color='black', linestyle='--', label='Vertical Line',linewidth=2.5)
plt.fill_between(df3[0], df3[6], where=(df3[0] <= x_value ), color='skyblue', alpha=0.4)
###
# Adding a double-sided arrow and text "thickness"
plt.annotate('', xy=(x_value , max(y_values3)/2), xytext=(0, max(y_values3)/2),
             arrowprops=dict(arrowstyle='<->', color='black', lw=2.5))
plt.text(x_value / 2, (max(y_values3)/2)*1.015, 'Sample Thickness', verticalalignment='bottom', horizontalalignment='center',
         color='black', fontsize=16)
###
plt.text(1,30,"$(YbTbGdErDy)_2Ti_2O_7$", fontsize=20, color='black')
plt.text(1,15,f" {round(desired_value3,2)} ± {round(std_deviation3,2)} Kev/nm",ha='left', fontsize=20, color='black')

plt.text(1, (max(y_values3))*1.015, '946 Mev Au',
         color='maroon', fontsize=20)

plt.xlabel(r"Depth ($\mu m$)", fontsize=20)
plt.ylabel("energy density(kev/nm)", fontsize=20)
plt.xlim(left=0)
plt.ylim(bottom=0,top = (max(df.iloc[:,6]))+10)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
#plt.grid(False)
#plt.savefig("(YbTbGdErDy)2Ti2O7_energyloss.svg")
plt.show()