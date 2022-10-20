# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import numpy as np
import matplotlib.pyplot as plt

Nominal_Power = np.arange(-57,-12,1) #dBm
V_trans = 5.65 #V
I_trans = 0.13 #A

V_rec = 12 #V
I_rec = 0.14 #A

Marker = 425.057 #MHz
Gain_1 = np.array([15.20,15.10,15.09,15.09,15.07,14.97,14.95,14.78,14.67,14.64,14.51,14.41,
                  14.20,13.86,13.61,13.38,12.99,12.56,12.14,11.57,11.03,10.35,9.791,9.005,
                  8.027,7.339,6.750,5.975,5.299,4.488,3.792,2.874,1.863,0.953,0.173]) #dB

print(np.shape(Gain_1))
print(np.shape(Nominal_Power))


# plt.scatter(Nominal_Power,Gain_1)
# plt.axhline(y=14.20,color='red')
# plt.axvline(x=Nominal_Power[12],color='red')
# plt.xlabel('Nominal Power [dBm]')
# plt.ylabel('Gain [dB]')
# plt.title('425 MHz')
#
# dB_comp_1 = Nominal_Power[12] # 1 dB compression point


# Added a 10 dB attenuator to port 1 of the VNA, so that I could "virtually" measure the Gain at a nominal power=-57 dBm.
# The VNA has a lower limit set at -47 dBm. The above step ensures that there is no increase in the Gain beyond -47 dBm (-48....-57)
# We need to check it, otherwise one might have wrong values for the 1 dB compression point.

# Nominal_Power_att = np.arange(-57,-47,1)
#
# Gain_att = np.array([5.15,5.15,5.11,5.05,5.05,5.05,5.03,5.00,4.956,4.885]) # very low stability and the gain should decrease by approx 10 dB
# # which it does...
# print(np.shape(Nominal_Power_att))
# print(np.shape(Gain_att))
#
# Gain_1 = Gain_1-10.16
# Gain_att = Gain_att.tolist()
# Gain_1 = Gain_1.tolist()
#
# Gain = Gain_att + Gain_1
# Gain = np.array(Gain)
#
# print(Gain)
#
# plt.scatter(Nominal_Power[0:32],Gain[0:32]) # approx final plot after 10 dB attenuation

# Phase vs Frequency. Will need to import data from VNA directly...


# See PyCharm help at https://www.jetbrains.com/help/pycharm/

# Marker at 600.03 MHz.
Nominal_Power_2 = np.arange(-47,-9,1)
Gain_2 = np.array([15.26,15.23,15.17,15.11,15.11,15.10,15.09,15.06,14.95,14.83,14.83,14.81,14.72,
                   14.67,14.56,14.42,14.38,14.12,13.99,13.66,13.13,12.84,12.37,11.55,11.21,10.36,
                   9.631,9.207,8.542,7.729,6.590,5.739,5.257,4.043,3.500,2.473,1.525,0.089]) #Much better stability

plt.scatter(Nominal_Power_2,Gain_2)
plt.axhline(y=14.26,color='red')
plt.axvline(x=Nominal_Power_2[16],color='red')
plt.xlabel('Nominal Power [dBm]')
plt.ylabel('Gain [dB]')
plt.title('600 MHz')

dB_comp_2 = Nominal_Power_2[16] # 1 dB compression point
print(dB_comp_2)
plt.show()