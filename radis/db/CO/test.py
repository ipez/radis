from radis import getMolecule
from radis.phys import cm2nm

CO_A = getMolecule("CO", 1, "A")
CO_B = getMolecule("CO", 1, "B")
tmp = CO_B.Erovib(1, 0) - CO_A.Erovib(0, 0)
tmp1 = CO_B.Erovib(1, 0) - CO_A.Erovib(1, 0)
tmp2 = CO_B.Erovib(0, 0) - CO_A.Erovib(0, 0)
tmp3 = CO_B.Erovib(0, 0) - CO_A.Erovib(1, 0)
tmp4 = CO_B.Erovib(0, 0) - CO_A.Erovib(2, 0)
tmp5 = CO_B.Erovib(1, 0) - CO_A.Erovib(4, 0)
tmp6 = CO_B.Erovib(0, 0) - CO_A.Erovib(3, 0)
tmp7 = CO_B.Erovib(0, 0) - CO_A.Erovib(4, 0)
print(cm2nm(tmp), cm2nm(tmp1), cm2nm(tmp2), cm2nm(tmp3), cm2nm(tmp4), cm2nm(tmp5), cm2nm(tmp6), cm2nm(tmp7))

from radis.levels.partfunc import PartFunc_Dunham

Qf_A = PartFunc_Dunham(CO_A)
print(Qf_A.at(T=3000))

from radis.lbl.base import linestrength_from_Einstein

A_b3_0 = (2.625e6, 4.021e6, 3.442e6, 2.159e6, 1.102e6, 4.811e5, 1.847e5, 6.312e4, 1.927e4, 5.252e3, 1.265e3, 2.626e2,
          4.389e1, 4.756e0, 8.036e-2, 1.739e-1, 3.496e-1, 3.098e-1, 2.077e-1, 1.208e-1, 6.365e-2)
A_b3_1 = (4.583e6, 9.020e5, 1.875e5, 1.557e6, 2.160e6, 1.745e6, 1.048e6, 5.145e5, 2.156e5, 7.882e4, 2.528e4, 7.077e3,
          1.699e3, 3.315e2, 4.444e1, 1.684e0, 9.006e-1, 2.884e0, 3.127e0, 2.425e0, 1.577e0)
A_B1_0 = (1.171e6, 2.321e6, 2.579e6, 2.104e6, 1.407e6, 8.190e5, 4.302e5, 2.087e5, 9.480e4, 4.063e4, 1.648e4, 6.312e3,
          2.271e3, 7.571e2, 2.286e2, 6.004e1, 1.271e1, 1.839e0, 1.099e-1, 6.119e-7)
A_B1_1 = (0.661e-2, 1.624e-1, 1.108e0, 1.400e0, 6.786e0, 2.836e1, 1.154e2, 2.667e2, 8.140e2, 1.228e3, 3.832e3, 3.166e3,
          1.209e4, 4.106e3, 2.777e4, 1.842e3, 4.654e4, 5.818e1, 5.269e4, 2.720e3, 3.259e4)

Qf_B = PartFunc_Dunham(CO_B)

CO_b3 = getMolecule("CO", 1, "b_")
CO_a3 = getMolecule("CO", 1, "a_")

Qf_b_ = PartFunc_Dunham(CO_b3)
Qf_a_ = PartFunc_Dunham(CO_a3)

from radis.api.hdf5 import DataFileManager

manager = DataFileManager("pytables")
dfA = manager.read("molecules_data_CO_iso1_A_levels.h5")
dfB = manager.read("molecules_data_CO_iso1_B_levels.h5")
dfa_ = manager.read("molecules_data_CO_iso1_a__levels.h5")
dfb_ = manager.read("molecules_data_CO_iso1_b__levels.h5")

E_values_A = []
df = dfA
# 遍历 v_value 从 0 到 20
for v_value in range(20):  # 因为 range 是前闭后开，所以需要到 21
    # 假设 j_value 始终为 0，如果需要变化，可以类似地遍历 j_value
    j_value = 0
    E_value = df.loc[(df['v'] == v_value) & (df['j'] == j_value), 'E'].values[0]
    E_values_A.append(E_value)

E_values_B = []
df = dfB
# 遍历 v_value 从 0 到 2
for v_value in range(2):  # 因为 range 是前闭后开，所以需要到 2
    # 假设 j_value 始终为 0，如果需要变化，可以类似地遍历 j_value
    j_value = 0
    E_value = df.loc[(df['v'] == v_value) & (df['j'] == j_value), 'E'].values[0]
    E_values_B.append(E_value)

E_values_a_ = []
df = dfa_
# 遍历 v_value 从 0 到 2
for v_value in range(19):  # 因为 range 是前闭后开，所以需要到 2
    # 假设 j_value 始终为 0，如果需要变化，可以类似地遍历 j_value
    j_value = 0
    E_value = df.loc[(df['v'] == v_value) & (df['j'] == j_value), 'E'].values[0]
    E_values_a_.append(E_value)

E_values_b_ = []
df = dfb_
# 遍历 v_value 从 0 到 2
for v_value in range(2):  # 因为 range 是前闭后开，所以需要到 2
    # 假设 j_value 始终为 0，如果需要变化，可以类似地遍历 j_value
    j_value = 0
    E_value = df.loc[(df['v'] == v_value) & (df['j'] == j_value), 'E'].values[0]
    E_values_b_.append(E_value)

wn = []
I = []

from itertools import islice

gu = 1
Ia = .986544E+00
T = 1000
Q = Qf_A.at(T) + Qf_B.at(T)
Eu = E_values_B[0]
lines = []
for A, El in islice(zip(A_B1_0, E_values_A), 10):
    line = linestrength_from_Einstein(A, gu, El, Ia, Eu - El, Q, T)
    lines.append(line)
    wn.append(Eu - El)
    I.append(line)
print(lines)

Eu = E_values_B[1]
lines_1 = []
for A, El in islice(zip(A_B1_1, E_values_A), 10):
    line = linestrength_from_Einstein(A, gu, El, Ia, Eu - El, Q, T)
    lines_1.append(line)
    wn.append(Eu - El)
    I.append(line)
print(lines_1)

Q = Qf_a_.at(T) + Qf_b_.at(T)
Eu = E_values_b_[0]
lines_b = []
for A, El in islice(zip(A_b3_0, E_values_a_), 10):
    line = linestrength_from_Einstein(A, gu, El, Ia, Eu - El, Q, T)
    lines_b.append(line)
    wn.append(Eu - El)
    I.append(line)
print(lines_b)

Eu = E_values_b_[1]
lines_b_1 = []
for A, El in islice(zip(A_b3_1, E_values_a_), 10):
    line = linestrength_from_Einstein(A, gu, El, Ia, Eu - El, Q, T)
    lines_b_1.append(line)
    wn.append(Eu - El)
    I.append(line)
print(lines_b_1)

import numpy as np
import astropy.units as u

wl = 1e7 / np.array(wn)
from radis import Spectrum

s = Spectrum.from_array(wl, np.array(I), 'radiance_noslit', wunit='nm', Iunit="cm-1")

s.sort()
s.normalize()

filename = "radiance_noslit_"+str(T)+".txt"
s.savetxt(filename, 'radiance_noslit', wunit='nm', Iunit='cm-1')

w_nm = s.get_wavelength()
print(w_nm)
# s.resample_even()
# s.resample_even(energy_threshold=1e10)
s.apply_slit(3, 'nm', shape="gaussian")

s.plot('radiance')
