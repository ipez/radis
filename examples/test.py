from radis import getMolecule
from radis.phys import cm2nm

CO_A = getMolecule("CO", 1, "A")
CO_B = getMolecule("CO", 1, "B")
tmp = CO_B.Erovib(1, 0)-CO_A.Erovib(0, 0)
tmp1 = CO_B.Erovib(1, 0)-CO_A.Erovib(1, 0)
tmp2 = CO_B.Erovib(0, 0)-CO_A.Erovib(0, 0)
tmp3 = CO_B.Erovib(0, 0)-CO_A.Erovib(1, 0)
tmp4 = CO_B.Erovib(0, 0)-CO_A.Erovib(2, 0)
tmp5 = CO_B.Erovib(1, 0)-CO_A.Erovib(4, 0)
tmp6 = CO_B.Erovib(0, 0)-CO_A.Erovib(3, 0)
tmp7 = CO_B.Erovib(0, 0)-CO_A.Erovib(4, 0)
print(cm2nm(tmp),cm2nm(tmp1),cm2nm(tmp2),cm2nm(tmp3),cm2nm(tmp4),cm2nm(tmp5),cm2nm(tmp6),cm2nm(tmp7))
