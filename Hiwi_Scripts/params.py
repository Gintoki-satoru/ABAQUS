a = 200
b = 200

'''# Isotropic material (Example 1)
E  = 68200.0     # MPa
nu = 0.326
h  = 1.27       # mm

Q = E / (1 - nu**2)
Q11 = Q
Q22 = Q
Q12 = nu * Q
Q66 = E / (2*(1+nu))

A11 = Q11 * h
A22 = Q22 * h
A12 = Q12 * h
A66 = Q66 * h
A16 = 0.0
A26 = 0.0

D11 = Q11 * h**3 / 12.0
D22 = Q22 * h**3 / 12.0
D12 = Q12 * h**3 / 12.0
D66 = Q66 * h**3 / 12.0
D16 = 0.0
D26 = 0.0'''

A11 = 1.921389451933432e+05
A22 = 1.921389451933432e+05
A12 = 5.731224511564325e+03
A66 = 14340
A16 = 0.0
A26 = 0.0

D11 = 1.2009e+06
D22 = 1.2009e+06
D12 = 3.5820e+04
D66 = 89625
D16 = 0.0
D26 = 0.0

K11 = 3.4488e+03
K22 = 3.4488e+03

meshSize = 4
lambdaEnd = 4.0
eval1LPF=2.09
eval2LPF=3.14
dir='D:\\abaqus_temp' 
evalAtQuater = True 
