"""
Coded by: Takuro TOKUNAGA, u1053018
Python version: 3.5.0
Last modified: 2017, April, 19
Separation of variables
"""
import os
import math
import cmath
import numpy

# parametres
# summation
n = 1
nmax = 100
imax = 100

# temperature
t1 = 0
t2 = 50
t3 = 100
tsteady = [0 for i in range(imax)]
t1steady = [0 for i in range(imax)]
t2steady = [0 for i in range(imax)]
t3steady = [0 for i in range(imax)]

# eigenvalue
ln1 = 0
ln2 = 0

# dimensions
w = 50*pow(10,-3)
l = 100*pow(10,-3)

# coordination
x =  0.092# change here
ymax = w
y = 0
dy = w/imax

# coefficient
c3 = (2*t3)/math.pi
c2 = (2*t2)/math.pi

# definition of function
def start():
	print("start")
	
def finish():
	print("finish")
				
# Main part
# Start of calculation
start()

for i in range(0, imax):
	for n in range(1, nmax):
		ln1 = (n*math.pi/l) # eigenvalue radian
		ln2 = (n*math.pi/w) # eigenvalue radian
		t1steady[i] = t1steady[i] + c3*((pow((-1),n+1)+1)/n)*math.sin(ln1*x)*math.sinh(ln1*y)/math.sinh(ln1*w)
		t2steady[i] = t2steady[i] + c3*((pow((-1),n+1)+1)/n)*math.sin(ln2*y)*math.sinh(ln2*x)/math.sinh(ln2*l)
		t3steady[i] = t3steady[i] + c2*(-math.cosh(ln2*l)/math.sinh(ln2*l))*((pow((-1),n+1)+1)/n)*math.sin(ln2*y)*\
					        ((-math.sinh(ln2*l)/math.cosh(ln2*l))*math.cosh(ln2*x)+math.sinh(ln2*x))
					      
	tsteady[i] = t1steady[i] + t2steady[i] + t3steady[i]
	y = y + dy

y = 0
# output for folder
f = open("./excel/analytical.txt","w")
f.write("y temp\n")
for i in range(0, imax):
	f.write(str(y))
	f.write(str(" "))
	f.write(str(tsteady[i]))
	f.write("\n")
	y = y + dy
f.close()

# End of calculation
finish()