from sympy import *
import numpy as np
init_printing(use_unicode=True)
T,t,T0,dT=symbols('T,t,T0,dT')

eq=-t/log(1-((T-T0)/dT))

T0der=diff(eq,T0)
tder=diff(eq,t)
dTder=diff(eq,dT)
Tder=diff(eq,T)

a=T0der.subs([(t,2),(T,51.5),(T0,50),(dT,2)])**2
b=tder.subs([(t,2),(T,51.5),(T0,50),(dT,2)])**2
c=dTder.subs([(t,2),(T,51.5),(T0,50),(dT,2)])**2
d=Tder.subs([(t,2),(T,51.5),(T0,50),(dT,2)])**2

x=0.08**2#sıcaklık
y=(2.5*10**-9)**2#zaman

uncert=sqrt(a*x+b*y+c*x+d*x)
print(uncert)
