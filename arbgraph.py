import matplotlib.pyplot as plt
import numpy as np
import scipy.special as spefun
import scipy.interpolate as interp


def fourier(x,T0,Tf,l):
    d1=max(x)-min(x)
    dx=d1/len(x)
    return dx*(T0-Tf)/d1
def r_eq(k1,k2,l1,l2,r_con,A):
    d1=max(l1)-min(l1)
    d2=max(l2)-min(l2)
    return (d1/(k1*A))+(r_con)+(d2/(k2*A))
def t_contact(T0,Tf,r_eq,l1,k1,A):
    if(len(l1)>1):
        d1=max(l1)-min(l1)
    else:
        d1=l1[0]
    return (T0-Tf)*d1/(k1*A*r_eq)
k1=2.500
k2=1.500
T1=100.
T2=150.
dist1=np.linspace(0,2.,101)
dist2=np.linspace(2,4.,101)
area=1.
rcon=1.
Tcon1=0.
Tcon2=0.
tot_res=r_eq(k1,k2,dist1,dist2,rcon,area)
dtcon1=t_contact(T1,T2,tot_res,dist1,k1,area)
Tcon1= (T1-dtcon1)
dtcon2=t_contact(T1,T2,tot_res,[rcon],1,1)
Tcon2= (Tcon1-dtcon2)
print(Tcon1)
print(Tcon2)
print(tot_res)
Temp1=[T1]
Temp2=[Tcon2]
grad1=fourier(dist1,T1,Tcon1,dist1)
grad2=fourier(dist2,Tcon2,T2,dist2)
print(grad1)
print(grad2)
if(np.round(grad1*k1,decimals=7)==np.round(grad2*k2,decimals=7)):
    print('check, q=%.7f'%(grad1*k1))
else:
    print('q1=%f,q2=%.7f' % (grad1*k1,grad2*k2))
for i in range(1,len(dist1)):
    Temp1.append(Temp1[i-1]-grad1)
for j in range(1,len(dist2)):
    Temp2.append(Temp2[j-1]-grad2)
print(len(dist1))
print(len(dist2))
plt.plot(dist1,Temp1,'r')
plt.plot(dist2,Temp2,'b')
distance=np.concatenate((dist1,dist2))
plt.plot(distance,Temp1+Temp2,'--')
plt.xlabel('Distance')
plt.ylabel('Temperature (K)')
plt.title('Temperature Graph of Contact')
plt.show()

f1=interp.interp1d(dist1, Temp1,fill_value='extrapolate')
f2=interp.interp1d(dist2, Temp2,fill_value='extrapolate')
x1=np.arange(min(dist1),max(dist1)+0.01,max(dist1)-min(dist1))
x2=np.arange(max(dist1),max(dist2)+0.011,max(dist2)-min(dist2))
print(np.mean(dist2))
y1=f1(x1)
y2=f2(x2)
y=np.concatenate((y1,y2))
print(y)
x=np.concatenate((x1,x2))
print(x)
x_ticks=['X1', 'Contact Pt.', 'Contact Pt.','X2']
y_ticks=['T1','T_Sur1','T_Sur2','T2']
plt.ylabel('Temperature (K)')
plt.xticks(x, x_ticks)
plt.yticks(y, y_ticks)
plt.xlabel('Distance')
plt.title('Temperature Graph of Contact')
plt.plot(x,y,'o-r')
plt.annotate('dT/dx = Q*R1', xy=(np.mean(dist1), np.mean(Temp1)), xytext=(0.1+np.mean(dist1)/2,1+abs(Tcon1/4)),textcoords='offset points',
             arrowprops=dict(arrowstyle='-'),# shrink=0.05),
             )
plt.annotate('dT/dx = Q*R2',   xy=(np.mean(dist2), np.mean(Temp2)), xytext=(0.1+np.mean(dist2)/2,1+abs(Tcon2/4)),textcoords='offset points',
             arrowprops=dict(arrowstyle='-'),
             )
plt.annotate('', xy=(max(dist1), Tcon1), xytext=(max(dist1), Tcon2),
             arrowprops=dict(arrowstyle="]-[",connectionstyle="bar,fraction=0.1"), #shrink=0.05),
             )
plt.text(s='deltaT/Q \n= R_contact/A',x=(np.mean(dist1)),y=(np.mean([Tcon1,Tcon2]))
             )
plt.grid(True)
plt.show()
##deltaT/Q = R_contact/A
