# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 11:48:40 2019

@author: cemsa
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import medfilt
from scipy.optimize import curve_fit

a=[]
b=[]           
c=[]
Voltage=[]
Tip=[]
Sample=[]


DAR=1000 #[Hz] Data Acqusition Rate
EF=0.500 #[Hz] Experiment Frequency
Sampling=(1/EF)*DAR
filename='24-04-2019-14.txt' #Put the filename with correct extension



#######Reading Data###########



A=filename.find('txt')#.find returns the starting point
#of the string you look for,if there is not it returns -1

if A > 0:

#USE IF IT IS A TXT
    with open(filename,'r') as f:
        rows=f.readlines() #Exclude the headers with [1:]

        for x in rows:
            Tip.append(x.split()[2])
            Sample.append(x.split()[3])
            Voltage.append(x.split()[1])

else:
#USE IF IT IS EXCEL
        
    df=pd.read_excel(filename)#[starting line:last line]

    S=df.shape
    Voltage=df['Voltage0']
    Tip=df['Temperature15']
    Sample=df['Temperature12']

        
###############################
    

for j in range(0,int(Sampling)): #Sum the same point in every cycle
    x=[]
    y=[]
    z=[]
    
    for i in range(j,len(Tip),int(Sampling)):
        
        x.append(float(Tip[i]))
        y.append(float(Sample[i]))
        z.append(float(Voltage[i]))
            
    a.append(sum(x)/len(x))#Tip
    b.append(sum(y)/len(y))#Sample
    c.append(sum(z)/len(z))#Voltage


        

ref1=0
ref2=0
ref3=0
count=0
tbma=[]#data To Be Moved
tbmb=[]
tbmc=[]
tol=(abs(abs(max(c))-abs(min(c))))/2


#####Determine where the Voltage level changes character

for u in range(len(c)):
     
    ref1=ref1+1
     
    if min(c)-tol < c[u] < min(c)+tol:
         ref3=ref3+1
         if ref1 != ref3:
             pindown=u
             ref3=ref1
  
    else:
         ref2=ref2+1
         if ref2 != ref1:
             pinup=u
             ref2=ref1


print('pindown=',pindown,'pinup=',pinup)

#############Extract the Data to be moved################
for p in range(len(c)):
    
    if p>pindown:
             tbma.append(a[p])
             tbmb.append(b[p])
             tbmc.append(c[p])
    
    
#######Patch data to the beginning of the list###########


tbma.extend(a[:pindown])
tbmb.extend(b[:pindown]) 
tbmc.extend(c[:pindown])       


Resultname=str(filename)
with open('sonuc.txt', "w") as n:
    for k in range(len(tbma)):
        n.write(str(tbma[k]))
        n.write("    ")
        n.write(str(tbmb[k]))
        n.write("    ")
        n.write(str(tbmc[k]))
        n.write("\n")

ared=medfilt(tbma,kernel_size=101)#a (Tip) data with reduced noise

plt.figure(figsize=(14,10))

#####FITTING THE CURVES########

t1=np.linspace(1,1000,1000)
t2=np.linspace(999,1999,1000)
aredrise=ared[:1000]
areddecay=ared[999:]

def f1(t,N1,tau1,c):
    """The model function for the rise of the curve"""
    return N1*(1-np.exp(-t/tau1))+c


def f2(t,N2,tau2,c):
    """Model function for the decay part of the curve"""
    return N2*(np.exp(-t/tau2))+c

popt1, pcov1 = params1, params_covariance1 =curve_fit(f1,t1,aredrise,maxfev=100000,p0=[40,30,40])
print('popt1=')
print(popt1)
print('pcov1=')
print(pcov1)


popt2, pcov2 = params2, params_covariance2 =curve_fit(f2,t2,areddecay,maxfev=100000,p0=[-55,63,47])
print('popt2=')
print(popt2)
print('pcov2=')
print(pcov2)
    
Risefit=f1(t1, *popt1)
Decayfit=f2(t2,*popt2)

with open('time_constants.txt',"a") as q:
     q.write(str(filename))
     q.write("    ")
     q.write(str(popt1[1]))
     q.write("    ")
     q.write(str(popt2[1]))
     q.write("\n")





#ax1=plt.subplot(3,1,1)
plt.plot(tbma,'orange') #tip temp
plt.hold('on')
plt.plot(ared,'red')
plt.plot(t1,Risefit,'black')
plt.plot(t2,Decayfit,'black')
#ax1.set_title('Tip')
#ax2=plt.subplot(3,1,2)
#plt.plot(tbmb,'blue')  #sample temp
#ax2.set_title('Sample')
#ax3=plt.subplot(3,1,3)
#plt.plot(tbmc,'green') #voltaj
#ax3.set_title('Voltage')


            




    
