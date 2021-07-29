import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import scipy.interpolate as interp
from scipy.signal import savgol_filter
from scipy.fftpack import rfft, irfft, rfftfreq
from nptdms import TdmsFile
import numpy as np
import glob
import os
from scipy.signal import medfilt


def f1(t,N1,tau1,c):
    """The model function for the rise of the curve"""
    return N1*(1-np.exp(-t/tau1))+c


def f2(t,N2,tau2,c):
    """Model function for the decay part of the curve"""
    return N2*(np.exp(-t/tau2))+c

files=glob.glob("/Users/ahmeduluca/Desktop/Cu/Cu-Process/22-06-2021/**/**/*.txt")
filename=""
ared=[]
taur=[]
taud=[]
time=[]
files.sort()
for i in files:
    if("Indentation" in i and "Temperature2" in i):

        ared=[]
        time=[]
        filename=os.path.dirname(i)
        with open(i,'r')as a:
            rows=a.readlines()
            for f in rows:
                ared.append(float(f.split()[1]))
                time.append(float(f.split()[0]))
#plt.plot(ared)
#plt.show()
        ared=medfilt(ared,kernel_size=21)#a (Tip) data with reduced noise
        mean=np.mean(ared)
        start=0
        for i in range(len(ared)):
            if(ared[i]<mean-0.6 and ared[-1]>ared[0]):
                start=i
                break
            elif(ared[i]>mean+0.3 and ared[-1]<ared[0]):
                start=i
                break
        print(start)
        print("mean = "+str(mean))
        for i in range(start,len(ared)):
            if(ared[i]>mean and ared[-1]>ared[0]):
                ared=ared[i:]
                time=time[i:]
                break
            elif(ared[i]<mean-0.3 and ared[-1]<ared[0]):
                ared=ared[i:]
                time=time[i:]
                break
#####FITTING THE CURVES########
        t1=np.linspace(0,2105,421)
        t2=np.linspace(2110,4215,421)
        t3=np.concatenate((t1,t2))
        print(i)
#        if("Step 5 Oscillation" in i):
#            ared=np.concatenate((ared[530:],ared[0:530]))
        areddecay=ared[0:421]
        aredrise=ared[421:]
        amax=max(ared)
        amin=min(ared)
        delta=amax-amin
        try:
            popt1, pcov1 = params1, params_covariance1 =curve_fit(f1,time,ared,maxfev=100000,p0=[delta,500,amin])
            print('popt1=')
            print(popt1)
            print('pcov1=')
            print(pcov1)
            Risefit=f1(np.array(time), *popt1)
            taur.append(popt1[1])
            with open(filename+"/timeconst.txt","w") as q:
                q.write(str(filename)+"\n"+"Rise TC \t Decay TC\n"+str(popt1[1])+"\n")
            plt.plot(time,ared,'red')
            plt.plot(time,Risefit,'black',label="TC ="+str(round(popt1[1]))+"s")
#          plt.plot(t1,Decayfit,'yellow',label="T_decay ="+str(round(popt2[1]))+"ms")
            plt.xlabel("Time (s)")
            plt.ylabel("Temperature (C)")
            plt.legend()
            plt.show()
        except:
            popt2, pcov2 = params2, params_covariance2 =curve_fit(f2,time,ared,maxfev=100000,p0=[-delta,500,amax])
            print('popt2=')
            print(popt2)
            print('pcov2=')
            print(pcov2)
            Decayfit=f2(time,*popt2)
            taud.append(popt2[1])
            with open(filename+"/timeconst.txt","w") as q:
                q.write(str(filename)+"\n"+"Decay TC\n"+"\t"+str(popt2[1])+"\n")
#ax1=plt.subplot(3,1,1)
#plt.plot(tbma,'orange') #tip temp
            plt.plot(time,ared,'red')
            plt.plot(time,Decayfit,'yellow',label="T_decay ="+str(round(popt2[1]))+"s")
            plt.xlabel("Time (s)")
            plt.ylabel("Temperature (C)")
            plt.legend()
            plt.show()
with open(filename+"/Alltimeconst.txt","w") as q:
    q.write("Rise\n")
    for i in range(len(taur)):
        q.write(str(taur[i])+"\n")
    q.write("Decay\n")
    for i in range(len(taud)):
        q.write(str(taud[i])+"\n")      
#plt.plot(taud,'r.',label="decay TC")
plt.plot(taur,'b.', label="Indent TCs")
plt.ylabel("tau (s)")
plt.xlabel("Step")
plt.legend()
plt.show()
#ax1.set_title('Tip')
#ax2=plt.subplot(3,1,2)
#plt.plot(tbmb,'blue')  #sample temp
#ax2.set_title('Sample')
#ax3=plt.subplot(3,1,3)
#plt.plot(tbmc,'green') #voltaj
#ax3.set_title('Voltage')
