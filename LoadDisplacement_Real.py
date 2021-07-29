import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import scipy.interpolate as interp
from scipy.interpolate import splev,splrep
from scipy.signal import savgol_filter
from scipy.fftpack import rfft, irfft, rfftfreq
from nptdms import TdmsFile
import numpy as np
import glob
import os
from scipy.signal import medfilt

"""
This code is for reading experimental data of Indenter/RC Thermal setup obtained by
NI-DAQ written/converted txt data files are being searched as txt. 
"""
## Calibration constants of Experiment ##
lc_mass_cal=4.75
grave=9.81
bit_to_load=lc_mass_cal*grave
loadcell_stiffness=536.5 ##mN/um
gage2um=10 ##um/V
###

load_files=[]
position_files=[]
temp_files=[]
sg_files=[]
lvsdis_file=[[]]

date=[[]]
n=0
a=[]
b=[]
c=[]
#t1=np.linspace(1,1000,1000)
#t2=np.linspace(999,1999,1000)

## Fit functions for experiment data
def rise(t,N1,tau1,c):
    """The model function for the rise of the curve"""
    return N1*(1-np.exp(-t/tau1))+c

def decay(t,N2,tau2,c):
    """Model function for the decay part of the curve"""
    return N2*(np.exp(-t/tau2))+c

def order4(x,a,b,c,d,e):
    return a*x**4+b*x**3+c*x**2+d*x+e

def linear(x,m,c):
    return m*x+c

def quad(x,a,b,c):
    return a*x**2+b*x+c
###

##Special Characters##
degsin=u"\N{DEGREE SIGN}"
degcel=degsin+"C"
deltasig="\u0394"
mikron="\u03BC"
###

## Main file of Expriment Data of text files -that converted from NI TDMS-
main="/Users/ahmeduluca/Desktop/Cu/Cu-Process/22-06-2021"
###

file_list=[f.path for f in os.scandir(main) if f.is_dir()]
file_list.sort()
expdates=[f.name for f in os.scandir(main) if f.is_dir()]
expdates.sort()
m=0
tcount=0
ch=0
starter=-1
final=-1
tpt0=[]
tpt1=[]
transit=0
t1=[]
t2=[]
tempo=[]
sg=[]
totaldisp=[]
totaltemp=[]
totalload=[]
totalpos=[]
## Search for Experiment Files that include Step Files ##
for i in file_list:
    load=[[]]
    position=[[]]
    voltage=[[]]
    temp=[[[],[]]]
    time=[[]]
    totaldisp=[]
    totaltemp=[]
    totalload=[]
    sub_dirs=[f.path for f in os.scandir(i) if f.is_dir()]
    sub_dirs.sort()
    stepnames=[f.name for f in os.scandir(i) if f.is_dir()]
    stepnames.sort()
    n=0
    
    ## Search for txt files of Load - Strain Gage - Temperatures ##
    for j in sub_dirs:
        tcount=0
        txt_files=glob.glob(os.path.join(j+'/*.txt'))
        ch=0
        for k in txt_files:
            if("Load" in k):
                load.append([])
                f=open(k,"r")
                for x in f:
                    load[n].append(float(x)*bit_to_load/1000)
            elif("Temperature" in k):
                if(ch==0 and n>0):
                    temp.append([[],[]])
                if(tcount==0):
                    time.append([])
                with open(k,'r') as f:
                    rows=f.readlines()
                    for x in rows:
                        if(tcount==0):
                            time[n].append(float(x.split()[0]))
                        temp[n][ch].append(float(x.split()[1]))
                tcount+=1
                ch+=1                    
            elif("Voltage" in k):
                voltage.append([])
                if(tcount==0):
                    time.append([])
                with open(k,'r') as f:
                    rows=f.readlines()
                    for x in rows:
                        if(tcount==0):
                            time[n].append(float(x.split()[0]))
                        voltage[n].append(float(x.split()[1])*gage2um)##um conversion
                tcount+=1
        ###
                
        ## Load & Displacement & Temperatures vs Time of each Step separately ##
        load[n]=medfilt(load[n],21)
        voltage[n]=medfilt(voltage[n],9)
        temp[n][1]=medfilt(temp[n][1],3)
        if(n==0):
            min_vol=min(voltage[0])
        for t in range(len(voltage[n])):
            voltage[n][t]=voltage[n][t]-min_vol
        orta=int(len(voltage[n])/len(load[n]))
        p2=[]
        for t in range(0,len(voltage[n]),orta):
            p2.append(np.mean(voltage[n][t:t+orta]))
        leng=len(p2)-1
        artik=leng%len(load[n])
        if(artik==0):
            artik=1
        out=int(leng/artik)
        bas=leng%artik
        for t in range(bas,leng,out):
            position[n]=np.concatenate((np.array(position[n]),p2[t:t+out-1]))
        inc=0
        gg=0
        for t in range(len(position[n])):
            position[n][t]=position[n][t]-load[n][t]/loadcell_stiffness
##        fig,ax1=plt.subplots()
##        fig.subplots_adjust(right=0.75)
##        ax2=ax1.twinx()
##        p1,=ax1.plot(time[n],voltage[n],'g.',zorder=4,label="Actuator Position %s = %.2f %sm"%(deltasig,round(max(voltage[n])-min(voltage[n]),2),mikron))
##        p2,=ax2.plot(np.linspace(0,time[n][-1],len(load[n])),load[n],'b.',zorder=5,label="Load %s = %.2f mN"%(deltasig,round(max(load[n])-min(load[n]),2)))
##        #ax1.plot(np.linspace(0,time[n][-1],len(position[n])),position[n],label="Pst")
##        ax1.set_xlabel('Time (s)')
##        ax1.set_ylabel('Position (%sm)'%(mikron), color='g')
##        ax2.set_ylabel('Load (mN)', color='b')
##        if("Indentation" in stepnames[n]):
##            ax3=ax1.twinx()
##            ax3.spines['right'].set_position(("axes",1.2))
##            p3,=ax3.plot(time[n],temp[n][1],'r.',ms=3,label="Tip Temperature %s = %.2f %s"%(deltasig,round(max(temp[n][1])-min(temp[n][1])),degcel))
##            ax3.set_ylabel('Temperature '+degcel, color='r')
##            leg=ax1.legend(handles=[p1,p2,p3], bbox_to_anchor=(0., 1.05, 1.24, .104), loc='lower left',
##               ncol=2, mode="expand", borderaxespad=0.)
##            leg.set_draggable('True')
##        else:
##            leg=ax1.legend(handles=[p1,p2], bbox_to_anchor=(0., 1.05, 1., .102), loc='lower left',
##               ncol=2, mode="expand", borderaxespad=0.)
##            leg.set_draggable('True')
##        ax2.grid("True",color='blue',which='both')
##        ax1.grid("True",color="black",axis='x')
##        plt.title(stepnames[n])
##        plt.show()
        ###

        ## Load vs Displacement Graph for each of Indent & Oscillation and corresponding Fits ##
        fark=len(load[n])-len(position[n])
        if(fark>0):
            load[n]=load[n][:-fark]
            plt.plot(position[n],load[n],'.')
        elif(fark==0):
            plt.plot(position[n],load[n],'.')
        else:
            position[n]=position[n][:fark]
            plt.plot(position[n],load[n],'.')
        loadav=np.mean(load[n])
        loadmin=min(load[n])
        loadmax=max(load[n])
        dload=loadmax-loadmin
        fitst=0
        fitend=0
        #if('Indentation' in stepnames[n]):
##            if(n+1<len(stepnames) or n<1):
##                for t in load[n]:
##                    fitst+=1
##                    if(t>loadmin+dload*0.025):
##                        break
##                for t in load[n]:
##                    fitend+=1
##                    if(t>loadmin+dload*0.9):
##                        break
##                curvest=0
##                for t in load[n]:
##                    curvest+=1
##                    if(t>loadmin+dload*0.1):
##                        break
##                try:## Quadratic fit to load vs disp -> F=h**2
##                    linpov,lincov=curve_fit(quad, position[n][fitst:fitend],load[n][fitst:fitend])
##                    linfit=quad(position[n][curvest:fitend],*linpov)
##                    plt.plot(position[n][curvest:fitend],linfit,label="(%.2f)$h^2$+(%.2f)h+(%.2f)"%(round(linpov[0],2),round(linpov[1],2),round(linpov[2],2)))
##                except:
##                    print("quadratic fit problem!")
##            if(n>0 and n+1==len(stepnames)):
##                for t in load[n]:
##                    fitend+=1
##                    if(t<loadmin+dload*0.75):
##                        break
##                try:## Linear fit to first 10% unloading load vs disp to obtain h_contact
##                    fitst=int(fitend*0.1)
##                    linpov,lincov=curve_fit(linear, position[n][fitst:fitend],load[n][fitst:fitend])
##                    linstop=-1
##                    for t in position[n][fitst:]:
##                        linstop+=1
##                        if (t<=(loadmin-linpov[1])/linpov[0]):
##                            break
##                    linfit=linear(position[n][fitst:linstop],*linpov)
##                    plt.plot(position[n][fitst:linstop],linfit,label="$h_{c}$=%f"%((loadmin-linpov[1])/linpov[0]))
##                except:
##                    print("Linear Fit Problem!")
##        plt.xlabel("Displacement (%sm)"%(mikron))
##        plt.ylabel("Load (mN)")
        plt.grid('True')
##        plt.title(stepnames[n]+" Load vs Displacement")
##        plt.legend()
##        plt.show()
###
        
        ## Total Displacement - Load - Temperature vs Time through all steps Graph of Experiment ## !!LABELS for separating steps
        if("Indentation" in stepnames[n]):
            for t in range(len(voltage[n])):
                totaldisp.append(voltage[n][t])
                totaltemp.append(temp[n][1][t])
            for t in range(len(load[n])):
                totalload.append(load[n][t])
            for t in range(len(position[n])):
                totalpos.append(position[n][t])
        n+=1
        position.append([])
##    fig,ax1=plt.subplots()
##    fig,ax1=plt.subplots()
##    fig.subplots_adjust(right=0.75)
##    ax2=ax1.twinx()
##    ax3=ax1.twinx()
##    ax3.spines['right'].set_position(("axes",1.2))
##    p3,=ax3.plot(0.005*np.arange(0,len(totaltemp)),totaltemp, 'r.',ms=3,label="Tip Temperature %s = %.2f %s"%(deltasig,round(max(totaltemp)-min(totaltemp),2),degcel))
##    p1,=ax1.plot(0.005*np.arange(0,len(totaldisp)),totaldisp,'g.',zorder=3,label="Actuator Position %s = %.2f %sm"%(deltasig,round(max(totaldisp)-min(totaldisp),2),mikron))
##    p2,=ax2.plot(0.1*np.arange(0,len(totalload)),totalload, 'b.',zorder=2,label="Load %s = %.2f mN"%(deltasig,round(max(totalload)-min(totalload),2)))
##    plt.title("Experiment: "+expdates[m])
##    ax1.set_xlabel('Time (s)')
##    ax1.set_ylabel('Position (%sm)'%(mikron), color='g')
##    ax2.set_ylabel('Load (mN)', color='b')
##    ax3.set_ylabel('Temperature '+degcel, color='r')
##    ax1.grid("True",color='green')
##    ax1.legend(handles=[p1,p2,p3])
##    plt.show()
    ### ++ Total Load vs Disp graph of Experiment through all steps!!
    plt.title(expdates[m]+" Load Displacement")
    plt.xlabel("Displacement ("+mikron+"m)")
    plt.ylabel("Load (mN)")
    plt.show()
    m+=1

