import pandas as pd
import numpy as np
import glob
import os
from scipy.special import iv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from io import open
import math
from numpy import inf
mainfile="/Volumes/AhmedUluca/Indenter_Data/main"

tipfiles=[f.path for f in os.scandir(mainfile) if f.is_dir()]

def trunConeConst(d,C,A,k,R_0=10,alpha=60):
    alpha=math.radians(alpha)
    return (C*A/(2*k))*(10**6/(np.pi*R_0**2+(d*R_0/np.cos(alpha))+d**2*np.tan(alpha)/np.cos(alpha)))
def generalArea(d,C,beta,h_c):
    return C/(d**2*h_c*beta*10**-12)
def trunConeArea(d,C,h_c,alpha=60,R_0=10):
    alpha=math.radians(alpha)
    return (C/h_c)*(10**12/(np.pi*R_0**2+(d*R_0/np.cos(alpha))+d**2*np.tan(alpha)/np.cos(alpha)))
def conic_etha(h, k, d, l):
    m=2*sqrt(h/(k*d))
    return (2/(m*l))*(iv(2,2*m*l)/iv(1,2*m*l))
def Spreading(contact_dia,k1,k2):
    return contact_dia*(2*k1*k2/(k1+k2))
def totalFun(d,C,A,B,m,n,d0,t0=2000):
    return t0+A*d+B*d**-1+C*d**-2
seqName="Contact"
processName="Indent Sweep 1"
expresults=[]
for tip in tipfiles:
    samplefiles=[f.path for f in os.scandir(tip) if f.is_dir()]
    for sample in samplefiles:
        procfiles=[f.path for f in os.scandir(sample) if f.is_dir()]
        for procedure in procfiles:
            expfiles=[f.path for f in os.scandir(procedure) if f.is_dir()]
            for experiment in expfiles:
                txts=glob.glob(experiment+"/*.txt")
                for text in txts:
                    if("riseTC_disp_load" in text):
                        with open(text,"r",encoding="utf-8") as f:
                            lines = f.readlines()[1:]
                            for line in lines:
                                expresults.append({"Tip":os.path.basename(tip),
"Sample":os.path.basename(sample),"Procedure":os.path.basename(procedure),
"Experiment":os.path.basename(experiment).replace("_"," "),"Sequence":"Contact","Depth":float(line.split()[1]),
"Load":float(line.split()[2]),"Time Constant":float(line.split()[0])})
                    elif("decayTC_disp_load" in text):
                        with open(text,"r",encoding="utf-8") as f:
                            lines = f.readlines()[1:]
                            for line in lines:
                                expresults.append({"Tip":os.path.basename(tip),
"Sample":os.path.basename(sample),"Procedure":os.path.basename(procedure),
"Experiment":os.path.basename(experiment).replace("_"," "),"Sequence":"Separation","Depth":float(line.split()[1]),
"Load":float(line.split()[2]),"Time Constant":float(line.split()[0])})
data=pd.DataFrame(expresults)
#data.to_excel(os.path.join(mainfile,"AllResults.xlsx"),engine='xlsxwriter')
#data.groupby("Sequence").plot.hexbin(x='Load',y='Time Constant',C='Depth',subplots='True')
#data.groupby("Sequence").plot.hexbin(x='Depth',y='Time Constant',C='Load',subplots='True')
#data.plot(x='Load',y='Time Constant')
#ax.clabel("Depth")
#data.groupby("Sample").plot.hexbin(x='Depth',y='Time Constant',C='Load',subplots='True')
tidy=data.groupby(["Experiment"])
##for name,group in tidy:
##    minval=min(group.Depth)
##    for i,row in data[data.Experiment==name].iterrows():
##        data["Depth"][i]-=minval
data.sort_values(by="Depth")
gurup=data.groupby(["Tip","Sample","Sequence","Procedure","Experiment"])
fig,axs=plt.subplots(2,3)
date0=""
for name,group in gurup:
    if(seqName in name and "Conic" in name and processName in name[3]):
        if("Al" in name[1]):
            axs[0,0].plot(group.Depth,group["Time Constant"],marker='o',ls=' ')
            axs[0,0].set_title(name[1]+" "+name[2])
            for process in group.Experiment:
                if(date0!=process):
                    date0=process
                    filtered=group[group.Experiment==process]
                    print(filtered.Depth,filtered["Time Constant"])
##                    try:
##                        params,covs = curve_fit(totalFun,filtered["Depth"],
##                                                (filtered["Time Constant"]),
##                                                maxfev=100000,p0=[1,1,1,-2,1,10],
##                                                bounds=((0,0,0,-inf,-inf,-4000),(inf,inf,inf,inf,inf,4000)),
##                                                check_finite=False)
##                        axs[0,0].plot(np.linspace(min(filtered.Depth),max(filtered.Depth)),totalFun(np.linspace(min(filtered.Depth),max(filtered.Depth)),*params))
##                        print("Al %s"%params)
##                    except:
##                        print("type error")                       
        elif("Au" in name[1]):
            axs[0,1].plot(group.Depth,group["Time Constant"],marker='o',ls=' ')
            axs[0,1].set_title(name[1]+" "+name[2])
            for process in group.Experiment:
                if(date0!=process):
                    date0=process
                    filtered=group[group.Experiment==process]
##                    try:
##                        params,covs = curve_fit(totalFun,filtered["Depth"],filtered["Time Constant"],maxfev=100000,p0=[1,1,1,-2,1,10],bounds=((0,0,0,-inf,-inf,-4000),(inf,inf,inf,inf,inf,4000)))
##                        axs[0,1].plot(np.linspace(min(filtered.Depth),max(filtered.Depth)),totalFun(np.linspace(min(filtered.Depth),max(filtered.Depth)),*params))
##                        print("Au %s"%params)
##                    except:
##                            print("type error") 
        elif("Cu" in name[1]):
            axs[0,2].plot(group.Depth,group["Time Constant"],marker='o',ls=' ')
            axs[0,2].set_title(name[1]+" "+name[2])
            for process in group.Experiment:
                if(date0!=process):
                    date0=process
                    filtered=group[group.Experiment==process]
##                    try:
##                        params,covs = curve_fit(totalFun,filtered["Depth"],filtered["Time Constant"],maxfev=100000,p0=[1,1,1,-2,1,10],bounds=((0,0,0,-inf,-inf,-4000),(inf,inf,inf,inf,inf,4000)))
##                        axs[0,2].plot(np.linspace(min(filtered.Depth),max(filtered.Depth)),totalFun(np.linspace(min(filtered.Depth),max(filtered.Depth)),*params))
##                        print("Cu %s"%params)
##                    except:
##                        print("type error") 
##    elif("Conic" in name):
##        if("Al" in name[1]):
##            axs[0,3].plot(group.Depth,group["Time Constant"],marker='o',ls=' ',label=str(name[3]))
##            axs[0,3].set_title(name[1]+" "+name[2])
##        elif("Au" in name[1]):
##            axs[0,4].plot(group.Depth,group["Time Constant"],marker='o',ls=' ',label=str(name[3]))
##            axs[0,4].set_title(name[1]+" "+name[2])
##        elif("Cu" in name[1]):
##            axs[0,5].plot(group.Depth,group["Time Constant"],marker='o',ls=' ',label=str(name[1]))
##            axs[0,5].set_title(name[1]+" "+name[2])
axs[0,0].set_xlabel("Depth (um)")
axs[0,0].set_ylabel("Time Constant (ms)")
axs[0,1].set_xlabel("Depth (um)")
axs[0,1].set_ylabel("Time Constant (ms)")
axs[0,2].set_xlabel("Depth (um)")
axs[0,2].set_ylabel("Time Constant (ms)")
##axs[0,3].set_xlabel("Depth (um)")
##axs[0,3].set_ylabel("Time Constant (ms)")
##axs[0,4].set_xlabel("Depth (um)")
##axs[0,4].set_ylabel("Time Constant (ms)")
##axs[0,5].set_xlabel("Depth (um)")
##axs[0,5].set_ylabel("Time Constant (ms)")
data.sort_values(by="Load")
gurup=data.groupby(["Tip","Sample","Sequence","Procedure","Experiment"])
for name,group in gurup:
    if(seqName in name and "Conic" in name and processName in name[3]):
        fig.suptitle(name[0]+" "+name[3])
        if("Al" in name[1]):
            axs[1,0].plot(group.Load,group["Time Constant"],marker='o',ls=' ')
            axs[1,0].set_title(name[1]+" "+name[2])
        elif("Au" in name[1]):
            axs[1,1].plot(group.Load,group["Time Constant"],marker='o',ls=' ')
            axs[1,1].set_title(name[1]+" "+name[2])
        elif("Cu" in name[1]):
            axs[1,2].plot(group.Load,group["Time Constant"],marker='o',ls=' ')
            axs[1,2].set_title(name[1]+" "+name[2])
##    elif("Conic" in name):
##        if("Al" in name[1]):
##            axs[1,3].plot(group.Load,group["Time Constant"],marker='o',ls=' ',label=str(name[3]))
##            axs[1,3].set_title(name[1]+" "+name[2])
##        elif("Au" in name[1]):
##            axs[1,4].plot(group.Load,group["Time Constant"],marker='o',ls=' ',label=str(name[3]))
##            axs[1,4].set_title(name[1]+" "+name[2])
##        elif("Cu" in name[1]):
##            axs[1,5].plot(group.Load,group["Time Constant"],marker='o',ls=' ',label=str(name[1]))
##            axs[1,5].set_title(name[1]+" "+name[2])
axs[1,0].set_xlabel("Load (mN)")
axs[1,0].set_ylabel("Time Constant (ms)")
axs[1,1].set_xlabel("Load (mN)")
axs[1,1].set_ylabel("Time Constant (ms)")
axs[1,2].set_xlabel("Load (mN)")
axs[1,2].set_ylabel("Time Constant (ms)")
##axs[1,3].set_xlabel("Load (mN)")
##axs[1,3].set_ylabel("Time Constant (ms)")
##axs[1,4].set_xlabel("Load (mN)")
##axs[1,4].set_ylabel("Time Constant (ms)")
##axs[1,5].set_xlabel("Load (mN)")
##axs[1,5].set_ylabel("Time Constant (ms)")
axs[0,0].legend()
axs[0,1].legend()
axs[0,2].legend()
##axs[0,3].legend()
##axs[0,4].legend()
##axs[0,5].legend()
axs[1,0].legend()
axs[1,1].legend()
axs[1,2].legend()
##axs[1,3].legend()
##axs[1,4].legend()
##axs[1,5].legend()
plt.show()
data.sort_values(by="Depth")
gurup=data.groupby(["Tip","Sample","Sequence","Procedure","Experiment"])
fig,axs=plt.subplots(2,3)
date0=""
for name,group in gurup:
    if(seqName in name and "Berkovich" in name and processName in name[3]):
        if("Al" in name[1]):
            axs[0,0].plot(group.Depth,group["Time Constant"],marker='o',ls=' ')
            axs[0,0].set_title(name[1]+" "+name[2])
            for process in group.Experiment:
                if(date0!=process):
                    date0=process
                    filtered=group[group.Experiment==process]
                    print(filtered.Depth,filtered["Time Constant"])
##                    try:
##                        params,covs = curve_fit(totalFun,filtered["Depth"],
##                                                (filtered["Time Constant"]),
##                                                maxfev=100000,p0=[1,1,1,-2,1,10],
##                                                bounds=((0,0,0,-inf,-inf,-4000),(inf,inf,inf,inf,inf,4000)),
##                                                check_finite=False)
##                        axs[0,0].plot(np.linspace(min(filtered.Depth),max(filtered.Depth)),totalFun(np.linspace(min(filtered.Depth),max(filtered.Depth)),*params))
##                        print("Al %s"%params)
##                    except:
##                        print("type error")                       
        elif("Au" in name[1]):
            axs[0,1].plot(group.Depth,group["Time Constant"],marker='o',ls=' ')
            axs[0,1].set_title(name[1]+" "+name[2])
            for process in group.Experiment:
                if(date0!=process):
                    date0=process
                    filtered=group[group.Experiment==process]
##                    try:
##                        params,covs = curve_fit(totalFun,filtered["Depth"],filtered["Time Constant"],maxfev=100000,p0=[1,1,1,-2,1,10],bounds=((0,0,0,-inf,-inf,-4000),(inf,inf,inf,inf,inf,4000)))
##                        axs[0,1].plot(np.linspace(min(filtered.Depth),max(filtered.Depth)),totalFun(np.linspace(min(filtered.Depth),max(filtered.Depth)),*params))
##                        print("Au %s"%params)
##                    except:
##                            print("type error") 
        elif("Cu" in name[1]):
            axs[0,2].plot(group.Depth,group["Time Constant"],marker='o',ls=' ')
            axs[0,2].set_title(name[1]+" "+name[2])
            for process in group.Experiment:
                if(date0!=process):
                    date0=process
                    filtered=group[group.Experiment==process]
##                    try:
##                        params,covs = curve_fit(totalFun,filtered["Depth"],filtered["Time Constant"],maxfev=100000,p0=[1,1,1,-2,1,10],bounds=((0,0,0,-inf,-inf,-4000),(inf,inf,inf,inf,inf,4000)))
##                        axs[0,2].plot(np.linspace(min(filtered.Depth),max(filtered.Depth)),totalFun(np.linspace(min(filtered.Depth),max(filtered.Depth)),*params))
##                        print("Cu %s"%params)
##                    except:
##                        print("type error") 
##    elif("Conic" in name):
##        if("Al" in name[1]):
##            axs[0,3].plot(group.Depth,group["Time Constant"],marker='o',ls=' ',label=str(name[3]))
##            axs[0,3].set_title(name[1]+" "+name[2])
##        elif("Au" in name[1]):
##            axs[0,4].plot(group.Depth,group["Time Constant"],marker='o',ls=' ',label=str(name[3]))
##            axs[0,4].set_title(name[1]+" "+name[2])
##        elif("Cu" in name[1]):
##            axs[0,5].plot(group.Depth,group["Time Constant"],marker='o',ls=' ',label=str(name[1]))
##            axs[0,5].set_title(name[1]+" "+name[2])
axs[0,0].set_xlabel("Depth (um)")
axs[0,0].set_ylabel("Time Constant (ms)")
axs[0,1].set_xlabel("Depth (um)")
axs[0,1].set_ylabel("Time Constant (ms)")
axs[0,2].set_xlabel("Depth (um)")
axs[0,2].set_ylabel("Time Constant (ms)")
##axs[0,3].set_xlabel("Depth (um)")
##axs[0,3].set_ylabel("Time Constant (ms)")
##axs[0,4].set_xlabel("Depth (um)")
##axs[0,4].set_ylabel("Time Constant (ms)")
##axs[0,5].set_xlabel("Depth (um)")
##axs[0,5].set_ylabel("Time Constant (ms)")
data.sort_values(by="Load")
gurup=data.groupby(["Tip","Sample","Sequence","Procedure","Experiment"])
for name,group in gurup:
    if(seqName in name and "Berkovich" in name and processName in name[3]):
        fig.suptitle(name[0]+" "+name[3])
        if("Al" in name[1]):
            axs[1,0].plot(group.Load,group["Time Constant"],marker='o',ls=' ')
            axs[1,0].set_title(name[1]+" "+name[2])
        elif("Au" in name[1]):
            axs[1,1].plot(group.Load,group["Time Constant"],marker='o',ls=' ')
            axs[1,1].set_title(name[1]+" "+name[2])
        elif("Cu" in name[1]):
            axs[1,2].plot(group.Load,group["Time Constant"],marker='o',ls=' ')
            axs[1,2].set_title(name[1]+" "+name[2])
##    elif("Conic" in name):
##        if("Al" in name[1]):
##            axs[1,3].plot(group.Load,group["Time Constant"],marker='o',ls=' ',label=str(name[3]))
##            axs[1,3].set_title(name[1]+" "+name[2])
##        elif("Au" in name[1]):
##            axs[1,4].plot(group.Load,group["Time Constant"],marker='o',ls=' ',label=str(name[3]))
##            axs[1,4].set_title(name[1]+" "+name[2])
##        elif("Cu" in name[1]):
##            axs[1,5].plot(group.Load,group["Time Constant"],marker='o',ls=' ',label=str(name[1]))
##            axs[1,5].set_title(name[1]+" "+name[2])
axs[1,0].set_xlabel("Load (mN)")
axs[1,0].set_ylabel("Time Constant (ms)")
axs[1,1].set_xlabel("Load (mN)")
axs[1,1].set_ylabel("Time Constant (ms)")
axs[1,2].set_xlabel("Load (mN)")
axs[1,2].set_ylabel("Time Constant (ms)")
##axs[1,3].set_xlabel("Load (mN)")
##axs[1,3].set_ylabel("Time Constant (ms)")
##axs[1,4].set_xlabel("Load (mN)")
##axs[1,4].set_ylabel("Time Constant (ms)")
##axs[1,5].set_xlabel("Load (mN)")
##axs[1,5].set_ylabel("Time Constant (ms)")
axs[0,0].legend()
axs[0,1].legend()
axs[0,2].legend()
##axs[0,3].legend()
##axs[0,4].legend()
##axs[0,5].legend()
axs[1,0].legend()
axs[1,1].legend()
axs[1,2].legend()
##axs[1,3].legend()
##axs[1,4].legend()
##axs[1,5].legend()
plt.show()
                    
