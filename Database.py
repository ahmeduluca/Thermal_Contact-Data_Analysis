import pandas as pd
import numpy as np
import glob
import os
from scipy.special import iv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import pearsonr, spearmanr
from io import open
import math
from numpy import inf,tanh,tan,sqrt
mainfile="/Volumes/AhmedUluca/Indenter_Data/SELECTIONS"

tipfiles=[f.path for f in os.scandir(mainfile) if f.is_dir()]

def noop(x):
    return x

def inverse(x,y=1):
    return y/x

def sqrtArea(d):
    return sqrt(TrunCont(d))

def lineer(x,m,c):
    return m*x+c

def MacroPic(x,m,n):
    return m*x**n

def cylinder_etha(h, k, d, l):
    m=2*sqrt(h/(k*d))
    return tanh((l+d/4)*m)/((l+d/4)*m)

def Spreading(contact_dia,k1,k2):
    return contact_dia*(2*k1*k2/(k1+k2))

def quad(d,c1):
#    d=d*10**(-6)
    return c1/((d)**2)

def oneOver(d,c1,c2=0):
#    d=d*10**(-6)
    return c1*(d+c2)**-1

def cylinderEff(d,c1,c2):
    return (c1/tanh(c2*(d+2.5e-6)))*(1+2.5e-6/d)      

def coneEff(d,c1,c2):
    return c1/(d**0.5*iv(2,c2*sqrt(d))/iv(1,c2*sqrt(d)))       

def comply(d0,l0,T,Cf=1.864,Tcal=23,cT=0.63): ##for calibration at 23 Celcius and experiment temperature T, compliance:
    if(T<23):
        T=23
    if(l0<0):
        l0=0
    return d0-l0*(Cf+(T-Tcal)*cT)*10**-9

def TrunProject(d,R_0=5.0e-6,alpha=60.):
    alpha=math.radians(alpha)
    return np.pi*(R_0+tan(alpha)*d)**2

def TrunCont(d,d_0=2.9e-6,R_0=5.0e-6,alpha=60):
    alpha=math.radians(alpha)
    return np.pi*((((d+d_0)**2*tan(alpha)-d_0*R_0)/np.cos(alpha))+R_0**2)

def Pressure(f,d):
    return f/TrunProject(d)

def Conductance(t,d,C=1.84E-5):
    return C/(t*TrunCont(d))



def choose(x):
    if(x=="lineer"):
        name="AreaFit1-berko"
        labelx="Contact Area ($m^2$)"
        labely="Inverse Time Constant (1/s)"
        parLin=[1,1]
        boundLin=((0,-inf),(inf,inf))
        return parLin,boundLin,lineer,name,labelx,labely
    elif(x=="inverse"):
        name="AreaFits-berko"
        labelx="Contact Area ($m^2$)"
        labely="Inverse Time Constant (1/s)"
        parLin=[1,1]
        boundLin=((0,-inf),(inf,inf))
        return parLin,boundLin,lineer,name,labelx,labely
    elif(x=="quad"):
#        name="QuadraticFits"
        name="QuadraticFits-area-berko"
#        labelx="Depth (m)"
        labelx="Square Root Contact Area (m)"
        labely="Time Constant (s)"
        parQuad=[1]
        boundQuad=((0),(inf))
        return parQuad,boundQuad,quad,name,labelx,labely
    elif(x=="ConeEff"):
 #       name="Bessel-Fits-berko"
        name="Bessel-Fits-area"
#        labelx="Depth (m)"
        labelx="Square Root Contact Area (m)"
        labely="Time Constant (s)"
        parEff=[1,1]
        boundEff=((0,0),(inf,inf))
        return parEff,boundEff,coneEff,name,labelx,labely
    elif(x=="CylinderEff"):
        name="Tanh-Fits-berko"
##        name="Tanh-Fits-area-berko"
        labelx="Depth (m)"
##        labelx="Square Root Contact Area (m)"
        labely="Time Constant (s)"
        parEff=[1,1]
        boundEff=((0,0),(inf,inf))
        return parEff,boundEff,cylinderEff,name,labelx,labely
    elif(x=="macro"):
        name="Pressure-Fits-berko"
        labelx="Pressure (Pa)"
        labely="Conductance (W/$m^2$K)"
        parMacro=[1,0.85]
        boundMacro=((0,0.0),(inf,3))
        return parMacro,boundMacro,MacroPic,name,labelx,labely
    elif(x=="1over"):
        name="ConstrictFits-berko"
##        name="ConstrictFits-area-berko"
        labelx="Depth (m)"
##        labelx="Square Root Contact Area (m)"
        labely="Time Constant (s)"
        par1over=[1]
        bound1over=((0),(100000))
        return par1over,bound1over,oneOver,name,labelx,labely

pars,bound,fonk,fitname,labelx,labely=choose("1over")
sigmaAbs=True
varian=np.inf
givenx=noop
giveny=noop

tipName="Berkovich-like"
seqName="Contact"
processName="Sweep"
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
"Experiment":os.path.basename(experiment).replace("_"," "),"Sequence":"Contact","Depth":float(line.split()[3])*10**-6,
"Load":float(line.split()[4])*10**-3,"Time Constant":float(line.split()[0])*0.001,"Sigma":float(line.split()[1])*0.001,"Temperature":float(line.split()[2])})
                    elif("decayTC_disp_load" in text):
                        with open(text,"r",encoding="utf-8") as f:
                            lines = f.readlines()[1:]
                            for line in lines:
                                if(float(line.split()[0])*0.001>0):
                                    expresults.append({"Tip":os.path.basename(tip),
    "Sample":os.path.basename(sample),"Procedure":os.path.basename(procedure),
    "Experiment":os.path.basename(experiment).replace("_"," "),"Sequence":"Separation","Depth":float(line.split()[3])*10**-6,
    "Load":float(line.split()[4])*10**-3,"Time Constant":float(line.split()[0])*0.001,"Sigma":float(line.split()[1])*0.001,"Temperature":float(line.split()[2])})
    data=pd.DataFrame(expresults)
for i in range (len(data['Depth'])):
    data['Depth'][i]=comply(data['Depth'][i],data['Load'][i],23)    
##data.to_excel(os.path.join(mainfile,"AllResults-3.xlsx"),engine='xlsxwriter')
#data.groupby("Sequence").plot.hexbin(x='Load',y='Time Constant',C='Depth',subplots='True')
#data.groupby("Sequence").plot.hexbin(x='Depth',y='Time Constant',C='Load',subplots='True')
#data.plot(x='Load',y='Time Constant')
#plt.show()
#ax.clabel("Depth")
#data.groupby("Sample").plot.hexbin(x='Depth',y='Time Constant',C='Load',subplots='True')
tidy=data.groupby(["Experiment"])
##for name,group in tidy:
##    minval=min(group.Depth)
##    for i,row in data[data.Experiment==name].iterrows():
##        data["Depth"][i]-=minval
#data.sort_values(by="Depth")
date0=""

                            ###########  CORRELEATIONS  ###############

##depCon=[]
##loadCon=[]
##tempCon=[]
##depSep=[]
##loadSep=[]
##tempSep=[]
##sigmaToCorr=[]
##
##gurup=data.groupby(["Tip","Sample","Experiment"])
##for name, group in gurup:
##    for exp in group.Experiment:
##        if(date0!=exp):
##            date0=exp
##            filtered=group[group.Experiment==exp]
##            if(group[group['Sequence']=='Contact']['Time Constant'].count()==group[group['Sequence']=='Separation']['Time Constant'].count()):               
##                for j in range (len(filtered.Sigma)):
##                    if(filtered.Sequence.values[j]=="Contact"):
##                        depCon.append(filtered.Depth.values[j])
##                        loadCon.append(filtered.Load.values[j])
##                        sigmaToCorr.append(100*filtered.Sigma.values[j]/filtered["Time Constant"].values[j])
##                        tempCon.append(filtered["Temperature"].values[j])
##                    else:
##                        depSep.append(filtered.Depth.values[j])
##                        loadSep.append(filtered.Load.values[j])
##                        tempSep.append(filtered["Temperature"].values[j])

##try:
##    print(pearsonr(np.array(depCon)-np.array(depSep),sigmaToCorr))
##except:
##    print("problem on depth correlation")
##try:
##    correl,p_val=spearmanr(np.array(depCon)-np.array(depSep),sigmaToCorr)
##    plt.plot(np.array(depCon)-np.array(depSep),sigmaToCorr,'.',label="Correlation = %.3f, p_0 = %.3e"%(correl,p_val))
##    plt.legend()
##    plt.xlabel(r'Amplitude ($\mu$m)')
##    plt.ylabel('% Uncertainty')
##    plt.show()
##except:
##    print("problem on depth correlation-spear")
##try:
##    print(pearsonr(np.array(loadCon)-np.array(loadSep),sigmaToCorr))
##except:
##    print("problem on load correlation")
##try:
##    correl,p_val=spearmanr(np.array(loadCon)-np.array(loadSep),sigmaToCorr)
##    plt.plot(np.array(loadCon)-np.array(loadSep),sigmaToCorr,'.',label="Correlation = %.3f, p_0 = %.3e"%(correl,p_val))
##    plt.legend()
##    plt.xlabel(r'$\Delta$Load (mN)')
##    plt.ylabel('% Uncertainty')
##    plt.show()
##except:
##    print("problem on load correlation-spear")
##try:
##    print(pearsonr(np.array(tempCon)-np.array(tempSep),sigmaToCorr))
##except:
##    print("problem on temp correlation")
##try:
##    correl,p_val=spearmanr(np.array(tempCon)-np.array(tempSep),sigmaToCorr)
##    plt.plot(np.array(tempCon)-np.array(tempSep),sigmaToCorr,'.',label="Correlation = %.3f, p_0 = %.3e"%(correl,p_val))
##    plt.legend()
##    plt.xlabel(r'$\Delta$Temperature ($^\circ$C)')
##    plt.ylabel('% Uncertainty')
##    plt.show()
##except:
##    print("problem on temp correlation-spear")

                        ####################################################
                        
gurup=data.groupby(["Tip","Sample","Sequence","Procedure","Experiment"])
ax=plt.gca()
fig,axs=plt.subplots(2,3)
fitresults=[]
AlSum=0
AuSum=0
CuSum=0
nAl=0
nAu=0
nCu=0
maxCu=0
maxAu=0
maxAl=0
minCu=np.inf
minAu=np.inf
minAl=np.inf
colors=[]
for name,group in gurup:
    if(seqName in name and tipName in name and processName in name[3]):
        if("Al" in name[1]):              
            for process in group.Experiment:
                if(date0!=process):
                    date0=process
                    filtered=group[group.Experiment==process]
                    try:
                        params,covs = curve_fit(fonk,givenx(filtered["Depth"]),
                                                giveny(filtered["Time Constant"]),sigma=filtered["Sigma"],absolute_sigma=sigmaAbs,maxfev=100000,p0=pars,bounds=bound)
                        if(covs[0,0]**0.5/params[0]<varian):
                            color=next(ax._get_lines.prop_cycler)['color']
                            colors.append(color)
                            axs[0,0].plot(np.linspace(min(givenx(filtered.Depth)),max(givenx(filtered.Depth))),fonk(np.linspace(min(givenx(filtered.Depth)),max(givenx(filtered.Depth))),*params),label=str(params),color=color)
                            print("Al %s"%params)
                            fitresults.append({"Experiment":filtered["Experiment"].values[0],"Tip":filtered.Tip.values[0],"Sample":filtered.Sample.values[0],
                                           "Procedure":filtered.Procedure.values[0],"Sequence":filtered.Sequence.values[0],"c1":params[0],"sigma1":covs[0,0]**0.5})#,"c2":params[1],"sigma2":covs[1,1]**0.5})#,"c3":params[2]})
    #                            axs[0,0].errorbar(givenx(group.Depth),giveny(group["Time Constant"]),yerr=group.Sigma,xerr=givenx(0.2*np.ones(len(group.Depth.values))),marker='o',ls=' ',color=color)
                            axs[0,0].plot(givenx(group.Depth),giveny(group["Time Constant"]),marker='o',ls=' ',color=color)
                            axs[0,0].set_title(name[1]+" "+name[2])
                            AlSum+=sum(group["Time Constant"])
                            nAl+=len(group.Depth.values)
                            if(min(group["Time Constant"])<minAl):
                                minAl=min(group["Time Constant"])
                            if(max(group["Time Constant"])>maxAl):
                                maxAl=max(group["Time Constant"])
                            axs[0,0].set_ylim([minAl-0.1,maxAl+0.1])
                        else:
                            print(covs[0,0]**0.5/params[0])
                    except:
                        print("type error")                       
        elif("Au" in name[1]):
            for process in group.Experiment:
                if(date0!=process):
                    date0=process
                    filtered=group[group.Experiment==process]
                    try:
                        params,covs = curve_fit(fonk,givenx(filtered["Depth"]),
                                                giveny(filtered["Time Constant"]),sigma=filtered["Sigma"],absolute_sigma=sigmaAbs,maxfev=100000,p0=pars,bounds=bound)
                        if(covs[0,0]**0.5/params[0]<varian):
                            color=next(ax._get_lines.prop_cycler)['color']
                            colors.append(color)
                            axs[0,1].plot(np.linspace(min(givenx(filtered.Depth)),max(givenx(filtered.Depth))),fonk(np.linspace(min(givenx(filtered.Depth)),max(givenx(filtered.Depth))),*params),label=str(params),color=color)
                            print("Au %s"%params)
                            fitresults.append({"Experiment":filtered["Experiment"].values[0],"Tip":filtered.Tip.values[0],"Sample":filtered.Sample.values[0],
                                           "Procedure":filtered.Procedure.values[0],"Sequence":filtered.Sequence.values[0],"c1":params[0],"sigma1":covs[0,0]**0.5})#,"c2":params[1],"sigma2":covs[1,1]**0.5})#,"c3":params[2]})
#                            axs[0,1].errorbar(givenx(group.Depth),giveny(group["Time Constant"]),yerr=group.Sigma,xerr=givenx(0.2*np.ones(len(group.Depth.values))),marker='o',ls=' ',color=color)
                            axs[0,1].plot(givenx(group.Depth),giveny(group["Time Constant"]),marker='o',ls=' ',color=color)
                            axs[0,1].set_title(name[1]+" "+name[2])
                            AuSum+=sum(group["Time Constant"])
                            nAu+=len(group.Depth.values)
                            if(min(group["Time Constant"])<minAu):
                                minAu=min(group["Time Constant"])
                            if(max(group["Time Constant"])>maxAu):
                                maxAu=max(group["Time Constant"])
                            axs[0,1].set_ylim([minAu-0.1,maxAu+0.1])
                        else:
                            print(covs[0,0]**0.5/params[0])
                    except:
                            print("type error") 
        elif("Cu" in name[1]):
            for process in group.Experiment:
                if(date0!=process):
                    date0=process
                    filtered=group[group.Experiment==process]
                    try:
                        params,covs = curve_fit(fonk,givenx(filtered["Depth"]),
                                                giveny(filtered["Time Constant"]),sigma=filtered["Sigma"],absolute_sigma=sigmaAbs,maxfev=100000,p0=pars,bounds=bound)
                        if(covs[0,0]**0.5/params[0]<varian):
                            color=next(ax._get_lines.prop_cycler)['color']
                            colors.append(color)
                            axs[0,2].plot(np.linspace(min(givenx(filtered.Depth)),max(givenx(filtered.Depth))),fonk(np.linspace(min(givenx(filtered.Depth)),max(givenx(filtered.Depth))),*params),label=str(params),color=color)
                            print("Cu %s"%params)
                            fitresults.append({"Experiment":filtered["Experiment"].values[0],"Tip":filtered.Tip.values[0],"Sample":filtered.Sample.values[0],
                                           "Procedure":filtered.Procedure.values[0],"Sequence":filtered.Sequence.values[0],"c1":params[0],"sigma1":covs[0,0]**0.5})#,"c2":params[1],"sigma2":covs[1,1]**0.5})#,"c3":params[2]})
#                            axs[0,2].errorbar(givenx(group.Depth),giveny(group["Time Constant"]),yerr=group.Sigma,xerr=givenx(0.2*np.ones(len(group.Depth.values))),marker='o',ls=' ',color=color)
                            axs[0,2].plot(givenx(group.Depth),giveny(group["Time Constant"]),marker='o',ls=' ',color=color)
#                            axs[0,2].plot(givenx(group.Depth),giveny(group["Time Constant"]),marker='o',ls=' ',color=color)
                            axs[0,2].set_title(name[1]+" "+name[2])
                            CuSum+=sum(group["Time Constant"])
                            nCu+=len(group.Depth.values)
                            if(min(group["Time Constant"])<minCu):
                                minCu=min(group["Time Constant"])
                            if(max(group["Time Constant"])>maxCu):
                                maxCu=max(group["Time Constant"])
                            axs[0,2].set_ylim([minCu-0.1,maxCu+0.1])
                        else:
                            print(covs[0,0]**0.5/params[0])
                    except:
                            print("type error") 
print("For %s Tip:"%(tipName))
if(nAl>0):
    print("Mean Tc for Al: %.3e"%(AlSum/nAl))
    print(nAl)
if(nAu>0):
    print("Mean Tc for Au: %.3e"%(AuSum/nAu))
    print(nAu)
if(nCu>0):    
    print("Mean Tc for Cu: %.3e"%(CuSum/nCu))
    print(nCu)
print("Max Tc for Al: %.3e"%(maxAl))
print("Max Tc for Au: %.3e"%(maxAu))
print("Max Tc for Cu: %.3e"%(maxCu))
print("Min Tc for Al: %.3e"%(minAl))
print("Min Tc for Au: %.3e"%(minAu))
print("Min Tc for Cu: %.3e\n"%(minCu))

axs[0,0].set_xlabel(labelx)
axs[0,0].set_ylabel(labely)
axs[0,1].set_xlabel(labelx)
axs[0,1].set_ylabel(labely)
axs[0,2].set_xlabel(labelx)
axs[0,2].set_ylabel(labely)
fits=pd.DataFrame(fitresults)

fits.to_excel(os.path.join(mainfile,"%s-v2.xlsx"%fitname),engine='xlsxwriter')

                                ##############################################
fitresults=[]
gurup=data.groupby(["Tip","Sample","Sequence","Procedure","Experiment"])
j=0
pars,bound,fonk,fitname,labelx,labely=choose("macro")
sigmaAbs=False
varian=np.inf
givenx=Pressure
giveny=Conductance
date0=""
palx=[]
paly=[]
paux=[]
pauy=[]
pcux=[]
pcuy=[]
for name,group in gurup:
    if(seqName in name and tipName in name and processName in name[3]):
        fig.suptitle(name[0]+" "+name[3])
        if("Al" in name[1]):
            for exp in group.Experiment:
                if(date0!=exp):
                    date0=exp
                    filtered=group[group.Experiment==exp]
                    try:
                        params,covs = curve_fit(fonk,givenx(group.Load,group.Depth), giveny(group["Time Constant"],group.Depth),sigma=filtered["Sigma"],absolute_sigma=sigmaAbs,maxfev=100000,p0=pars,bounds=bound)
                        if(covs[0,0]**0.5/params[0]<varian):
                            axs[1,0].plot(np.linspace(min(givenx(group.Load,group.Depth)),max(givenx(group.Load,group.Depth))),fonk(np.linspace(min(givenx(group.Load,group.Depth)),max(givenx(group.Load,group.Depth))),*params),label=str(params),color=colors[j])
        ##                        axs[1,0].errorbar(group.Load,group["Time Constant"],yerr=group.Sigma,xerr=2.5*np.ones(len(group.Load.values)),marker='o',ls=' ')
                            axs[1,0].plot(givenx(group.Load,group.Depth),giveny(group["Time Constant"],group.Depth),marker='o',ls=' ',color=colors[j])
        #                        axs[1,0].set_title(name[1]+" "+name[2])
                            fitresults.append({"Experiment":filtered["Experiment"].values[0],"Tip":filtered.Tip.values[0],"Sample":filtered.Sample.values[0],
                                           "Procedure":filtered.Procedure.values[0],"Sequence":filtered.Sequence.values[0],"c1":params[0],"sigma1":covs[0,0]**0.5,"c2":params[1],"sigma2":covs[1,1]**0.5})#,"c3":params[2]})
#
                    except:
                        print("Load Fit error")
        elif("Au" in name[1]):
            for exp in group.Experiment:
                if(date0!=exp):
                    date0=exp
                    filtered=group[group.Experiment==exp]
                    try:
                        params,covs = curve_fit(fonk,givenx(group.Load,group.Depth), giveny(group["Time Constant"],group.Depth),sigma=filtered["Sigma"],absolute_sigma=sigmaAbs,maxfev=100000,p0=pars,bounds=bound)
                        if(covs[0,0]**0.5/params[0]<varian):
                            axs[1,1].plot(np.linspace(min(givenx(group.Load,group.Depth)),max(givenx(group.Load,group.Depth))),fonk(np.linspace(min(givenx(group.Load,group.Depth)),max(givenx(group.Load,group.Depth))),*params),label=str(params),color=colors[j])
##                        axs[1,0].errorbar(group.Load,group["Time Constant"],yerr=group.Sigma,xerr=2.5*np.ones(len(group.Load.values)),marker='o',ls=' ')
                            axs[1,1].plot(givenx(group.Load,group.Depth),giveny(group["Time Constant"],group.Depth),marker='o',ls=' ',color=colors[j])
#                        axs[1,0].set_title(name[1]+" "+name[2])
                            fitresults.append({"Experiment":filtered["Experiment"].values[0],"Tip":filtered.Tip.values[0],"Sample":filtered.Sample.values[0],
                                           "Procedure":filtered.Procedure.values[0],"Sequence":filtered.Sequence.values[0],"c1":params[0],"sigma1":covs[0,0]**0.5,"c2":params[1],"sigma2":covs[1,1]**0.5})#,"c3":params[2]})
# 
                    except:
                        print("Load Fit error")
        elif("Cu" in name[1]):
            for exp in group.Experiment:
                if(date0!=exp):
                    date0=exp
                    filtered=group[group.Experiment==exp]
                    try:
                        params,covs = curve_fit(fonk,givenx(group.Load,group.Depth), giveny(group["Time Constant"],group.Depth),sigma=filtered["Sigma"],absolute_sigma=sigmaAbs,maxfev=100000,p0=pars,bounds=bound)
                        if(covs[0,0]**0.5/params[0]<varian):
                            axs[1,2].plot(np.linspace(min(givenx(group.Load,group.Depth)),max(givenx(group.Load,group.Depth))),fonk(np.linspace(min(givenx(group.Load,group.Depth)),max(givenx(group.Load,group.Depth))),*params),label=str(params),color=colors[j])
        ##                        axs[1,0].errorbar(group.Load,group["Time Constant"],yerr=group.Sigma,xerr=2.5*np.ones(len(group.Load.values)),marker='o',ls=' ')
                            axs[1,2].plot(givenx(group.Load,group.Depth),giveny(group["Time Constant"],group.Depth),marker='o',ls=' ',color=colors[j])
        #                        axs[1,0].set_title(name[1]+" "+name[2])
                            fitresults.append({"Experiment":filtered["Experiment"].values[0],"Tip":filtered.Tip.values[0],"Sample":filtered.Sample.values[0],
                                           "Procedure":filtered.Procedure.values[0],"Sequence":filtered.Sequence.values[0],"c1":params[0],"sigma1":covs[0,0]**0.5,"c2":params[1],"sigma2":covs[1,1]**0.5})#,"c3":params[2]})
# 
                    except:
                        print("Load Fit error")
        j+=1
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
axs[1,0].set_xlabel(labelx)
axs[1,0].set_ylabel(labely)
axs[1,1].set_xlabel(labelx)
axs[1,1].set_ylabel(labely)
axs[1,2].set_xlabel(labelx)
axs[1,2].set_ylabel(labely)

palx=givenx(data[(data.Tip==tipName) & (data.Sample=="Al")].Load.values,
            data[(data.Tip==tipName) & (data.Sample=="Al")].Depth.values)
paly=giveny(data[(data.Tip==tipName) & (data.Sample=="Al")]["Time Constant"].values,
            data[(data.Tip==tipName) & (data.Sample=="Al")].Depth.values)
paux=givenx(data[(data.Tip==tipName) & (data.Sample=="Au")].Load.values,
            data[(data.Tip==tipName) & (data.Sample=="Au")].Depth.values)
pauy=giveny(data[(data.Tip==tipName) & (data.Sample=="Au")]["Time Constant"].values,
            data[(data.Tip==tipName) & (data.Sample=="Au")].Depth.values)
pcux=givenx(data[(data.Tip==tipName) & (data.Sample=="Cu")].Load.values,
            data[(data.Tip== tipName) & (data.Sample=="Cu")].Depth.values)
pcuy=giveny(data[(data.Tip==tipName) & (data.Sample=="Cu")]["Time Constant"].values,
            data[(data.Tip==tipName) & (data.Sample=="Cu")].Depth.values)
try:
    ptot,ctot=curve_fit(fonk,palx,paly,absolute_sigma=sigmaAbs,maxfev=100000,p0=pars,bounds=bound)
    print(ptot)
    axs[1,0].plot(np.linspace(min(palx),max(palx)),fonk(np.linspace(min(palx),max(palx)),*params),label=str(params))
except:
    print("error on Al pressure fit")
try:
    ptot,ctot=curve_fit(fonk,palx,paly,absolute_sigma=sigmaAbs,maxfev=100000,p0=pars,bounds=bound)
    print(ptot)
    axs[1,1].plot(np.linspace(min(paux),max(paux)),fonk(np.linspace(min(paux),max(paux)),*params),label=str(params))
except:
    print("error on Au pressure fit")
try:
    ptot,ctot=curve_fit(fonk,paux,pauy,absolute_sigma=sigmaAbs,maxfev=100000,p0=pars,bounds=bound)
    print(ptot)
    axs[1,2].plot(np.linspace(min(pcux),max(pcux)),fonk(np.linspace(min(pcux),max(pcux)),*params),label=str(params))
except:
    print("error on Cu pressure fit")
plt.show()
fits1=pd.DataFrame(fitresults)
fits1.to_excel(os.path.join(mainfile,"%s-v2.xlsx"%fitname),engine='xlsxwriter')


                    ###########################################################
##
##data.sort_values(by="Depth")
##gurup=data.groupby(["Tip","Sample","Sequence","Procedure","Experiment"])
##
##maxAlfit=max(fits[fits.Sample=="Al"].c2/fits[fits.Sample=="Al"].c1)
##minAlfit=min(fits[fits.Sample=="Al"].c2/fits[fits.Sample=="Al"].c1)
##meanAlfit=np.mean(fits[fits.Sample=="Al"].c2/fits[fits.Sample=="Al"].c1)
##plt.errorbar(np.arange(0,len(fits[fits.Sample=="Al"].c1)),fits[fits.Sample=="Al"].c2/fits[fits.Sample=="Al"].c1,marker=".",yerr=np.array((fits[fits.Sample=="Al"].sigma2**2+(fits[fits.Sample=="Al"].sigma1/fits[fits.Sample=="Al"].c1**2)**2)**0.5),label="Al mean=%.5e"%meanAlfit,ls=" ")
##maxAlfit=max(fits[fits.Sample=="Au"].c2/fits[fits.Sample=="Au"].c1)
##minAlfit=min(fits[fits.Sample=="Au"].c2/fits[fits.Sample=="Au"].c1)
##meanAlfit=np.mean(fits[fits.Sample=="Au"].c2/fits[fits.Sample=="Au"].c1)
##plt.errorbar(np.arange(0,len(fits[fits.Sample=="Au"].c1)),fits[fits.Sample=="Au"].c2/fits[fits.Sample=="Au"].c1,marker=".",yerr=np.array((fits[fits.Sample=="Au"].sigma2**2+(fits[fits.Sample=="Au"].sigma1/fits[fits.Sample=="Au"].c1**2)**2)**0.5),label="Au mean=%.5e"%meanAlfit,ls=" ")
##maxAlfit=max(fits[fits.Sample=="Cu"].c2/fits[fits.Sample=="Cu"].c1)
##minAlfit=min(fits[fits.Sample=="Cu"].c2/fits[fits.Sample=="Cu"].c1)
##meanAlfit=np.mean(fits[fits.Sample=="Cu"].c2/fits[fits.Sample=="Cu"].c1)
##plt.errorbar(np.arange(0,len(fits[fits.Sample=="Cu"].c1)),fits[fits.Sample=="Cu"].c2/fits[fits.Sample=="Cu"].c1,marker=".",yerr=np.array((fits[fits.Sample=="Cu"].sigma2**2+(fits[fits.Sample=="Cu"].sigma1/fits[fits.Sample=="Cu"].c1**2)**2)**0.5),label="Cu mean=%.5e"%meanAlfit,ls=" ")

            ###########################################################
maxAlfit=max(1/fits[fits.Sample=="Al"].c1)
minAlfit=min(1/fits[fits.Sample=="Al"].c1)
meanAlfit=np.mean(1/fits[fits.Sample=="Al"].c1)
plt.errorbar(np.arange(0,len(fits[fits.Sample=="Al"].c1)),1/fits[fits.Sample=="Al"].c1,marker=".",yerr=np.array(fits[fits.Sample=="Al"].sigma1/fits[fits.Sample=="Al"].c1**2),label="Al mean=%.5e"%meanAlfit,ls=" ")
maxAufit=max(1/fits[fits.Sample=="Au"].c1)
minAufit=min(1/fits[fits.Sample=="Au"].c1)
meanAufit=np.mean(1/fits[fits.Sample=="Au"].c1)
plt.errorbar(np.arange(0,len(fits[fits.Sample=="Au"].c1)),1/fits[fits.Sample=="Au"].c1,marker=".",yerr=np.array(fits[fits.Sample=="Au"].sigma1/fits[fits.Sample=="Au"].c1**2),label="Au mean=%.5e"%meanAufit,ls=" ")
maxCufit=max(1/fits[fits.Sample=="Cu"].c1)
minCufit=min(1/fits[fits.Sample=="Cu"].c1)
meanCufit=np.mean(1/fits[fits.Sample=="Cu"].c1)
plt.errorbar(np.arange(0,len(fits[fits.Sample=="Cu"].c1)),1/fits[fits.Sample=="Cu"].c1,marker=".",yerr=np.array(fits[fits.Sample=="Cu"].sigma1/fits[fits.Sample=="Cu"].c1**2),label="Cu mean=%.5e"%meanCufit,ls=" ")
###########################################################

plt.legend()
#plt.ylabel("Fit Parameter (h/C)")
plt.ylabel(r"Thermal Diffusivity per Volume ($m^{-1}s^{-1}$)")
plt.show()
                            ###############################################

maxAlfit=max(fits1[fits1.Sample=="Al"].c2)
minAlfit=min(fits1[fits1.Sample=="Al"].c2)
meanAlfit=np.mean(fits1[fits1.Sample=="Al"].c2)
plt.errorbar(np.arange(0,len(fits1[fits1.Sample=="Al"].c2)),fits1[fits1.Sample=="Al"].c2,marker=".",yerr=np.array(fits1[fits1.Sample=="Al"].sigma1),label="Al mean=%.5e"%meanAlfit,ls=" ")
maxAufit=max(fits1[fits1.Sample=="Au"].c2)
minAufit=min(fits1[fits1.Sample=="Au"].c2)
meanAufit=np.mean(fits1[fits1.Sample=="Au"].c2)
plt.errorbar(np.arange(0,len(fits1[fits1.Sample=="Au"].c2)),fits1[fits1.Sample=="Au"].c2,marker=".",yerr=np.array(fits1[fits1.Sample=="Au"].sigma1),label="Au mean=%.5e"%meanAufit,ls=" ")
maxCufit=max(fits1[fits1.Sample=="Cu"].c2)
minCufit=min(fits1[fits1.Sample=="Cu"].c2)
meanCufit=np.mean(fits1[fits1.Sample=="Cu"].c2)
plt.errorbar(np.arange(0,len(fits1[fits1.Sample=="Cu"].c2)),fits1[fits1.Sample=="Cu"].c2,marker=".",yerr=np.array(fits1[fits1.Sample=="Cu"].sigma1),label="Cu mean=%.5e"%meanCufit,ls=" ")
plt.legend()
plt.ylabel("Fit Parameter c2")
plt.show()
