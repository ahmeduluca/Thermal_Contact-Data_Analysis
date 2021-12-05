import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from io import open

mainfile="/Volumes/AhmedUluca/Indenter_Data/main"

tipfiles=[f.path for f in os.scandir(mainfile) if f.is_dir()]

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
data.to_excel(os.path.join(mainfile,"AllResults.xlsx"),engine='xlsxwriter')
data.groupby("Sequence").plot.hexbin(x='Load',y='Time Constant',C='Depth',subplots='True')
data.groupby("Sequence").plot.hexbin(x='Depth',y='Time Constant',C='Load',subplots='True')
#data.plot(x='Load',y='Time Constant')
#ax.clabel("Depth")
plt.show()

                    
