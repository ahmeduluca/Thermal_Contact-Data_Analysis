import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import glob

load_files=[]
position_files=[]
lvsdis_file=[[]]
load=[[]]
position=[[]]
date=[[]]
n=0
bit_to_load=4.75*9.81
loadcell_stiffness=536.5

file_list=glob.glob("/Users/ahmeduluca/Desktop/FusedSilicaIndents/Fused Silica Indentation-29-04/**/*.txt")
    
for i in file_list:
    if("Load" in i):
        date.append([])
        load_files.append(i)
        if(len(i)>36):
            date[n].append(i[-36:-9])
        else:
            date[n].append(i[-29:-9])            
        n+=1
    elif("Position" in i):
        position_files.append(i)
n=0        
for k in range(len(load_files)):
    load.append([])
    position.append([])
    f=open(load_files[k],"r")
    for x in f:
        load[k].append(float(x)*bit_to_load/1000)
  ##      plt.plot(load[k])
    for j in position_files:
        if(len(date[k])!=0):
            if(date[k][0] in j):
                g=open(j,"r")
                n=0
                try:
                    for x in g:
                        position[k].append((float(x)-load[k][n]/loadcell_stiffness)*1000)
##                        position[k].append((float(x))*1000)
                        n+=1
                    plt.plot(position[k],load[k], label=date[k])
                    plt.legend()
                    plt.show()
                    continue
                except:
                    print("error in"+load_files[k])
                    continue
                
plt.xlabel('Displacement (nm)')
plt.ylabel('Force (mN)')
plt.legend()
plt.show()
