import matplotlib.animation as ani
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import glob
main="/Users/ahmeduluca/Desktop/download/process1/22-10-2021_14-54-01"
txtMains=glob.glob(os.path.join(main+'/*.txt'))
file_list=[f.path for f in os.scandir(main) if f.is_dir()]
file_list.sort()
stepnames=[f.name for f in os.scandir(main) if f.is_dir()]
stepnames.sort()
insteps=[]
time=[]
load=[]
temp=[]
depth=[]
inco=0
arr=[0]
osco=0
cc=0
time0=0
tt=0
osc_dloads=[]
osc_uloads=[]
count=0
count1=0
for i in txtMains:
    if("decayTC" in i):
        f=open(i,"r",encoding='utf-8')
        rows=f.readlines()[1:]
        for x in rows:
           osc_dloads.append(float(x.split()[2]))
    elif("riseTC" in i):
        f=open(i,"r",encoding='utf-8')
        rows=f.readlines()[1:]
        for x in rows:
           osc_uloads.append(float(x.split()[2]))           
for i in file_list:
    txt_files=glob.glob(os.path.join(i+'/*.txt'))
    cc+=1
    for j in txt_files:            
        if("ForceDisplacement" in j):
            f=open(j,"r",encoding='utf-8')
            rows=f.readlines()[2:]
            count=0
            for x in rows:
                if(count1 %50==0):                    
                    time.append(float(x.split()[0])+time0)
                    depth.append(float(x.split()[1]))
                    load.append(float(x.split()[2]))
                    temp.append(float(x.split()[3]))
                    inco+=1
                count1+=1
            time0=time[-1]
            arr.append(inco)
            insteps.append(stepnames[cc-1])
        if("AverageResults" in j):
            f=open(j,"r",encoding='utf-8')
            rows=f.readlines()[2:]
            tt=0
            osco=0
            for x in rows:
                tt+=1
                time.append(0.05*tt+time0)
                depth.append(float(x.split()[2]))
                temp.append(float(x.split()[1]))
                if(osco < int(len(rows)/2)):
                    load.append(osc_dloads[int(count/2)])
                else:
                    load.append(osc_uloads[int(count/2)])
                osco+=1
            count+=1
            time0=time[-1]
            arr.append(osco)
            insteps.append(stepnames[cc-1])
            
print(time0)                
colors = ['red', 'green', 'blue', 'orange',"cyan","purple","black"]
fig = plt.figure()
#plt.xticks(rotation=45, ha="right", rotation_mode="anchor") #rotate the x-axis values
#plt.subplots_adjust(bottom = 0.2, top = 0.9) #ensuring the dates (on the x-axis) fit in the screen
plt.ylabel('Temperature')
plt.xlabel('Not Scale Time')
m=0
ax=fig.add_subplot()
ax.tick_params(bottom=False)
load_line,=ax.plot([],[],'o-',lw=2)
x=[]
y=[]
def buildmebarchart(i):
    if(i==0):
        x.clear()
        y.clear()
    x.append(time[i])
    y.append(temp[i])
    load_line.set_data(x, y)
    return load_line,
print(np.arange(len(time))[-1])
animator = ani.FuncAnimation(fig, buildmebarchart, len(time), interval =5, save_count=100,frames=1000, blit=True)
plt.xticks(labels=None)
##plt.plot(time,load)
##plt.show()
animator.save(main+"/animation.gif")


