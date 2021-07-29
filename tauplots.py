# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 17:50:24 2019

@author: cemsa
"""

import matplotlib.pyplot as plt
import numpy as np

G_r=[]#Graphene Rise
G_d=[]#Graphene_Decay
Si_r=[]
Si_d=[]

with open('time_constants_Si.txt','r') as f:
        rows1=f.readlines() [1:]

        for x in rows1:
            Si_r.append(x.split()[1])
            Si_d.append(x.split()[2])
            
with open('time_constants_G.txt','r') as p:
        rows2=p.readlines() [1:]

        for x in rows2:
            G_r.append(x.split()[1])
            G_d.append(x.split()[2])



x=np.linspace(1,len(G_r),8)

y=[]
def pm(x,y,ort):#plotmean(pm) gives a array with the dimensions of
    #the plotted x data which has the mean value for every row
    y=np.arange(len(x))
    for i in y:
        y[i]=np.mean(ort)
    return y


#####SUBPLOTS########   

#plt.figure(figsize=(14,14))
#
#ax1=plt.subplot(2,2,1)
#plt.scatter(x, Si_r,c='r')
#plt.hold('on')
#plt.plot(x,pm(x,y,Si_r),c='black')
#ax1.set_title('Silicon_rise')
#ax1.set_ylim([40, 200])
#ax1.legend([np.mean(Si_r)])
#
#ax2=plt.subplot(2,2,2)
#ax2.scatter(x, G_r,c='r')
#plt.plot(x,pm(x,y,G_r),c='black')
#ax2.set_title('Graphene_rise')
#ax2.set_ylim([40, 200])
#ax2.legend([np.mean(G_r)])
#
#ax3=plt.subplot(2,2,3)
#ax3.scatter(x, Si_d,c='r')
#plt.plot(x,pm(x,y,Si_d),c='black')
#ax3.set_title('Silicon_decay')
#ax3.set_ylim([100, 220])
#ax3.legend([np.mean(Si_d)])
#
#ax4=plt.subplot(2,2,4)
#ax4.scatter(x, G_d,c='r')
#plt.plot(x,pm(x,y,G_d),c='black')
#ax4.set_title('Graphene_decay')
#ax4.set_ylim([100, 220])
#ax4.legend([np.mean(G_d)])


#####ALL IN ONE########



y1=pm(x,y,Si_r)
y2=pm(x,y,G_r)
y3=pm(x,y,Si_d)
y4=pm(x,y,G_d)




plt.figure(figsize=(10,10))
plt.ylabel('Time_Constant [miliseconds]')
plt.xlabel('Experiment Number \n \n (Horizontal lines are the mean value of the corresponding data set)')
           

plt.scatter(x, Si_r,c='r',label='Silicon Rise')
plt.hold('on')
plt.scatter(x, G_r,c='b',label='Graphene Rise')
plt.scatter(x, Si_d,c='g',label='Silicon Decay')
plt.scatter(x, G_d,c='black',label='Graphene Decay')
plt.plot(x,y1,c='r')
plt.plot(x,y2,c='b')
plt.plot(x,y3,c='g')
plt.plot(x,y4,c='black')
plt.legend()

