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
"""
This code is for reading experimental data of Indenter/RC Thermal setup obtained by
NI-DAQ written to TDMS-- there are also txt data (suchas Loadcell 24BitADC HX711)
So,'Load' and 'Displacement' files are being searched also as txt. 
"""
bit_to_load=4.75*9.81
loadcell_stiffness=536.5
load_files=[]
position_files=[]
lvsdis_file=[[]]
load=[[]]
position=[[]]
date=[[]]
n=0
tdms_file_list=glob.glob("/Users/ahmeduluca/Desktop/Al/11-06-2021/**/*.tdms")
txt_file_list=glob.glob("/Users/ahmeduluca/Desktop/2020-cals/approach_cal/fused_silica/*.txt")
##for i in txt_file_list:
##    if("Load" in i):
##        date.append([])
##        load_files.append(i)
##        date[n].append(i[-29:-9])
##        n+=1
##    elif("Position" in i):
##        position_files.append(i)
##n=0        
##for k in range(len(load_files)):
##    load.append([])
##    position.append([])
##    f=open(load_files[k],"r")
##    for x in f:
##        load[k].append(float(x)*bit_to_load/1000)
##    for j in position_files:
##        if(len(date[k])!=0):
##            if(date[k][0] in j):
##                g=open(j,"r")
##                n=0
##                for x in g:
##                    position[k].append((float(x)-load[k][n]/loadcell_stiffness)*1000)
## ##                   position[k].append((float(x))*1000)
##                    n+=1
##                plt.plot(position[k],load[k], label=date[k])
##                continue
##plt.xlabel('Displacement (nm)')
##plt.ylabel('Force (mN)')
##plt.legend()
##plt.show()
sayg=0
sayc=0
freq=[]
spect=[]
cut_sign=[]
tdms_file_list.sort()
for i in tdms_file_list:
    try:
        with TdmsFile.open(i) as tdms0:
            sayg=(len(tdms0.groups()))
            sayc=0
            for group in tdms0.groups():
                sayc=len(group.channels()) if len(group.channels())>sayc else sayc
            fig,axs=plt.subplots(sayc,sayg)
            fig.subplots_adjust(hspace=0.1)
            fig.subplots_adjust(wspace=0.5)
            sayg=0
            for group in tdms0.groups():
                sayc=0
                for channel in group.channels():
                    times=channel.time_track()
                    t_new=np.arange(times[0],times.max(),max(times)/len(times))
                    try:
                        if(len(tdms0.groups())==1 and len(group.channels())>1):
                            axs[sayc].set_ylabel(channel.name)
                            axs[sayc].plot(times,channel,'b')
                            axs[0].set_title(group.name)
##                            tck=interp.splrep(times, channel)
                            if (sayc==-1):
                                fft_signal=rfft(channel[:])
                                freq=rfftfreq(len(channel[:]), d=times[1]-times[0])
##                                if(freq[0]==0):
##                                    fft_signal=fft_signal[1:]
##                                    freq=freq[1:]
                                print(freq)
                                spect=freq**2
                                cut_f_sign=fft_signal.copy()
                                maxf=np.where(abs(cut_f_sign)==max(abs(cut_f_sign)))
                                print(cut_f_sign)
                                cut_f_sign[(spect>spect.max()*(1-np.exp(-100/len(channel[:]))))]=0
                                cut_sign=irfft(cut_f_sign)
                                max_sign=np.sin(2*times*np.abs(freq[maxf[0]-1])*np.pi)+np.sin(2*times*np.abs(freq[maxf[0]])*np.pi)+np.sin(2*times*np.abs(freq[maxf[0]+1])*np.pi)
#                                f00=savgol_filter(channel, 9,5)
                                axs[sayc].plot(times,channel,'b',times,cut_sign,'r--',times,max_sign,'y*')
                                axs[sayc+1].plot(freq,np.abs(cut_f_sign),'ro',freq,np.abs(fft_signal),'b')
                                print(freq[maxf]*times[-1])
                                print(times[-1])
                                print(freq[maxf],maxf.index)
                                axs[0].set_title(group.name)
                            elif(sayc==len(group.channels())-1):
                                axs[sayc].set_xlabel('Time')
                            sayc+=1
                        elif(len(tdms0.groups())==1 and len(group.channels())==1):
                            axs.set_ylabel(channel.name)
                            tck=interp.splrep(times, channel)
                            f00=interp.splev(t_new, tck)
                            axs.plot(t_new,f00,'r')
                            axs.set_xlabel('Time')
                            axs.set_title(group.name)
                            sayc+=1
                        elif(len(tdms0.groups())>1 and len(group.channels())==1):
                            axs[sayg].set_ylabel(channel.name)
                            tck=interp.splrep(times, channel)
                            f00=interp.splev(t_new, tck)
                            axs[sayg].plot(times,channel,'r')
                            axs[sayg].set_title(group.name)
                            axs[sayg].set_xlabel('Time')
                            sayc+=1
                        else:
                            axs[sayc,sayg].set_ylabel(channel.name)
                            tck=interp.splrep(times, channel)
                            f00=interp.splev(t_new, tck)
                            axs[sayc,sayg].plot(t_new,f00,'r')
                            if (sayc==0):
                                axs[0,sayg].set_title(group.name)
                            elif(sayc==len(group.channels())-1):
                                axs[sayc,sayg].set_xlabel('Time')
                            sayc+=1
                    except KeyError as err:
                        print(err," : "+group.name+" in "+tdms0.properties["name"])
                        continue
                    except IndexError as err:
                        print(err," : "+group.name+" in "+tdms0.properties["name"])
                        continue
                    except:
                        print("Some error in plotting : " + tdms0.properties["name"])
                        continue
                sayg+=1
            fig.suptitle(os.path.basename(os.path.dirname(i)))
            plt.show()
            tdms0.close()
            sayc=0
            sayg=0
    except KeyError as err:
        print(err," : ",i)
        continue
    except:
        print('some error in loop, file:',i)
        continue
#        break
##        plt.show()
