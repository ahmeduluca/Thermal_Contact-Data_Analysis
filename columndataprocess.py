import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import glob
import pandas as pd
import matplotlib.gridspec as gridspec
import scipy.interpolate as interp
from scipy.signal import savgol_filter
from scipy.fftpack import rfft, irfft, rfftfreq

file1="/Users/ahmeduluca/Desktop/CRANN-INTERNSHIP/thermal_tabletop/"
file2="/Users/ahmeduluca/Desktop/püton/thesis/thermalexp/"
#file_list=glob.glob("/Users/ahmeduluca/Desktop/CRANN-INTERNSHIP/thermal_tabletop/thermal3/*.txt")
#txt_file_list=glob.glob("/Volumes/AhmedUluca/YEDEK_24-04-2020/CEM_GRAPHENE/IFP/CEM_GRAPHENE/**/*.txt")
txt_file_list=glob.glob("/Users/ahmeduluca/Desktop/Al/Process/**/**/*.txt")
dir_list=glob.glob("/Users/ahmeduluca/Desktop/Al/Process/*/")
step_list=glob.glob("/Users/ahmeduluca/Desktop/Al/Process/**/*/")
txt_file_list.sort()
    
##i=0
##def get_col(arr, col):
##    return np.transpose(arr)[col]
####unpack=True ile get_col olmadan da çalışılabilir
##for k in file_list:
##    try:
##        data=np.genfromtxt(k,delimiter="\t")
##        if(data.ndim<2):
##           data=np.genfromtxt(k)
##        if(data.ndim>1 and get_col(data,0)[0]!="nan"):
##            for i in range(data.shape[1]):
##                dt=get_col(data,i)
##                plt.plot(dt,".-",label="col"+str(i))
##        elif(data.ndim==1 and data.size>0 and data[0]!="nan"):
##            plt.plot(data,".-",label="tek")
##            print(data[0])
##        else:
##            data=np.genfromtxt(k)##pandas vs deneme ',' ondalık  gibi yazımlar için
##            plt.plot(data,".-",label="lambda")
##
##        plt.title(k.replace(file2,""))
##        plt.legend()
##        plt.show()
##    except:
##        print("error in : ", k)
##        continue

##for k in txt_file_list:
##    try:
##        data0,data1,data2,data3 = np.genfromtxt(k,delimiter="\t", unpack=True)
##        plt.plot(data0,data1,data0,data2,data0,data3)
##        plt.title(k.replace(file2,""))
##        plt.legend()
##        plt.show()
##    except:
##        print("error in : ", k)
##        continue

for k in txt_file_list:
    
    try:
        data0,data1 = np.genfromtxt(k, unpack=True)
        plt.plot(data0,data1)
 #       plt.title(str(os.path.basename(k)))
        plt.legend()
        plt.show()
        fft_signal=rfft(data1)
        freq=rfftfreq(len(data1), d=data0[1]-data0[0])
        if(freq[0]==0):
            fft_signal=fft_signal[1:]
            freq=freq[1:]
        spect=freq**2
        plt.plot(spect,fft_signal)
        plt.show()
        cut_f_sign=fft_signal.copy()
        maxf=np.where(abs(cut_f_sign)==max(abs(cut_f_sign)))
        print(maxf[0])
        cut_f_sign[(spect>100)]=0
        cut_sign=irfft(cut_f_sign)
        plt.plot(data0, cut_sign)
        plt.show()
        max_sign=np.sin(2*data0*np.abs(freq[maxf[0]-1])*np.pi)+np.sin(2*data0*np.abs(freq[maxf[0]])*np.pi)+np.sin(2*data0*np.abs(freq[maxf[0]+1])*np.pi)
        plt.plot(data0, max_sign)
        plt.show()
#
    except:
        print("error in : ", k)
        continue
    
