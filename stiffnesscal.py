import tkinter as tk
import tkinter.messagebox as messagetk
import numpy
import scipy
import matplotlib.pyplot as plt
from numpy import pi, linspace, inf, array
top=tk.Tk()
top.title("Stiffness  Calculator")
def takemasval(Event):
    try:
       mass=float(Event.widget.get())*0.001
       calc[0]=mass
       mc=1
       a[0]=1
       cal(calc[0],calc[1],calc[2],a[0],a[1],a[2])
    except ValueError as err:
        messagetk.showerror("Value Error", err)
def takecntval(Event):
    try:
       count=float(Event.widget.get())
       calc[1]=count
       cc=1
       a[1]=1
       cal(calc[0],calc[1],calc[2],a[0],a[1],a[2])
    except ValueError as err:
        messagetk.showerror("Value Error", err)
def takestfval(Event):
    try:
       stiff=float(Event.widget.get())
       calc[2]=stiff
       sc=1
       a[2]=1
       cal(calc[0],calc[1],calc[2],a[0],a[1],a[2])
    except ValueError as err:
        messagetk.showerror("Value Error", err)
def cal(mass, cnt, stiff, mc, cc, sc):
    if(sc==1 and cc==1):
        a[0]=0
        a[1]=0
        a[2]=0
        mass=1000*stiff*(cnt*dx)/g
        E0.delete(0,tk.END)
        E0.insert(0,str(mass))
        return stiff*(cnt*dx)/g
    elif(sc==1 and mc==1):
        sc=0
        mc=0
        a[0]=0
        a[1]=0
        a[2]=0
        count=mass*g/(stiff*dx)
        A1.delete(0,tk.END)
        A1.insert(0,str(count))
        return mass*g/(stiff*dx)
    elif(mc==1 and cc==1):
        mc=0
        cc=0
        a[0]=0
        a[1]=0
        a[2]=0
        stiff=mass*g/(cnt*dx)
        R1.delete(0,tk.END)
        R1.insert(0,str(stiff))
        return mass*g/(cnt*dx)
    else:
        return 0
g=9.81
count=0.
dx=655.e-9
mass=0.
stiff=0.
a=[0,0,0]
calc=[0.,0.,0.]
L0=tk.Label(top,text="Mass(g):",width=8, height=2)
L0.pack(fill=tk.Y,side=tk.TOP,anchor='w')
E0=tk.Entry(top,bd=4,width=4)
E0.bind('<Return>',takemasval)
E0.pack(fill=tk.X,side=tk.TOP,anchor='w')
A0=tk.Label(top,text="Pattern Counts:",width=12, height=2)
A0.pack(fill=tk.Y,side=tk.TOP,anchor='w')
A1=tk.Entry(top,bd=4,width=5)
A1.pack(fill=tk.X,side=tk.TOP,anchor='w')
A1.bind('<Return>',takecntval)
R0=tk.Label(top,text="Stiffness (N/m):",width=10, height=2)
R0.pack(fill=tk.Y,side=tk.TOP,anchor='w')
R1=tk.Entry(top,bd=4,width=5)
R1.pack(fill=tk.X,side=tk.TOP,anchor='w')
R1.bind('<Return>',takestfval)
top.minsize(width=500,height=200)
top.mainloop()

      
    
