import csv
#import numpy as np
#import matplotlib.pyplot as plt
b=[ ]
c=[ ]
d=[ ]
e=[ ]
f=[ ]
g=[ ]
tip=[ ]
sample=[ ]
medium=[ ]
with open("acthermaltest20hz1Vac.txt") as a:
    lines=a.readlines()
    for x in lines:
        tip.append(x.split()[1])
        sample.append(x.split()[2])
        medium.append(x.split()[3])
    for j in range (0,100):
        for i in range(j,len(tip),100):
            b.append(float(tip[i]))
            d.append(float(sample[i]))
            e.append(float(medium[i]))
        c.append(sum(b)/len(b))
        f.append(sum(d)/len(d))
        g.append(sum(e)/len(e))
print(c)
print(f)
print(g)
with open("tip_aver_20hz1v.txt", "w") as n:
    for k in range(len(c)):
        n.write(str(c[k]))
with open("sample_aver_20hz1v.txt", "w") as m:
    for l in range(len(c)):
        m.write(str(c[k]))
with open("medium_aver_20hz1v.txt", "w") as o:
    for p in range(len(c)):
        o.write(str(c[k]))
#table=[[ ],[ ],[ ]]
#for l in range(len(c)):
#    for m in range(len(c)):
#        for o in range(len(c)):
#            table.append(g[o])
#        table.append(f[m])
#    table.append(c[l])
#print (table)
#with open("test1.csv", "w") as csvfile:
#               writer= csv.writer(csvfile)
#               [writer.writerow(r) for r in table]
#count=np.arange(0,len(c))
#plt.plot(c,count,"r--",f,count,"bs",g,count,"g")


