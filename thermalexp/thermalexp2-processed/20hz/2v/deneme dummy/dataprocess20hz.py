#import csv
#import numpy as np
#import matplotlib.pyplot as plt

c=[ ]
tip=[ ]
data=open('dummy.txt', 'r')
lines= data.readlines()
for x in lines:
    tip.append(x.replace(',', '.'))
for j in range (0,90):
    b=[ ]
    for i in range(j,len(tip),90):
        b.append(float(tip[i]))
    c.append(sum(b)/len(b))
print(c)
with open("dummy_aver2p2.txt", "w") as n:
    for k in range(len(c)):
        n.write(str(c[k]))
        n.write ("                                                                        ")
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


