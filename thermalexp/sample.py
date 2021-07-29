a= open("aluminum2.txt", "r")
b= [ ]
c= [ ]
for j in range (0,10):
    i=j+1
    for i in range(0,13800,10):
        b.append(float(a.readline()))
    c.append(sum(b)/len(b))
print(c)
n=open("aluminum_aver.txt", "w")
for k in range(0,10):
    n.write(str(c[k]))
n.close()
