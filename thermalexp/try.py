a= open("termoc1.txt", "r")
b= [ ]
c= [ ]
for j in range (0,10):

    for i in range(j,13800,10):
        b.append(float(a.readline()))
    c.append(sum(b)/len(b))
print(c)
n=open("thermocouple_aver1.txt", "w")
for k in range(0,10):
    n.write(str(c[k]).rjust)
n.close()

