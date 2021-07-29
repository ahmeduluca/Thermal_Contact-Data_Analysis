a= open("tip02.txt", "r")
b= [ ]
c= [ ]
for j in range (0,50):

    for i in range(j,331200,50):
        b.append(float(a.readline()))
    c.append(sum(b)/len(b))
print(c)
n=open("sample2_aver2.txt", "w")
for k in range(0,50):
    n.write(str(c[k]))
n.close()
