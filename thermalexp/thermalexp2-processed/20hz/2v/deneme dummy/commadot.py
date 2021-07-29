y=[ ]
data=open('dummy.txt', 'r')
lines= data.readlines()
for x in lines:
    y.append(x.replace(',', '.'))
print(y)
