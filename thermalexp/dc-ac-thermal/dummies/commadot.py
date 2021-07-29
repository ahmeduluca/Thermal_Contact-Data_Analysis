with open('acthermaltest2_40hz_1600mvac_46mvst.txt', 'r') as data:
    plaintext = data.read()
    plaintext = plaintext.replace(',', '.')
n=open("acthermaltest2_40hz_1600mvac_46mvst.txt", "w")    
n.write(plaintext)
n.close
