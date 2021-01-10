# Mean Residence Time and E Curve

t = [0, 5, 10, 15, 20, 25, 30, 35]
c = [0, 3, 5, 5, 4, 2, 1, 0]
dtI = 5
totalNum = 0
totalDeno = 0
ECurve = []

def NumDeno():
    global totalDeno, totalNum
    for i in range(len(t)):
        tINum = t[i] * c[i] * dtI
        totalNum += tINum
        tIDeno = c[i] * dtI
        totalDeno += tIDeno

def eCurve():
    global ECurve
    for i in range(len(c)):
        e = c[i]/totalDeno
        ECurve.append(e)
    

NumDeno()
eCurve()
print(totalNum/totalDeno, totalDeno)
print(ECurve)


'''
tI0 = (t[0]*c[0]*dtI)#/(c[0]*dtI)
tI1 = (t[1]*c[1]*dtI)#/(c[1]*dtI)
tI2 = (t[2]*c[2]*dtI)#/(c[2]*dtI)
tI3 = (t[3]*c[3]*dtI)#/(c[3]*dtI)
tI4 = (t[4]*c[4]*dtI)#/(c[4]*dtI)
tI5 = (t[5]*c[5]*dtI)#/(c[5]*dtI)
totaltI = tI0 + tI1 + tI2 + tI3 + tI4 + tI5

TI0 = (c[0]*dtI)
TI1 = (c[1]*dtI)
TI2 = (c[2]*dtI)
TI3 = (c[3]*dtI)
TI4 = (c[4]*dtI)
TI5 = (c[5]*dtI)
TITotal = TI0 + TI1 + TI2 + TI3 + TI4 + TI5

total = totaltI/TITotal

print(totaltI)
'''
