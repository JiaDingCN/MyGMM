#GMMResult是一个纯文本文件，里面放的是Main.cpp运行后输出的分类结果，具体见Sample
dataFile=open("/home/jiading/dataGMM.txt")
labelFile=open("/home/jiading/Desktop/GMMResult")
X=[]
Y=[]
label=[]
dataNumber=0
for line in dataFile.readlines():
    tempVector=line.split(" ")
    if(len(tempVector[0])==0):
        X.append(float(tempVector[2]))
        Y.append(float(tempVector[3]))
    else:
        X.append(float(tempVector[1]))
        Y.append(float(tempVector[2]))
    dataNumber+=1
for line in labelFile.readlines():
    label.append(int(line))
import matplotlib.pyplot as plt
for i in range(dataNumber):
    if(label[i]==0):
        plt.scatter(X[i],Y[i],c='b')
    else:
        plt.scatter(X[i],Y[i],c='g')
plt.show()