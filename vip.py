import pandas as pd
import math
import numpy as np
df = pd.read_csv ('Vipvalues.csv')
Seq =[]

data = np.genfromtxt("Vipvalues.csv", dtype= str, encoding = None, delimiter=",")


Seq=np.transpose(data)[0]
vIP=np.transpose(data)[1]


print('Please provide a non empty nucleotide sequence:')
x = input()
INPUT=list(x)
lenin=len(INPUT)
a=0;
for SEI in INPUT:
    if SEI != 'A' and SEI != 'C' and SEI != 'G' and SEI != 'T' and SEI != 'M' and SEI != '\n':
        print("Please provide a sequence with only standard characters for nucleotides : C (Cytosine), G (Guanine), A (Adenine), T (Thymine), M (5-Methylcytosine))");
        a=1;
        break; 

if a == 0:
    if lenin < 7:
       idx =  np.where(Seq == x)
       if(len(idx[0])>0):
           print('The vIP value of "',x,'"is',vIP[idx[0][0]])
       if(len(idx[0]) == 0):
           DL=8.34424*np.exp(0.56074*(pow(lenin, -0.145481)-1))
           if (lenin==5):
               DLMENO1=8.34424*np.exp(0.56074*(pow(lenin-1, -0.145481)-1))
               primo=x[:-1]
               primovip=np.where(Seq == primo)  
               UNO=vIP[primovip[0][0]]
               second=x[1:]
               seconvip=np.where(Seq == second)
               DUE=vIP[seconvip[0][0]]
               PP = DL+(0.3620744553889976 * (float(UNO)-DLMENO1)+0.588373793863788*(float(DUE)-DLMENO1))
               print(PP)
           if (lenin==6):
               DLMENO1=8.34424*np.exp(0.56074*(pow(lenin-1, -0.145481)-1))
               DLMENO2=8.34424*np.exp(0.56074*(pow(lenin-2, -0.145481)-1))
               primo=x[:-2]
               primovip=np.where(Seq == primo)
               UNO=vIP[primovip[0][0]]
               second=x[1:-1]
               seconvip=np.where(Seq == second)
               DUE=vIP[seconvip[0][0]]
               terzo=x[2:]
               terzovip=np.where(Seq == terzo)
               TRE=vIP[terzovip[0][0]]
               PP = DL + 0.3620744553889976 * 0.3620744553889976 * float(UNO) + 0.588373793863788*0.588373793863788*float(TRE)+2*0.3620744553889976 * 0.588373793863788 * float(DUE)-(0.3620744553889976+0.588373793863788)*(0.3620744553889976+0.588373793863788)*DLMENO2
               print('The vIP value of "',x,'" is',PP)
    if (lenin>=7):
       CECCA=lenin-3;
       POLLO = ["" for l in range(CECCA)]
       RISULTATO=np.zeros(CECCA,float)
       for i in range(CECCA):
           POLLO[i]=x[i:i+4]
           fame=np.where(Seq == POLLO[i])           
           RISULTATO[i]=math.factorial(CECCA-1)/(math.factorial(i)*math.factorial(CECCA-i-1))*pow(0.3620744553889976,CECCA-1-i)*pow(0.588373793863788,i)*float(vIP[fame[0][0]])
       FINAL=8.34424*np.exp(0.56074*(pow(lenin, -0.145481)-1))+np.sum(RISULTATO)-pow(0.950448,CECCA-1)*8.34424*np.exp(0.56074*(pow(lenin-CECCA+1, -0.145481)-1)) 
       print('The vIP value of "',x,'" is',FINAL) 
