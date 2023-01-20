import pandas as pd
import math
import numpy as np
df = pd.read_csv ('Vipvalues.csv')
Seq =[]

def Reverse(lst):
    new_lst = lst[::-1]
    return new_lst


data = np.genfromtxt("Vipvalues.csv", dtype= str, encoding = None, delimiter=",")


Seq=np.transpose(data)[0]
vIP=np.transpose(data)[1]


print('Please provide a non empty nucleotide sequence:')
x = input()


print('Please choose: S (single strand) or D (double strand)')
z = input()

if z!= 'S' and z!= 'D':
    print("Please chose S or D")
    z = input() 	
    if z!= 'S' and z!='D':    
        print("Compute single strand vIP")  
   
INPUT=list(x)

REVINPUT=Reverse(INPUT)


for i in range(len(REVINPUT)):
  
    # replace A
    if REVINPUT[i] == 'A':
        REVINPUT[i] = 't'
  
    # replace T
    if REVINPUT[i] == 'T':
        REVINPUT[i] = 'a'

    # replace C
    if REVINPUT[i] == 'C':
        REVINPUT[i] = 'g'
  
    # replace G
    if REVINPUT[i] == 'G':
        REVINPUT[i] = 'c'


    # replace M
    if REVINPUT[i] == 'M':
        REVINPUT[i] = 'g'


    REVINPUT[i] = REVINPUT[i].upper()

y=''.join(REVINPUT)


lenin=len(INPUT)
a=0;
for SEI in INPUT:
    if SEI != 'A' and SEI != 'C' and SEI != 'G' and SEI != 'T' and SEI != 'M' and SEI != '\n':
        print("Please provide a sequence with only standard characters for nucleotides : C (Cytosine), G (Guanine), A (Adenine), T (Thymine), M (5-Methylcytosine))");
        a=1;
        break; 

if z=='S':
    if a == 0:
        if lenin < 6:
           idx =  np.where(Seq == x)
           if(len(idx[0])>0):
               print('The vIP value of the single stranded DNA "',x,'" is',round(float(vIP[idx[0][0]]),3),'eV')
        if (lenin>=6):
           CECCA=lenin-3;
           POLLO = ["" for l in range(CECCA)]
           RISULTATO=np.zeros(CECCA,float)
           FACTOR=np.zeros(CECCA,float)       
           for i in range(CECCA):
               POLLO[i]=x[i:i+4]
               fame=np.where(Seq == POLLO[i])          
               FACTOR[i]=pow(lenin-i,2)/(pow(lenin-i,2)-2) 
               RISULTATO[i]=math.factorial(CECCA-1)/(math.factorial(i)*math.factorial(CECCA-i-1))*pow(0.4009940,CECCA-1-i)*pow(1-0.4009940,i)*float(vIP[fame[0][0]])
           FINAL=8.34424*np.exp(0.56074*(pow(lenin, -0.145481)-1))+np.prod(FACTOR)*np.sum(RISULTATO)-np.prod(FACTOR)*8.34424*np.exp(0.56074*(pow(lenin-CECCA+1, -0.145481)-1)) 
           print('The vIP value of the single stranded DNA "',x,'" is',round(FINAL,3),'eV')


if z=='D':
    if a == 0:
        if lenin < 6:
           idx =  np.where(Seq == x)
           idy =  np.where(Seq == y)
           if(len(idx[0])>0):
               print('The vIP value of the double stranded DNA "',x,'" is', round((float(vIP[idx[0][0]])+float(vIP[idy[0][0]]))/2,3),'eV')
        if (lenin>=6):
           CECCA=lenin-3;
           POLLO = ["" for l in range(CECCA)]
           RISULTATO=np.zeros(CECCA,float)
           FACTOR=np.zeros(CECCA,float)
           for i in range(CECCA):
               POLLO[i]=x[i:i+4]
               fame=np.where(Seq == POLLO[i])
               FACTOR[i]=pow(lenin-i,2)/(pow(lenin-i,2)-2)
               RISULTATO[i]=math.factorial(CECCA-1)/(math.factorial(i)*math.factorial(CECCA-i-1))*pow(0.4009940,CECCA-1-i)*pow(1-0.4009940,i)*float(vIP[fame[0][0]])
           FINAL1=8.34424*np.exp(0.56074*(pow(lenin, -0.145481)-1))+np.prod(FACTOR)*np.sum(RISULTATO)-np.prod(FACTOR)*8.34424*np.exp(0.56074*(pow(lenin-CECCA+1, -0.145481)-1))
           for i in range(CECCA):
               POLLO[i]=y[i:i+4]
               fame=np.where(Seq == POLLO[i])
               FACTOR[i]=pow(lenin-i,2)/(pow(lenin-i,2)-2)
               RISULTATO[i]=math.factorial(CECCA-1)/(math.factorial(i)*math.factorial(CECCA-i-1))*pow(0.4009940,CECCA-1-i)*pow(1-0.4009940,i)*float(vIP[fame[0][0]])
           FINAL2=8.34424*np.exp(0.56074*(pow(lenin, -0.145481)-1))+np.prod(FACTOR)*np.sum(RISULTATO)-np.prod(FACTOR)*8.34424*np.exp(0.56074*(pow(lenin-CECCA+1, -0.145481)-1))
           print('The vIP value of the double stranded DNA "',x,'" is',round((FINAL1+FINAL2)/2,3),'eV') 
