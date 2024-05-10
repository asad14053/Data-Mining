import random
import string
import sys

# To print all possible combination and write it to a file
sys.stdout = open('output_code1.txt', 'w')

def decoder():
    secret1 = ("cldgiJB$3sb<)?se?gx|vD2d[sP\>1)J8k`us-U\Hl{nBsJc!pBf-HS_*7Jv][*/Bl;&.Kv.V5uoRx9!9kIIev]wF<Z]B~Xe~Ml*81OxTs,Y>U9B0B3_X=B(SV8_2y#0f4bkd\
C$op)W&h6_>CushH9")

    # The possible 1st secret messages are:
         # Jesus_loves_you  key1:  10 key2:  5

    secret2 = ("3_\iW*JHc`(MTBhI3-\
aQMxrKaJh03~L^vs?ZXZ=26p]B+*RCf1m;<kp>^*px$v?$69a&;mg:kRKyi4JRE:86/DNjn)(9.y7CBP>7`%N0.ll\
HzT@|:~y\#5HPK?^DL1iRluL_uJmE#p+dm{\gY-\
'Sx67nXc.kd371e^PsRRrJ{0a}@@*$<rzU[jTj_d^9h(dI6J.^{_w+WiW`Q&Bz5ApY'|BmXsSq-.q;WO,7eApWc\
%BXQ(k$`s\oVa[5S9+j=?M'-<d6]irtwv0b&?16")
    
    # The possible 2nd secret messages are:
        # Happy.Tiger.Year  key1:  17 key2:  7

        # Rai  key1:  91 key2:  79
        # hlw  key1:  92 key2:  14
        # W0W  key1:  99 key2:  4
        # QLA  key1:  105 key2:  19
        # Mic  key1:  106 key2:  20
        # 2k6  key1:  116 key2:  39
        # Cs  key1:  116 key2:  47
        # mr  key1:  116 key2:  50

    # choose input between secret1 or, secret2
    secret = secret1
    l = len(secret)
    for key1 in range(1, l):
       
        for key2 in range(1, key1):
            real=""

            # for each key1 and key2 check full encrypted message 
            for i in range(l):
                if i % key1 == key2:
                # message1: key1 = 10, key2 = 5
                # message2: key1 = 17, key2 = 7
                    real+= secret[i]
            print(real," key1: ", key1, "key2: ", key2)

decoder()