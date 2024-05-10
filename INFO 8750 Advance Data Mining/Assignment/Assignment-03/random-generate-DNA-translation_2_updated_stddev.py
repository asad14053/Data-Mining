# @Copyright-
#     Compiled by: Md Asaduzzaman Jabin
#     Ph.D. Student, GRA, WaveLab
#     School of ECE, UGA
# For class homework 2 Info 8750

#Task-01
import random  
import string
import sys
import math
import statistics

# To print all possible combination and write it to a file
#sys.stdout = open('output_code1.txt', 'w')

# mean calculation
def mean(a, n):
    sm = 0
    for i in range(n):
        sm+=a[i]
    return sm/n


# manual stardard deviation calculation
def stdev(m, n):
    sum_x = 0
    sum_x2 = 0
    for i in range(n):
        sum_x += m[i]
        sum_x2 += m[i]**2
    stddev = (sum_x2/n-(sum_x/n)**2)**0.5
    return stddev

#DNA translation Alogorithm
def random_dna():

    # To print all possible combination and write it to a file
    #sys.stdout = open('output_code1.txt', 'w')

    # Algorithm
    # 1) if we notice to the DNA translation list below , we will find 64 unique DNA unit sequence
    # 2) generate random number between 1-1000 and map it through 1 - 64 to applly this dictionaries to get your protien sequence
    # 3) map the sequence and check 'ATG' as starting, 'TAA', 'TAG' or, 'TGA' for ending sequence
    # 4) find average length and standard deviation of protien sequence

    # Translation dictionaries
    ############################################### Generating Sequence ##################################################
    dna2aa = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "stop", "TAG": "stop",
        "TGT": "C", "TGC": "C", "TGA": "stop", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
        
        # mapped random number into sequence using dictionaries

    ran2dna = ["TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT",
               "TGC", "TGA", "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC",
                "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA",
                "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG",
                "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"]
    
    comb_std = []
    comb_mean = []
    
    # Generate 1000 DNA sequence
    for k in range (1000):
            
        # Task -02: randomly generate 1000 lenth DNA sequence but, we are considering 3 sequence at a time, so loop will continue 333*3 --> 999+1 = 1000
        full_seq = ''
        full_seq += random.choice(['A','G','T', 'C']) #just generate 1 nuclic acid
        for i in range (1,333):
            n = random.randint(1, 1001)
            #print(n)

            # Task-03: since we have 64 sequence to map from random number
            modulo = int(n % 64)
            full_seq = full_seq + (ran2dna[modulo]).upper() #convert into uppercase

        print("Random Sequence No: ",k+1," - ", full_seq)

        ################################################### Translate the sequence #####################################

        #virtual translation
        a  = full_seq
        l = len(a)
        pro = ''
        seq =''
        no_pro = 0
        stop_coden = ['TAA','TAG','TGA']
        start_pos = 0
        max_len = 0
        max_seq = ''
        avg_len = []
        len_list =[]
        ln  = 0

        
        for j in range(l-3): # iterate untill the end of the length 
            c3 = a[j]+a[j+1]+a[j+2]

            if c3 == 'ATG':   #for us to locate where all 'ATG' are
                start_pos = j
            
            # Task-04 Start Iteratig while we found a ATG 
            for i in range(start_pos,l-3,3): #Start Iteratig while we found a ATG 
                c3 = a[i]+a[i+1]+a[i+2]

                #if c3 == 'TAA' or c3 == 'TAG' or c3 == 'TGA':
                if c3 in stop_coden:
                    #print(c3)
                    break
                else:
                    pro+= dna2aa[c3]
                    seq+=c3

            # Task-04 Find max sequence
            ln = len(pro)
            if(max_len<ln):
                max_len = ln
                max_seq = pro

            # Task-05 Find average Protein length 
            avg_len.append(ln)
            len_list.append(ln)
            #print(avg_len)

            no_pro += 1
            # print("Number of protien: ", no_pro)
            # print("DNA sequence: ", seq)
            # print("Protein length: ", len(pro))
            # print("Protien Sequence: ", pro)
            # print()
            pro = ''
            seq = ''

        print()
        print("Max Protien length: ", max_len)
        print("Max Protien sequence: ", max_seq)

        # Task-05 Find average Protein length
        mn = math.ceil(mean(avg_len, no_pro))
        print("Average Protien length: ", mn)

        # Task-06 standard deviation of protein length
        st = math.ceil(stdev(len_list, no_pro))
        print("Protein length Standard Deviation: ", st)


        # store each local mean and std to find total 1000 DNA sequence mean and standard sequence 
        comb_mean.append(mn)
        comb_std.append(st)
    
    print()
    # Task-07 Find global average Protein length
    mn = math.ceil(mean(comb_mean, 1000))
    print("Global Average Protien length: ", mn)

    # Task-08 global standard deviation of protein length
    st = math.ceil(stdev(comb_std, 1000))
    print("Global Protein length Standard Deviation: ", st)


random_dna()