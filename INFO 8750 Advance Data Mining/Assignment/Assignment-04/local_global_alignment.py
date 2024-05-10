# @Copyright-
#     Compiled by: Md Asaduzzaman Jabin
#     Ph.D. Student, GRA, WaveLab
#     School of ECE, UGA
# For class homework 3 Info 8750

import random  
import string
import sys
import math
import statistics

#import installed packages 
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats 
from numpy.random import seed
from numpy.random import randn
from numpy.random import normal
from scipy.stats import ttest_1samp

# match bonus
match = 3
# mismatch penalty
mismatch = -3
#gap penalty
gap = -2 

#store global scores
st_global = []
#store local scores
st_local = []

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

#Mapped random number into sequence using dictionaries --> 64 Sequence
ran2dna = ["TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT",
            "TGC", "TGA", "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC",
            "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA",
            "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG",
            "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"]

#Write it to a file
#sys.stdout = open('output_code1.txt', 'w')

#Mean calculation
def mean(a, n):
    sm = 0
    for i in range(n):
        sm+=a[i]
    return sm/n


#Manual stardard deviation calculation
def stdev(m, n):
    sum_x = 0
    sum_x2 = 0
    for i in range(n):
        sum_x += m[i]
        sum_x2 += m[i]**2
    stddev = (sum_x2/n-(sum_x/n)**2)**0.5
    return stddev

#Define function for random DNA sequence of length 'N'
#Function definition:
##Generate_random1000(len_of_sequence)
def gen_ran1000(N):
            
    #randomly generate 1000 lenth DNA sequence but, we are considering 3 sequence at a time, so loop will continue 333*3 --> 999+1 = 1000
    full_seq = ''
    N = int(N)
    mod = int(N%3)
    if(mod):
        while(mod):
            full_seq += random.choice(['A','G','T', 'C']) #just generate few nuclic acid that does not fit with mod 3
            mod-=1
    for i in range (1,int(N/3.0)):
        n = random.randint(1, 1000)
        #print(n)

        #since we have 64 sequence to map from random number
        modulo = int(n % 64)
        full_seq = full_seq + (ran2dna[modulo]).upper() #convert into uppercase

    return full_seq


#Define function for T-test
#Function definition:
## T-test(data_sample, mean, standard_deviation, name_of_alignment)
def T_test(sample, avg, std, name):
    
    if(name == 'Global'):
        print ('-------------------------------------------T-test for Global alignment-----------------------------------------------------')
    else:
        print ('-------------------------------------------T-test for Local alignment-----------------------------------------------------')


    # T-statistic close to zero represents the lowest evidence against the hypothesis
    # P-value is the percentage probability of the t-statistic to have occurred by chance
    # alpha = ( p-value /2 ) always should be > 0.05, A significance level is the percentage probability of rejecting a true null hypothesis. This is also called alpha
    # Null hypothesis:  H0 : sample mean <= hypothetical mean
    # Reject the null hypothesis if p-value <= alpha

    t_stat, p_value = ttest_1samp(sample, popmean= avg+10)    # pop_mean >= sample mean
    alpha = 0.05
    print("T-statistic value: ", t_stat)  
    print("P-Value: ", p_value)

    print("Null hypothesis: Sample mean <= Hypothetical mean")
    print()
    print("Significance Decision: ")

    if(p_value<= alpha):
        print ("Reject the hypothesis")
    else:
        print("Failed to reject the hypothesis")
    
    if(t_stat<0):
        print("Direction of the sample mean extreme")
    else:
        print("Direction of the sample mean is normal, and has effect on the difference between sample and population means")
    print()

#Define function for showing translated protien
#Function definition:
##protien_translation (sequence)
def translation(full_seq):
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
    ln  = 0

    
    for j in range(l-3): # iterate untill the end of the length 
        c3 = a[j]+a[j+1]+a[j+2]

        if c3 == 'ATG':   #for us to locate where all 'ATG' are
            start_pos = j
        
        #Start Iteratig while we found a ATG 
        for i in range(start_pos,l-3,3): #Start Iteratig while we found a ATG 
            c3 = a[i]+a[i+1]+a[i+2]

            #if c3 == 'TAA' or c3 == 'TAG' or c3 == 'TGA': stop here
            if c3 in stop_coden:
                #print(c3)
                break
            else:
                pro+= dna2aa[c3]
                seq+=c3
    return pro


#Define function for showing alignment score
#Function definition:
##Show_alignment_score (sequence 1, sequence 2, nof row, nof column, list)
def show_score(seq1, seq2, a,b,l):
    for i in range(a):
        if i < 1:  # i == 0
            print('     ',end="")
            for j in range(b-1):
                print("  ",seq1[j],end="")
            print()
        if i < 1:
            print(' ',end="")
        else:
            print(seq2[i-1],end="")
        for j in range(b):
            if l[i][j] < 0 or l[i][j] > 9:
                print(" ",l[i][j],end="")
            else:
                print("  ",l[i][j],end="")
        print()

#Define function for showing alignment direction
#Function definition:
##Show_alignment_direction ( sequence 1, sequence 2, nof row, nof column, list)
def show_direction(seq1, seq2, a,b,l):
    for i in range(a):
        if i < 1:  # i == 0
            print('     ',end="")
            for j in range(b-1):
                print("  ",seq1[j],end="")
            print()
        if i < 1:
            print(' ',end="")
        else:
            print(seq2[i-1],end="")
        for j in range(b):
                print("  ",l[i][j],end="")
        print()

#Define function for calculating local alignment score
#Function definition:
##calculating_local_alignment_score (sequence 1, sequence 2)
def local_alignment(seq1, seq2):

    l1 = len(seq1)+1
    l2 = len(seq2)+1
    score = [[0]*l1 for i in range(l2)]
    dir = [['X']*l1 for i in range(l2)]
    max_score = -1
    x = l1 - 1
    y = l2 - 1
    for i in range(l2):
        if i == 0:
            for j in range(1,l1):
                score[0][j] = 0 #score[0][j-1]+gap
        else:
            for j in range(l1):
                if j == 0:
                    score[i][0] = 0 #score[i-1][0]+gap
                else:
                    if seq1[j-1] == seq2[i-1]: # match
                        dig = score[i-1][j-1]+match
                    else:
                        dig = score[i-1][j-1]+mismatch
                    top = score[i-1][j] + gap
                    left = score[i][j-1]+gap
                    score[i][j] = max(dig,top,left,0)
                    if score[i][j] >= max_score:
                        max_score = score[i][j]
                        x = j
                        y = i
                    if score[i][j] == dig:
                        dir[i][j] = 'D'
                    elif score[i][j] == top:
                        dir[i][j] = 'T'
                    elif score[i][j] == left:
                        dir[i][j] = 'L'

    #Print score matrix
    print("Local Score Matrix: ")
    show_score(seq1, seq2, l2, l1, score)

    #Print direction matrix
    print("Show Direction Matrix: ")
    show_direction(seq1,seq2,l2,l1,dir)

    final_score = score[y][x]
    align1 = ''
    align2 = ''
    align3 = ''
    while(score[y][x]>0):
        if dir[y][x] == 'D':
            align1 = seq1[x-1]+align1
            align3 = seq2[y-1]+align3
            if seq1[x-1] == seq2[y-1]:
                align2 = '|'+align2
            else:
                align2 = ' '+align2
            x-=1
            y-=1
        elif dir[y][x] == 'T':
            align1 = '-'+align1
            align3 = seq2[y-1]+align3
            align2 = ' '+align2
            y-=1
        elif dir[y][x] == 'L':
            align1 = seq1[x-1]+align1
            align3 = '-'+align3
            align2 = ' '+align2
            x-=1

    #Print local alignment
    print("Final Local Alignment")
    print(align1)
    print(align2)
    print(align3)

    #print(max_score)
    return max_score


#Define function for calculating global alignment score
#Function definition:
##calculating_global_alignment_score (sequence 1, sequence 2)
def global_alignment(seq1, seq2):

    l1 = len(seq1)+1
    l2 = len(seq2)+1
    score = [[0]*l1 for i in range(l2)]
    dir = [['X']*l1 for i in range(l2)]
    for i in range(l2):
        if i == 0:
            for j in range(1,l1):
                score[0][j] = score[0][j-1]+gap
        else:
            for j in range(l1):
                if j == 0:
                    score[i][0] = score[i-1][0]+gap
                else:
                    if seq1[j-1] == seq2[i-1]: # match
                        dig = score[i-1][j-1]+match
                    else:
                        dig = score[i-1][j-1]+mismatch
                    top = score[i-1][j] + gap
                    left = score[i][j-1]+gap
                    score[i][j] = max(dig,top,left)
                    if score[i][j] == dig:
                        dir[i][j] = 'D'
                    elif score[i][j] == top:
                        dir[i][j] = 'T'
                    elif score[i][j] == left:
                        dir[i][j] = 'L'
    
    # #Print score matrix
    print("Global Score Matrix: ")
    show_score(seq1, seq2, l2, l1, score)

    #Print direction matrix
    print("Show Direction Matrix: ")
    show_direction(seq1,seq2,l2,l1,dir)

    x = l1-1
    y = l2-1
    final_score = score[y][x]
    align1 = ''
    align2 = ''
    align3 = ''
    while(x>0 and y>0):
        if dir[y][x] == 'D':
            align1 = seq1[x-1]+align1
            align3 = seq2[y-1]+align3
            if seq1[x-1] == seq2[y-1]:
                align2 = '|'+align2
            else:
                align2 = ' '+align2
            x-=1
            y-=1
        elif dir[y][x] == 'T':
            align1 = '-'+align1
            align3 = seq2[y-1]+align3
            align2 = ' '+align2
            y-=1
        elif dir[y][x] == 'L':
            align1 = seq1[x-1]+align1
            align3 = '-'+align3
            align2 = ' '+align2
            x-=1
    
    #Print global alignment
    print("Final Global Alignment")
    print(align1)
    print(align2)
    print(align3)
    #print(final_score)
    return final_score


def main():

    # Test sequence
    seq1 = 'GCATGCT'
    seq2 = 'GATTACA'

    # Real sequence
    # seq1 = 'gcctggactccacatccgtcctggactgttgagcgcgcagaccagaggcggttgaggaccagtggtgaggaacggccgaggcggcgtctgagcgggtctccggagttcagcatgcgtgagtgtatctctatccacgtggggcaggcaggtgtccagatcggcaatgcctgctgggaactgtactgccttgaacatggcattcagcctgacggtcagatgccaagcgacaaaaccattggcggcggggacgactcattcaacacattcttcagtgagactggagccggcaagcacgtgcccagggcagtgtttgtggacctggagcccactgtggtggatgaggtgcgcacgggaacctaccggcagctttttcacccagagcagctgatcactggaaaggaagatgcagccaataattatgccagaggccactacaccatcggcaaagagattgtcgacctggtcctggatcgaatccgaaagctggccgatctgtgcacgggactgcagggcttcctcatcttccacagctttggaggaggcacagggtctgggtttgcatcgctgctgatggagcggctttcagtggactatggcaagaagtccaagctggagtttgccatctacccagccccccaggtttctacagcggtcgtggagccttacaactccatcctgaccacgcacaccaccctagagcattccgactgtgctttcatggtggataacgaagccatctacgacatctgccggcgcaacctggatattgaacgtcccacatacaccaacctcaatcgtctgattgggcagattgtgtcgtccattacagcctccctgaggtttgatggcgccctgaatgtggacttaacagaattccagaccaacctggtgccataccctcgcatccacttcccactggcaacctacgccccggtcatctcagctgagaaggcataccatgagcagctgtcagtggcagagatcaccaatgcttgcttcgagccagccaatcagatggtcaagtgtgaccctcgccacggcaaatacatggcctgctgcatgttgtaccggggggatgtggttcccaaagatgtcaacgcggctattgcaaccatcaagacaaagcgcaccatccagtttgtagattggtgtccgactggatttaaggtgggtattaactaccagcctcccactgtggtccctgggggagacctggccaaagtgcagcgggccgtgtgcatgctgagcaataccacggccatcgcagaggcctgggcccgcctggaccacaaatttgacctcatgtacgccaagcgagcctttgtgcattggtacgtgggagaaggcatggaggaaggggagttctccgaggcccgggaggacctggcagcgctggagaaggactatgaagaggtgggcgtggattccgtggaagcagaggcagaagaaggggaggagtactgagcgcatgggtctggctggcggccgtccatttatgtcttccccaccattggaaataaaggatatattattaaaatttctagactgaggc'
    # seq1 = seq1.upper()

    # seq2 = 'agtgcgttacttacctcgactcttagcttgtcggggacggtaaccgggacccggtgtctgctcctgtcgccttcgcctcctaatccctagccactatgcgtgagtgcatctccatccacgttggccaggctggtgtccagattggcaatgcctgctgggagctctactgcctggaacacggcatccagcccgatggccagatgccaagtgacaagaccattgggggaggagatgactccttcaacaccttcttcagtgagacgggcgctggcaagcacgtgccccgggctgtgtttgtagacttggaacccacagtcattgatgaagttcgcactggcacctaccgccagctcttccaccctgagcagctcatcacaggcaaggaagatgctgccaataactatgcccgagggcactacaccattggcaaggagatcattgaccttgtgttggaccgaattcgcaagctggctgaccagtgcaccggtcttcagggcttcttggttttccacagctttggtgggggaactggttctgggttcacctccctgctcatggaacgtctctcagttgattatggcaagaagtccaagctggagttctccatttacccagcaccccaggtttccacagctgtagttgagccctacaactccatcctcaccacccacaccaccctggagcactctgattgtgccttcatggtagacaatgaggccatctatgacatctgtcgtagaaacctcgatatcgagcgcccaacctacactaaccttaaccgccttattagccagattgtgtcctccatcactgcttccctgagatttgatggagccctgaatgttgacctgacagaattccagaccaacctggtgccctacccccgcatccacttccctctggccacatatgcccctgtcatctctgctgagaaagcctaccatgaacagctttctgtagcagagatcaccaatgcttgctttgagccagccaaccagatggtgaaatgtgaccctcgccatggtaaatacatggcttgctgcctgttgtaccgtggtgacgtggttcccaaagatgtcaatgctgccattgccaccatcaaaaccaagcgcagcatccagtttgtggattggtgccccactggcttcaaggttggcatcaactaccagcctcccactgtggtgcctggtggagacctggccaaggtacagagagctgtgtgcatgctgagcaacaccacagccattgctgaggcctgggctcgcctggaccacaagtttgacctgatgtatgccaagcgtgcctttgttcactggtacgtgggtgaggggatggaggaaggcgagttttcagaggcccgtgaagatatggctgcccttgagaaggattatgaggaggttggtgtggattctgttgaaggagagggtgaggaagaaggagaggaatactaattatccattccttttggccctgcagcatgtcatgctcccagaatttcagcttcagcttaactgacagacgttaaagctttctggttagattgttttcacttggtgatcatgtcttttccatgtgtacctgtaatatttttccatcatatctcaaagtaaagtcattaacatcaaaa'
    # seq2 = seq2.upper()

    #Load sequence from file 
    #f1 = open('./a.txt', 'r')
    #seq1 = f1.readline()
    #f1.close()

    #f2 = open('./b.txt', 'r')
    #seq2 = f2.readline()
    #f2.close()

    print("Sequence 1:", seq1)
    print()
    print("Sequence 2:",seq2)
    print()

    #Task 1: Align 2 sequence with global alignment
    print("---------------------------Task 1: Global sequence calculation-------------------------------")
    gscore = global_alignment(seq1,seq2)
    print('Score of the global alignment:', gscore)
    print()

    #Task 2: Align 2 sequence with local alignment
    print("---------------------------Task 2: Local sequence calculation-------------------------------")
    lscore = local_alignment(seq1,seq2)
    print('Score of the local alignment:', lscore)
    print()

    #Task 3: Translate both sequence
    print("---------------------------Task 3: Translate both sequence-------------------------------")

    pro1 = translation(seq1)
    pro2 = translation(seq2)
    print("Translated sequence 1:",pro1)
    print("Translated sequence 2:",pro2)
    print()

    #Task 4: Align 2 Translated sequence
    print("---------------------------Task 4: Align 2 Translated sequence-------------------------------")
    
    #Align 2 Translated sequence with global alignment
    progscore = global_alignment(pro1,pro2)
    print('Score of the translated global alignment:', progscore)

    #Align 2 Translated sequence with local alignment
    prolscore = local_alignment(pro1,pro2)
    print('Score of the translated local alignment:', prolscore)
    print()

    #Task 5: Generated sequences at least 1,000 times and check global and local
    print("---------------------------Task 5: Generated sequences at least 1,000 times and check global and local-------------------------------")

    #generate DNA sequence of len 1000
    N = 1000
    ln1 = len(seq1)
    gen = gen_ran1000(ln1*3)
    #translation
    gen1 = translation(gen)

    ln2 = len(seq2)

    #generate 1000 DNA sequence
    for i in range (N):
        #generate DNA sequence of len seq2
        gen = gen_ran1000(ln2*3)
        #translation
        gen2 = translation(gen)

        print("Generated sequence 1: ", gen1)
        print("Generated sequence 2: ", gen2)
        print()

        #Align 2 Translated sequence with global alignment
        progscore = global_alignment(gen1,gen2)
        print('Score of the translated global alignment ',i+1,':', progscore)
        st_global.append(progscore)

        #Align 2 Translated sequence with local alignment
        prolscore = local_alignment(gen1,gen2)
        print('Score of the translated local alignment ',i+1,':', prolscore)
        st_local.append(prolscore)
        print()

    #Task 6: Average and standard deviation of alignment scores
    print("---------------------------Task 6: Average and standard deviation of alignment scores-------------------------------")

    avg_global = mean(st_global,N)
    avg_local = mean(st_local,N)

    std_global = stdev(st_global,N)
    std_local = stdev(st_local,N)

    print("Average global alignment score: ", avg_global)
    print("Average local alignment score: ", avg_local)

    print("Stddev global alignment score: ", std_global)
    print("Stddev local alignment score: ", std_local)
    print()

    #Task 7: Significance
    print("---------------------------Task 7: Significance-------------------------------")

    # Define labels, positions, bar heights and error bar heights
    #ref: https://problemsolvingwithpython.com/06-Plotting-with-Matplotlib/06.07-Error-Bars/
    #https://thedatascientist.com/how-to-do-a-t-test-in-python/

    #------------------------------------------------T-test-----------------------------------------------------------

    T_test(st_global, avg_global, std_global, 'Global')
    T_test(st_local, avg_local, std_local, 'Local')
    print()

    
    #-----------------------------------------------------Error Bar plot---------------------------------------------------------
    labels = ['Global', 'Local']
    x_pos = np.arange(len(labels))
    CTEs = [avg_global, avg_local]
    error = [std_global, std_local]

    # Build the plot
    fig, ax = plt.subplots()
    ax.bar(x_pos, CTEs,yerr=error,align='center',alpha=0.5,ecolor='black',capsize=10)
    ax.set_ylabel('Alignment (Error)')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels)
    ax.set_title('Error bar in alignment')
    ax.yaxis.grid(True)

    # Save the figure and show
    plt.tight_layout()
    plt.savefig('Bar_plot_with_error_bars.png')
    plt.show()



main()