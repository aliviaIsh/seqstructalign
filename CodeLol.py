
# Import statements
import numpy as np
import csv
from Alignment import *
from Sequence import *

out = 'pos_out_5.csv'
i_weight = 1.5
sim_weight = 2
size_weight = 3
hydro_weight = 1
gap_weight = 3
len_weight = -2

# put fasta files into csv format -------------------------------------------------------------------------------------
def conversion(infile):
    name = infile
    outname = infile + "_out.csv"
    f = open(name, 'r')

    # calculate the number of lines in the txt file
    Lines = f.readlines()
    count_lines = 0
    count_seq = 0

    # Strips the newline character
    for line in Lines:
        count_lines += 1
        i = line.count(">")
        count_seq += i

    # make an empty 2d list with proper identifiers
    rows, cols = (count_seq + 1, 9)
    arr = [[0 for i in range(cols)] for j in range(rows)]

    arr[0][0] = "Score Total"
    arr[0][1] = "NCBI Identifier"
    arr[0][2] = "Sequence"
    arr[0][3] = 'Identity Score'
    arr[0][4] = 'Similarity Score'
    arr[0][5] = 'Size Score'
    arr[0][6] = 'Hydropathy Score'
    arr[0][7] = 'Gap Penalty'
    arr[0][8] = "Length Penalty"


    h = 1
    num1 = 0

    temp = ""
    for line in Lines:
        if line[0] == ">":
            num1 = line.find(" ")
            s = line[num1 + 1:]
            nline = line[1:num1]
            num2 = line.find("[")
            aline = line[num1 + 1: num2 -1]
            num3 = line.find("]")
            oline = line[(num2 + 1): (num3)]

            if nline.find("GLN") == -1:

                arr[h][1] = nline
                arr[h][2] = oline

                h += 1
                temp = ""
        else:
            temp = temp + line
            temp = temp.strip("\n")
            arr[h-1][2] = temp

            arr[h-1][8] = findLen(arr[h-1][2])

    x = len(arr)
    i = 0
    while i < x:
        string = str(arr[i][2])
        s2 = string.replace("-", "O")
        arr[i][2] = s2
        i = i + 1

    a = np.array(arr)

    return a
# Termination of the conversion of .txt to csv ----------------------------------------------------------------------

# get hydropathy of residues
def get_hpathy(res):
    if res == 'I':
        return 4.5
    elif res == 'V':
        return 4.2
    elif res == 'L':
        return 3.8
    elif res == 'F':
        return 2.8
    elif res == 'C':
        return 2.5
    elif res == 'M':
        return 1.9
    elif res == 'A':
        return 1.8
    elif res == 'W':
        return -0.9
    elif res == 'G':
        return -0.4
    elif res == 'T':
        return -0.7
    elif res == 'S':
        return -0.8
    elif res == 'Y':
        return -1.3
    elif res == 'P':
        return -1.6
    elif res == 'H':
        return -3.2
    elif res == 'N':
        return -3.5
    elif res == 'D':
        return -3.5
    elif res == 'Q':
        return -3.5
    elif res == 'E':
        return -3.5
    elif res == 'K':
        return -3.9
    elif res == 'R':
        return -4.5
    else:
        return 0

# Scoring methods/matrixes ------------------------------------------------------------------------------------------
# High number = very similar   Low number = very dissimilar
# Score based on identity MAKE CERTAIN RESIDUE IDENTITIES WORTH MORE
def score_i(align_1, align_2):
    rows, cols = (2, len(align_1))
    r = [[0 for i in range(cols)] for j in range(rows)]
    p= 0
    for i in align_1:
        r[0][p] = str(i)
        p += 1

    p = 0
    for i in align_2:
        r[1][p] = str(i)
        p += 1

    length = len(r[1])
    i = 0
    score = 0
    while i < length:
        if r[0][i] == r[1][i]:
            score += 1
        i += 1

    return score

# Score based on chemistry
def score_c(align_1, align_2):
    rows, cols = (2, len(align_1))
    r = [[0 for i in range(cols)] for j in range(rows)]
    p = 0
    for i in align_1:
        r[0][p] = str(i)
        p += 1

    p = 0
    for i in align_2:
        r[1][p] = str(i)
        p += 1

    length = len(r[1])
    i = 0
    score = 0
    while i < length:
        # if both residues are positively charged
        if str(r[0][i]) in "RHK" and str(r[1][i]) in "RHK":
            score += 3
        # if both residues are negatively charge
        elif str(r[0][i]) in "DE" and str(r[1][i]) in "DE":
            score += 3
        # if both residues are OH
        elif str(r[0][i]) in "ST" and str(r[1][i]) in "ST":
            score += 2
        # if both residues are amines
        elif str(r[0][i]) in "NQ" and str(r[1][i]) in "NQ":
            score += 2
        # if both residues are sulfur friends
        elif str(r[0][i]) in "CM" and str(r[1][i]) in "CM":
            score += 2
        # if both residues are smallish and hydrophobic
        elif str(r[0][i]) in "GAVIL" and str(r[1][i]) in "GAVIL":
            score += 1
        # if both residues are largeish and hydrophobic
        elif str(r[0][i]) in "FYW" and str(r[1][i]) in "FYW":
            score += 1

        i += 1

    return score

# Score based on size
def score_s(align_1, align_2):
    rows, cols = (2, len(align_1))
    r = [[0 for i in range(cols)] for j in range(rows)]
    p = 0
    for i in align_1:
        r[0][p] = str(i)
        p += 1

    p = 0
    for i in align_2:
        r[1][p] = str(i)
        p += 1

    length = len(r[1])
    i = 0
    score = 0
    while i < length:
        # if both residues are 189-228 A
        if str(r[0][i]) in "FWY" and str(r[1][i]) in "FWY":
            score += 3
        # if both residues are 162-174 A
        elif str(r[0][i]) in "ILMKR" and str(r[1][i]) in "ILMKR":
            score += 3
        # if both residues are 138-154 A
        elif str(r[0][i]) in "VHEQ" and str(r[1][i]) in "VHEQ":
            score += 2
        # if both residues are 108-117 A
        elif str(r[0][i]) in "CPTDN" and str(r[1][i]) in "CPTDN":
            score += 2
        # if both residues are 60-90 A
        elif str(r[0][i]) in "AGS" and str(r[1][i]) in "AGS":
            score += 2
        # left overs
        else:
            score -= 0.5

        i += 1

    return score

# Score based on hydropathy using actual hydropathy indexes
def score_h(align_1, align_2):
    rows, cols = (2, len(align_1))
    r = [[0 for i in range(cols)] for j in range(rows)]
    p = 0
    for i in align_1:
        r[0][p] = str(i)
        p += 1

    p = 0
    for i in align_2:
        r[1][p] = str(i)
        p += 1

    length = len(r[1])
    i = 0
    score = 0
    while i < length:
        score += (9 - abs((float(get_hpathy(r[0][i]))) - (float(get_hpathy(r[1][i])))))/20
        i += 1
    return score

# Score based on gaps ADD METHOD TO MAKE LINES OF ----- CORRESPOND TO MORE POINT LOSS
def score_g(align_1, align_2):
    val_1 = align_1.count("-")
    val_2 = align_2.count("-")

    return -1*(val_1 + val_2)

def findLen(str):
    counter = 0
    for i in str:
        counter += 1
    return counter

# Scoring total + weights and polishing off the file
def totalscores(a):
# iterate through each sequence
    length = len(a)
    S = str(a[1][2])
    i = 2
    while i < length:
        print(str(i/length * 100) + "% Complete")
        T = str(a[i][2])
        # construct alignment
        value = allignSequences(S, T)
        string = str(value[1])
        parse = string.find(":")
        s_align = string[:parse]
        t_align = string[parse + 1:]

        # send them off to the scoring methods
        a[i][3] = float(score_i(s_align, t_align))
        a[i][4] = float(score_c(s_align, t_align))
        a[i][5] = float(score_s(s_align, t_align))
        a[i][6] = float(score_h(s_align, t_align))
        a[i][7] = float(score_g(s_align, t_align))

        l = findLen(T)
        if l > 600:
            a[i][8] = 200
        else:
            a[i][8] = 0

        a[i][0] = (i_weight*float(a[i][3])) + (sim_weight*float(a[i][4])) + (size_weight*float(a[i][5])) + (hydro_weight*float(a[i][6])) + (gap_weight*float(a[i][7]) + (len_weight*float(a[i][8])))
        i += 1

    with open(out, 'w', newline='') as file:
        mywriter = csv.writer(file, delimiter=',')
        mywriter.writerows(a)
# Close files