
import os
import glob
import sys


itr = 0
nucount = 0
templist = []

global user_choice 
#global gapscore
global protien_list1 
global protien_list2 

f = " "
matrix_alphabet = []

path = "./P1SubMatrices"

selected_MATRIX = []



def cls():
    os.system('cls' if os.name=='nt' else 'clear')

def set_sequenceA():
    filesource1 = open("./P1AASeqs/sequenceA1.txt", "r")

def set_sequenceB():
    filesource1 = open("./P1AASeqs/sequenceA1.txt", "r")
def set_sequenceC():
    filesource1 = open("./P1AASeqs/sequenceA1.txt", "r")

def setas_NUCLEOTIDE():
    filesource = open("./P1SubMatrices/AAnucleoPP.txt", "r")
    return populate_MATRIX(filesource)

def setas_BLOSUM():
    filesource = open("./P1SubMatrices/BLOSUM62.txt", "r")
    return populate_MATRIX(filesource)

def setas_HP():
    filesource = open("./P1SubMatrices/HP.txt", "r")
    return populate_MATRIX(filesource)

def setas_PAM():
    filesource = open("./P1SubMatrices/PAM250-scores.txt", "r")
    return populate_MATRIX(filesource)


def find_MAX(x,y,sub_value):

    

    d = int(sub_value) + int(global_MATRIX[x-1][y-1])
    h = int(-gapscore) + int(global_MATRIX[x][y-1])
    v = int(-gapscore ) + int(global_MATRIX[x-1][y])

    maxvalue = 0

    if(d >= h and d >= v):
        maxvalue = d
        global_MAP[x][y] = 'd'
        semilocal_MAP[x][y] = 'd'
        


    elif(h >= d and h >= v):
        maxvalue = h
        global_MAP[x][y] = 'h'
        semilocal_MAP[x][y] = 'h'
        

   
    else:
        
        maxvalue = v
        global_MAP[x][y] = 'v'
        semilocal_MAP[x][y] = 'v'
        

    #print(d)
    #print(h)
    #print(v)


    return maxvalue

def find_localMAX(x,y,sub_value):

    

    d = int(sub_value) + int(global_MATRIX[x-1][y-1])
    h = int(-gapscore) + int(global_MATRIX[x][y-1])
    v = int(-gapscore ) + int(global_MATRIX[x-1][y])

    maxvalue = 0

    if(d >= h and d >= v and d >= 0):
        maxvalue = d
        global_MAP[x][y] = 'd'
        semilocal_MAP[x][y] = 'd'
        local_MAP[x][y] = 'd'


    elif(h >= d and h >= v and h >= 0):
        maxvalue = h
        global_MAP[x][y] = 'h'
        semilocal_MAP[x][y] = 'h'
        local_MAP[x][y] = 'h'

    elif(0 >= d and 0 >= h and 0 >= v):
        maxvalue = 0


    else:
        
        maxvalue = v
        global_MAP[x][y] = 'v'
        semilocal_MAP[x][y] = 'v'
        local_MAP[x][y] = 'v'

    #print(d)
    #print(h)
    #print(v)


    return maxvalue

def find_localMAX_float(x,y,sub_value):

    

    d = float(sub_value) + float(global_MATRIX[x-1][y-1])
    h = float(-gapscore) + float(global_MATRIX[x][y-1])
    v = float(-gapscore ) + float(global_MATRIX[x-1][y])

    maxvalue = 0

    if(d >= h and d >= v and d >= 0):
        maxvalue = d
        global_MAP[x][y] = 'd'
        semilocal_MAP[x][y] = 'd'
        local_MAP[x][y] = 'd'


    elif(h >= d and h >= v and h >= 0):
        maxvalue = h
        global_MAP[x][y] = 'h'
        semilocal_MAP[x][y] = 'h'
        local_MAP[x][y] = 'h'

    elif(0 >= d and 0 >= h and 0 >= v):
        maxvalue = 0


    else:
        
        maxvalue = v
        global_MAP[x][y] = 'v'
        semilocal_MAP[x][y] = 'v'
        local_MAP[x][y] = 'v'

    #print(d)
    #print(h)
    #print(v)


    return float(maxvalue)



def find_MAX_float(x,y,sub_value):

    

    d = float(sub_value) + float(global_MATRIX[x-1][y-1])
    h = float(-gapscore) + float(global_MATRIX[x][y-1])
    v = float(-gapscore ) + float(global_MATRIX[x-1][y])

    maxvalue = 0

    if(d >= h and d >= v):
        maxvalue = d
        global_MAP[x][y] = 'd'
        semilocal_MAP[x][y] = 'd'
        local_MAP[x][y] = 'd'


    elif(h >= d and h >= v):
        maxvalue = h
        global_MAP[x][y] = 'h'
        semilocal_MAP[x][y] = 'h'
        local_MAP[x][y] = 'h'


    else:
        maxvalue = v
        global_MAP[x][y] = 'v'
        semilocal_MAP[x][y] = 'v'
        local_MAP[x][y] = 'v'
        

    #print(d)
    #print(h)
    #print(v)


    return float(maxvalue)

    
def calculate_cell(x,y,sub_matrix):

    pos1 = 0
    pos2 = 0
    #print(x,y)
    if(int(x) < len(protien_list1)):
        #print(seq1[x], sub_matrix[0][x])     

        while(str(protien_list1[x]) != str(sub_matrix[0][pos1])):        
            
            pos1 = pos1 + 1

    if(int(y) < len(protien_list2)):
        #print(seq2[y], sub_matrix[0][y])
        while(protien_list2[y] != sub_matrix[0][pos2]):        
            
            pos2 = pos2 + 1

    #print(sub_matrix[pos1+1][pos2], glob_matrix[x-1][y-1], glob_matrix[x-1][y], glob_matrix[x][y-1])
    
    if(matrixType == "hp"):
        return find_MAX_float(x,y,sub_matrix[pos1+1][pos2])    
    else:
        return find_MAX(x,y,sub_matrix[pos1+1][pos2])


def calculate_localcell(x,y,sub_matrix):

    pos1 = 0
    pos2 = 0
    #print(x,y)
    if(int(x) < len(protien_list1)):
        #print(seq1[x], sub_matrix[0][x])     

        while(str(protien_list1[x]) != str(sub_matrix[0][pos1])):        
            
            pos1 = pos1 + 1

    if(int(y) < len(protien_list2)):
        #print(seq2[y], sub_matrix[0][y])
        while(protien_list2[y] != sub_matrix[0][pos2]):        
            
            pos2 = pos2 + 1

    #print(sub_matrix[pos1+1][pos2], glob_matrix[x-1][y-1], glob_matrix[x-1][y], glob_matrix[x][y-1])
    
    if(matrixType == "hp"):
        return find_localMAX_float(x,y,sub_matrix[pos1+1][pos2])    
    else:
        return find_localMAX(x,y,sub_matrix[pos1+1][pos2])


def find_globalalignment():
    maxrow = len(protien_list1) - 1
    maxcol = len(protien_list2) - 1
    max_x = maxrow
    max_y = maxcol
    itr = 0
    MAXcell = -999

    alignment1 = ""
    alignment2 = ""

    nxtMAP = global_MAP[max_x][max_y]
    nxtCELL = global_MATRIX[max_x][max_y]

    #print(nxtMAP, nxtCELL)

    for i in range(maxrow):
        for j in range(maxcol):
            if(nxtMAP == 'd'):
                #print("Current MAP: ", global_MAP[max_x][max_y], "List 1 letter:", protien_list1[max_x], "List 2 letter:", protien_list2[max_y], global_MATRIX[max_x][max_y])
                
                alignment1 = alignment1 + protien_list1[max_x]
                alignment2 = alignment2 + protien_list2[max_y]
                
                max_x = max_x - 1
                max_y = max_y - 1
                nxtMAP = global_MAP[max_x][max_y]
                nxtCELL = global_MATRIX[max_x][max_y]

                #print("Next MAP: ", global_MAP[max_x][max_y], "List 1 letter:", protien_list1[max_x], "List 2 letter:", protien_list2[max_y], global_MATRIX[max_x][max_y], "\n")

                
                i = i + 1

            elif(nxtMAP == 'h'):
                #print("Current MAP: ", global_MAP[max_x][max_y], "List 1 letter:", protien_list1[max_x], "List 2 letter:", protien_list2[max_y], global_MATRIX[max_x][max_y])

                alignment1 = alignment1 + protien_list1[max_x]
                alignment2 = alignment2 + "_"
                
                
                max_x = max_x - 1
                nxtMAP = global_MAP[max_x][max_y]
                nxtCELL = global_MATRIX[max_x][max_y]

                #print("Next MAP: ", global_MAP[max_x][max_y], "List 1 letter:", protien_list1[max_x], "List 2 letter:", protien_list2[max_y], global_MATRIX[max_x][max_y], "\n")
                j = j + 1
                
                #j = j - 1

            elif(nxtMAP == 'v'):
                #print("Current MAP: ", global_MAP[max_x][max_y], "List 1 letter:", protien_list1[max_x], "List 2 letter:", protien_list2[max_y], global_MATRIX[max_x][max_y])
                
                alignment1 = alignment1 + "_"
                alignment2 = alignment2 + protien_list2[max_y]

                max_y = max_y - 1
                nxtMAP = global_MAP[max_x][max_y]
                nxtCELL = global_MATRIX[max_x][max_y]
                i = i + 1
                #print("Next MAP: ", global_MAP[max_x][max_y], "List 1 letter:", protien_list1[max_x], "List 2 letter:", protien_list2[max_y], global_MATRIX[max_x][max_y], "\n")                

            else:

                if(global_MAP[max_x][max_y] != None):

                    #print("Current MAP: ", global_MAP[max_x][max_y], "List 1 letter:", protien_list1[max_x], "List 2 letter:", protien_list2[max_y], global_MATRIX[max_x][max_y])
                    
                    alignment1 = alignment1 + "_"
                    alignment2 = alignment2 + "_"

                    max_y = max_y - 1
                    nxtMAP = global_MAP[max_x][max_y]
                    nxtCELL = global_MATRIX[max_x][max_y]
                    i = i + 1
                    #print("Next MAP: ", global_MAP[max_x][max_y], "List 1 letter:", protien_list1[max_x], "List 2 letter:", protien_list2[max_y], global_MATRIX[max_x][max_y], "\n")                
                


    align1 = ""
    align2 = ""

    print("Global Alignment:")

    print(alignment1, "\n", alignment2)
   
    #print(global_MAP[max_x][max_y], global_MATRIX[max_x][max_y])
    

    

def find_localalignment():
    print()

def find_semilocalalignment():
    maxrow = len(protien_list1) - 1
    maxcol = len(protien_list2) - 1
    max_x = 0
    max_y = 0
    itr_row = 0
    itr_col = 0

    alig1 = ""
    alig2 = ""
    
    MAXcell = -999

    while(itr_row <= maxrow):
        while(itr_col <= maxcol):
            cur = global_MATRIX[itr_row][maxcol]

            if(global_MATRIX[maxrow][itr_col] >= MAXcell):
                MAXcell = global_MATRIX[maxrow][itr_col]
                max_x = maxrow
                max_y = itr_col

            if(global_MATRIX[itr_row][maxcol] >= MAXcell):
                MAXcell = global_MATRIX[itr_row][maxcol]
                max_x = itr_row
                max_y = maxcol

            itr_col = itr_col + 1

        itr_row = itr_row + 1


    for i in range(maxrow):
        for j in range(maxcol):

            nxt_MAP = global_MAP[max_x][max_y]

            
            

            if(nxt_MAP == 'd'):
                alig1 = alig1 + protien_list1[max_x]
                alig2 = alig2 + protien_list2[max_y]

                max_x = max_x - 1
                max_y = max_y - 1
                

                i = i + 1

            elif(nxt_MAP == 'h'):

                alig1 = alig1 + "_"
                

                max_x = max_x - 1
                j = j + 1
                
            
            elif(nxt_MAP == 'v'):
                alig2 = alig2 + "_"

                max_y = max_y - 1
                i = i + 1

   

    print("Semi-Local Alignment:")
    print(alig1, "\n", alig2)

    

    

def populate_MATRIX(f): 

    i = 0

    title = f.readline()
    
    alphabet = f.readline()

    for letter in alphabet:
        if(letter != "," and  letter != "\n"):
            matrix_alphabet.append(letter)
            #print(letter)
    
    active_matrix = []
    active_matrix.clear()
    
    active_matrix.append(alphabet.strip("\n,").split(","))

    line = f.readline()

    while (line != ""): 

        active_matrix.append(line.strip("\n").split(","))
        line = f.readline()

    f.close()

    matrix_length = 0
    matrix_length = len(selected_MATRIX)

    itr = 0
    tempitr = 0
    itrscore = 0

    for x in range(matrix_length):
        for y in range(matrix_length):
            if(selected_MATRIX[itr] == "\n"):
                itr = itr+1
            #print(selected_MATRIX[itr], itrscore)
        
            itr = itr+1
            tempitr = tempitr+1
            itrscore = itrscore+1

    matrix_filename = str(f.name)
    #print(matrix_filename)

    user_choice = ""
    templist = []
    if(matrix_filename == "./P1SubMatrices/HP.txt"):
        for x in range(matrix_length):
            for y in range(matrix_length):
                #for x in range(matrix_length):
                if(x == 0):
                    templist.append(float(selected_MATRIX[y]))
                    #print("Make it here")
                else:
                    templist.append(float(selected_MATRIX[y+(20*x)]))
            #print("Temp: ", templist)
            
            active_matrix.append(templist)

            #Clears Temporary List
            templist = []

        for x in selected_MATRIX:
            if(x == "\n"):
                x.split("/n")
    else:
        for x in range(matrix_length):
            for y in range(matrix_length):
                #for x in range(matrix_length):
                if(x == 0):
                    templist.append(int(selected_MATRIX[y]))
                    #print("Make it here")
                else:
                    templist.append(int(selected_MATRIX[y+(20*x)]))        
            #print("Temp: ", templist)
            
            active_matrix.append(templist)

            #Clears Temporary List
            templist = []

        for x in selected_MATRIX:
            if(x == "\n"):
                x.split("/n")
        
    print(active_matrix)
    print(active_matrix[1][1]+ active_matrix[1][4])
    #wait = input()
    return active_matrix

def gapscore_choice():
    cls()
    print("Set a gap score by typing a number:")
    gap_scoreInput = abs(int(input()))

        
    return gap_scoreInput




def seq_setup(user_choice, read_MATRIX):

    global gapscore
    gapscore = gapscore_choice()
    
    cls()

    print("\n\nGap Score:", -gapscore )

    global protien_list1
    protien_list1 = []

    global protien_list2
    protien_list2 = []

    

    if(user_choice == "1"):

        if(matrixType == "neculeotide"):
            
            print("Enter the 1st sequence of the file: ")
            seq1 = input()

            print("Enter the 2nd sequence of the file: ")
            seq2 = input()

            
            seqA = open("./P1AASeqs/" + seq1, "r")
            seqB = open("./P1AASeqs/" + seq2, "r")
            
        else:
            seqA = open("./P1AASeqs/sequenceA1.txt", "r")
            seqB = open("./P1AASeqs/sequenceA2.txt", "r")

        
    elif(user_choice == "2"):

        if(matrixType == "neculeotide"):
            print("Enter the 1st sequence of the file: ")
            seq1 = input()

            print("Enter the 2nd sequence of the file: ")
            seq2 = input()

            
            seqA = open("./P1AASeqs/" + seq1, "r")
            seqB = open("./P1AASeqs/" + seq2, "r")
            
        else:
            seqA = open("./P1AASeqs/sequenceB1.txt", "r")
            seqB = open("./P1AASeqs/sequenceB2.txt", "r")

    elif(user_choice == "3"):

        if(matrixType == "neculeotide"):
            seqA = open("./P1AASeqs/Neuc1.txt", "r")
            seqB = open("./P1AASeqs/Neuc2.txt", "r")
            
        else:
            seqA = open("./P1AASeqs/sequenceC1.txt", "r")
            seqB = open("./P1AASeqs/sequenceC2.txt", "r")
        
        

    protine_title1 = seqA.readline()

    protine_title2 = seqB.readline()

    
    protien1 = seqA.read(1)
    protien2 = seqB.read(1)

    protien_list1.insert(0, " ")
    protien_list2.insert(0, " ")

    while(protien1 != "" and protien1 != "\n"):
    
        protien_list1.append(protien1)
        protien1 = seqA.read(1)
    
    while(protien2 != "" and protien2 != "\n"):
    
        protien_list2.append(protien2)
        protien2 = seqB.read(1)

    
    

    print("\n", protine_title1, protien_list1)
    print("\n", protine_title2, protien_list2)

    

    global global_MATRIX
    global semilocal_MATRIX
    global local_MATRIX

    global global_MAP
    global semilocal_MAP
    global local_MAP

    
    global_MAP = [[None]*len(protien_list2) for i in range(len(protien_list1))]
    semilocal_MAP = [[None]*len(protien_list2) for i in range(len(protien_list1))]
    local_MAP = [[None]*len(protien_list2) for i in range(len(protien_list1))]
    
    global_MATRIX = [[None]*len(protien_list2) for i in range(len(protien_list1))]
    semilocal_MATRIX = [[None]*len(protien_list2) for i in range(len(protien_list1))]
    local_MATRIX = [[None]*len(protien_list2) for i in range(len(protien_list1))]
   
    matrix_length = len(global_MATRIX[0]) 

    itrcount = 0
    xcount = 0
    ycount = 0

    
    if(matrixType == "hp"):

        for x in range(len(protien_list1)):
            for y in range(len(protien_list2)):
                
                if(x == 0):
                    global_MATRIX[x][y] = float(-gapscore * y)
                    semilocal_MATRIX[x][y] = 0
                    local_MATRIX[x][y] = 0 
                    

                elif(y == 0):
                    global_MATRIX[x][y] = float(-gapscore * x)
                    semilocal_MATRIX[x][y] = 0
                    local_MATRIX[x][y] = 0 

                else:
                    #print(protien_list1[x])
                    #print(protien_list2[y])
                    global_MATRIX[x][y] = float(calculate_cell(x,y,read_MATRIX))
                    semilocal_MATRIX[x][y] = float(calculate_cell(x,y,read_MATRIX))
                    local_MATRIX[x][y] = float(calculate_localcell(x,y,read_MATRIX))

    else:
        for x in range(len(protien_list1)):
            for y in range(len(protien_list2)):
                
                if(x == 0):
                    global_MATRIX[x][y] = -gapscore * y
                    semilocal_MATRIX[x][y] = 0
                    local_MATRIX[x][y] = 0 
                    

                elif(y == 0):
                    global_MATRIX[x][y] = -gapscore * x
                    semilocal_MATRIX[x][y] = 0
                    local_MATRIX[x][y] = 0 

                else:
                    #print(protien_list1[x])
                    #print(protien_list2[y])
                    global_MATRIX[x][y] = calculate_cell(x,y,read_MATRIX)
                    semilocal_MATRIX[x][y] = calculate_cell(x,y,read_MATRIX)
                    local_MATRIX[x][y] = calculate_localcell(x,y,read_MATRIX) 

        
    print("\n\t  ", "Global Matrix: ")

    for i in range(len(global_MATRIX)):
        print(global_MATRIX[i])

    print("\n\t  ", "Global MAP: ")

    for i in range(len(global_MAP)):
        print(global_MAP[i])

    print("\n\t  ", "Semi Local Matrix: ")

    for i in range(len(semilocal_MATRIX)):
        print(semilocal_MATRIX[i])

    print("\n\t  ", "Semi-local MAP: ")

    for i in range(len(semilocal_MATRIX)):
        print(semilocal_MAP[i])

    print("\n\t  ", "Local Matrix: ")
    
    for i in range(len(local_MATRIX)):
        print(local_MATRIX[i])

    print("\n\t  ", "local MAP: ")

    for i in range(len(local_MATRIX)):
        print(local_MAP[i])

    #print(global_MATRIX)

    xint = 0
    yint = 0

def type_selection(read_MATRIX):

    cls()
    gap_score = 0

    if(matrixType == "neculeotide"):
        print("Chooose a Neculeotide Sequence")
        print(" 1: Sequence 1 \n 2: Sequence 2 \n 3: Sequence 3 ")
        user_choice = input()

    else:

        user_choice = ""
        print("\nChoose a SEQUENCE\n")
        print(" 1: Sequence A \n 2: Sequence B \n 3: Sequence C \n")
        user_choice = input()

    seq_setup(user_choice,read_MATRIX)

    find_semilocalalignment()

    find_globalalignment()

    wait = input()


def menuPrompt():
    

    cls()

    global matrixType
    matrixType = ""
    

    user_choice = ""
    print("\n\t Choice a MATRIX\n")
    print(" 1: Neculeotide \n 2: BLOSUM \n 3: HP \n 4: PAM \n")
    user_choice = input()

    if user_choice == "1":
        matrixType = "neculeotide"
        print(user_choice)
        type_selection(setas_NUCLEOTIDE())

    elif user_choice == "2":
        matrixType = "blosum"
        print(user_choice)
        type_selection(setas_BLOSUM())        

    elif user_choice == "3":
        matrixType = "hp"
        print(user_choice)
        type_selection(setas_HP())

    elif user_choice == "4":
        matrixType = "pam"
        print(user_choice)
        type_selection(setas_PAM())
    
    
while(1):
    menuPrompt()

print(path)
