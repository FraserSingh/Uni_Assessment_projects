#!/usr/bin/python3
#Author B217754
#code which partially implements Smith-Waterman local alignment with Affine Gaps from first principles
#This code is adapted from and uses variable names from SmithWaterman_go.py where possible from BA Lecture 4

import argparse

parser = argparse.ArgumentParser(description='Aligning sequences...')
parser.add_argument('seq1',action="store",help="First sequence")
parser.add_argument('seq2',action="store",help="Second sequence")
parser.add_argument('extend_gap',action="store",help="Penalty for gap extension",type =int) #NOTE take gap extension penalties from the user
parser.add_argument('open_gap',action="store",help="Penalty for gap creation",type =int) #NOTE take gap creation penalties from the user
parser.add_argument('seqmatch',action="store",help="Score for matching bases, +1 recommended for DNA",type =int,default=+1) #NOTE take base match score from user
parser.add_argument('seqmismatch',action="store",help="Score for non-matching bases, -1 recommended for DNA",type =int,default=-1) #NOTE take base mimatch score from user

margs = parser.parse_args()
print(margs)

#Scoring system 
seqmatch =margs.seqmatch
seqmismatch=margs.seqmismatch
h=margs.open_gap 
g=margs.extend_gap 

#NOTE Ensures that extension penalty is smaller than opening penalty
if h>g:
    print("Penalty for opening a gap must be more positive than penalty for extending a gap")
    exit()
else:
    #matrix defintion needed here, traditionally, I would have used np.matrices, but we have to make the code compatible with BA-supplied code
    #make the matrices' overall structure NOTE this is a major change from standard SW
    def initialise_matrices():
        global rows, cols
        mats=[0,0,0] #first position refers to matrix M, second Ix, third Iy
        matrices_ = [[mats for col in range(cols+1)] for row in range(rows+1)]
        # set Ix and Iy first column as negative infinity 
        for row in range(0,rows+1):
            matrices_[row][0]=[0,float('-inf'),float('-inf')]
        # set Ix and Iy first row as negative infinity 
        for col in range(0,cols+1):
            matrices_[0][col]=[0,float('-inf'),float('-inf')]
        return matrices_

    def calc_score(matrices_,row, col): #perform recursion calculation for the current cell
        global seq1,seq2,seqmatch,seqmismatch
        #NOTE show which letters are being compared and their position in the sequence
        print("seq1:",seq1[col- 1]," seq2: "+seq2[row - 1],"\trow:",row," column:",col)
        #score takes 1 if there is no match, -1 if there is a match, previously called 'sc'
        match_status= seqmatch if seq1[col-1]==seq2[row-1] else seqmismatch 

        MatchXY=(matrices_[row-1][col-1][0]+match_status) #(cell one up/left in M)
        InsertX=(matrices_[row-1][col-1][1]+match_status) #(cell one up/left in Ix)
        InsertY=(matrices_[row-1][col-1][2]+match_status) #(cell one up/left in Iy)
        M_val=max(0,MatchXY,InsertX,InsertY)

        OpenGapX=(matrices_[row-1][col][0]+g+h) #(cell one left in M)+gap_ext+gap_opn 
        ExtendGapX=(matrices_[row-1][col][1]+g) #(cell one left in Ix)+gap_ext
        Ix_val=max(0,OpenGapX,ExtendGapX)

        OpenGapY=(matrices_[row][col-1][0]+g+h) #(cell one up in M)+gap_ext+gap_opn
        ExtendGapY=(matrices_[row][col-1][2]+g) #(cell one up in Iy)+gap_ext
        Iy_val=max(0,OpenGapY,ExtendGapY)
        #update the matrices' 
        vals=[M_val,Ix_val,Iy_val]
        print(f"vals = {vals}")
        return vals 

    #Generate scores for each cell by calling calc_score (builds the initial scoring matrix used for traceback)
    def build_matrix(matrices_):

        global rows, cols
        #matrices=[matrices_]
        #for every matrix, calculate the scores
        #for matrix in matrices:
            #for every cell, calculate the score and store it in the empty mymatrix
        for row in range(1, rows+1):
            for col in range(1, cols+1):
                    matrices_[row][col] = calc_score(matrices_, row, col)
        return matrices_

    #gets the max value from the built matrix
    def get_max(matrices_):
        #set up blank variables which update
        max=matrices_[0][0][0] #added matrix index
        mrow=0
        mcol=0
        mmatrix=0

        #define number of rows and columns
        rows = len(matrices_)
        cols = len(matrices_[0])
        matrixes=len(matrices_[0][0])

        #Check every cell in scoring matrix, if the cell's value (produced by function ____) is higher than the top left cell, update the highest score and record the cell coordinates
        #From the bottom right, so that identical alignments are found
        for row in reversed(range(1, rows)):
            for col in reversed(range(1, cols)):
                for matrix in range(0,matrixes):
                    if matrices_[row][col][matrix]>max:
                        max=matrices_[row][col][matrix]
                        mrow=row
                        mcol=col
                        mmatrix=matrix

        #return the M Ix and Iy scores of the cell with the highest score
        return [mrow,mcol,mmatrix] #used to return maxv

    #NOTE First instance of maxv is the coordinates for the highest scoring number across all matrices, provided by get_max()
    #Takes the matrix, the cell of importance (with the current highest score) and returns the coordinates of the origin cell
    def traceback_affine(matrices_,maxv,state): 
        global h,g #import gap penalties
        state = state
        mrow=maxv[0]
        mcol=maxv[1]
        val=max(matrices_[mrow][mcol]) #the score in question (highest in the cell) 
        match_status= seqmatch if seq1[mcol-1] == seq2[mrow-1] else seqmismatch #NOTE refers to position one up/left diagonally of the highest scoring match (by referencing seqeunces highest scoring letters)
        print(f"Current cell's scores: {matrices_[mrow][mcol]} at coordinates: {mcol,mrow}")

    #deal with being in top row ()
        if mrow==0:
            state='Iy'
            val=matrices_[mrow][mcol-1][2]
            maxv=[mrow,mcol-1,2]
    #deal with being in first column
        elif mcol==0:
            state='Ix'
            val=matrices_[mrow-1][mcol][1]
            maxv=[mrow-1,mcol,1]
        elif val==0:
            return
    #Determine origin of M coordinate, 3 potential origins for M, each from cell on up left diagonal
        if state == 'M':
            if val == matrices_[mrow-1][mcol-1][0] + match_status:
                maxv=[mrow-1,mcol-1,0] 
                state='M'
                val=matrices_[mrow-1][mcol-1][0]# MatchXY, (cell one up/left in M)
                print('MatchXY')
            elif val == matrices_[mrow-1][mcol-1][1] + match_status: 
                #[mrow-1][mcol-1][1] index refers to the Ix score, so the state variable is updated to reflect that we are moving into that matrix
                state='Ix'
                maxv=[mrow-1,mcol-1,1] #the Ix coordinate
                val=matrices_[mrow-1][mcol-1][1]# InsertX, (cell one up/left in Ix)
                print('InsertX')                
            elif val == matrices_[mrow-1][mcol-1][2] + match_status: 
                #[mrow-1][mcol-1][2] index refers to the Iy score, so the state variable is updated to reflect that we are moving into that matrix
                state='Iy'
                maxv=[mrow-1,mcol-1,2] #the Iy coordinate
                val=matrices_[mrow-1][mcol-1][2]# InsertY, (cell one up/left in Iy)
                print('InsertY')
        #Dealing with cases of equivalent origins, meaning matrix M is preferred over insertions
            if val == matrices_[mrow-1][mcol-1][0] + match_status and val == matrices_[mrow-1][mcol-1][1] + match_status:
                state='M'
                maxv=[mrow-1,mcol-1,0]
                val=matrices_[mrow-1][mcol-1][0]
                print('equiv.1')
            elif val == matrices_[mrow-1][mcol-1][0] + match_status and val == matrices_[mrow-1][mcol-1][2] + match_status:
                state='M'
                maxv=[mrow-1,mcol-1,0]
                val=matrices_[mrow-1][mcol-1][0]
                print('equiv.2')

    #Determine origin of Ix coordinate, 2 potential origins for Ix
        elif state == 'Ix':
            if val == matrices_[mrow-1][mcol][0] + g+h: 
                #[mrow-1][mcol][0] is M score, so 'state' reflects this
                state = 'M'
                maxv= [mrow-1,mcol,0]
                val=matrices_[mrow-1][mcol][0] # OpenGapX, (cell one left in M)+gap_ext+gap_opn
                print('OpenGapx')

            elif val == matrices_[mrow-1][mcol][1] + g:
                #[mrow-1][mcol][1] is Ix score, so 'state' reflects this
                state = 'Ix'
                maxv= [mrow-1,mcol,1]
                val=matrices_[mrow-1][mcol][1] # ExtendGapX, (cell one left in Ix)+gap_ext
                print('ExtendGapX')

    #Determine origin of Iy coordinate
        elif state == 'Iy':
        #2 Potential origins for Iy
            if val == matrices_[mrow][mcol-1][0] + g+h:# OpenGapY=(cell one up in M)+gap_ext+gap_opn
                state = 'M'
                maxv= [mrow,mcol-1,0]
                val=matrices_[mrow][mcol-1][0] 
                print('OpenGapY')

            elif val == matrices_[mrow][mcol-1][2] + g: # ExtendGapY= (cell one up in Iy)+gap_ext
                state = 'Iy'
                maxv= [mrow,mcol-1,2]
                val=matrices_[mrow][mcol-1][2]
                print('ExtendGapY')

        return maxv,state,val 


    #print out the visual alignment of letters and traceback of the best scoring alignment
    def print_traceback(matrices_):
        #this will print as expected with internal gaps
        print("Building traceback...")
        #print coordinates of the highest scoring cell as a list of two numbers
        maxv=get_max(matrices_) #part of the get_max() function is to print the 'max score:<>'
        print(f'starting score{maxv}') #print cooridnates of highest score

        #traverse the matrix to find the traceback elements
        #if more than one path just pick one
        topstring=""
        midstring=""
        bottomstring=""

        #here we need to store the first (old) position so we can track if it is an insertion of deletion
        old_maxv=maxv

        #add the alignment elements('-' for gaps)
        search=True
        #define 'state'based on number of highest value returned by get_max
        #state will be used to tell the traceback when gaps have been opened or extended
        if maxv[2]==0: 
            state= 'M'
        elif maxv[2]==1:
            state='Ix'
        elif maxv[2]==2:
            state='Iy'

        while(search):
            #NOTE if the coordinates of the highest score in the matrix give value 0, stop running traceback becasue the highest score is the top left cell
            maxv,state,highest=traceback_affine(matrices_,maxv,state)
               
            if(state=='M'):
                topstring+=seq1[maxv[1]]
                bottomstring +=seq2[maxv[0]]
            elif(state=='Ix'):
                topstring += "-"
                bottomstring += seq2[maxv[0]]
            elif(state=='Iy'):
                topstring += seq1[maxv[1]]
                bottomstring += "-"   

            #NOTE Add pipe if the letter in seq1 is the same as seq2 and this is not the previous value considered
            if (seq1[maxv[1]] == seq2[maxv[0]]) & (old_maxv[1]!=maxv[1]) & (old_maxv[0] != maxv[0]):
                midstring += "|"
            else:
                midstring += " "

            old_maxv = maxv

            if maxv[0]==0 and maxv[1]==0 : #if in the first row or first column, add gaps in seq in the other sequence until the rest fo the string is printed          
                search=False
                continue
            elif highest==0 and maxv[0]!=0 and maxv[1]!=0:
                topstring+=seq2[:maxv[0]][::-1]
                bottomstring+=seq1[:maxv[1]][::-1]
                search=False
                continue
                

        print(topstring[::-1])
        print(midstring[::-1])
        print(bottomstring[::-1])

    def print_matrix(matrices_):
        global rows, cols,seq1,seq2
        s1="  " +seq1
        s2=" "+seq2
        
        print(f"Dimensions: r= {rows}, c= {cols}")

        for a in range(0,cols+2):
            print(s1[a], end ="")
            print("\t\t", end ="")
        print("\n",end="")

        for i in range(0, rows+1):
            print(s2[i], end ="")
            print(" \t", end ="")
            for j in range(0, cols+1):
                print(f"\t{(matrices_[i][j])}",end="")
            print("\n",end="")

    # Use Biopython for real projects that read fasta!!!
    def read_fasta_filename(filename):
        seq = ""

        with open(filename, 'r') as filehandle:

            for line in filehandle:

                if(line[0]== ">"):
                    continue
                seq= seq+ line

            #delete any whitespaces & quotes
            import re
            pattern = re.compile(r'\s+')
            seq = re.sub(pattern, '', seq)

            seq = seq.replace('\"', '')

            return seq

#Load the sequences, converts fasta file to a string variable
    seq1 = read_fasta_filename(margs.seq1)
    seq2 = read_fasta_filename(margs.seq2)
    cols=len(seq1)
    rows=len(seq2)

    print("Sequence1: ",seq1);
    print("Sequence2: ",seq2);

    #create the 3 matrices, returns the 3
    matrices_=initialise_matrices() #NOTE makes empty matrix in the shape of sequences
    matrices_=build_matrix(matrices_) #NOTE calls score calculation for each position in the pre-built matrix
    print_matrix(matrices_)  #print out the (best scoring path from the) SW matrix
    print_traceback(matrices_) #NOTE converts scoring into the alignment by comparing ______
