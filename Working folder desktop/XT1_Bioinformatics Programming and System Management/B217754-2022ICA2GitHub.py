#!/usr/bin/python3

#ICA2 Protein Query Programme (P.Q.P)
#Written by: B217754-2022
#Purpose:Retrieve sequences for a protein family of interest, perform preliminary conservation assessment, motif identification, amino acid property tagging.

#############################################################################
#############################################################################
################                HOUSEKEEPING               ##################
#############################################################################
#############################################################################

#Importing the required modules
import os, subprocess, pandas as pd, time

#Defining a Quality of Life function which delays output, giving the user time to read the text above before status outputs for the next section go to screen.
def ReadingTime(x): 
    for i in range(x): #creates a line of ..... and waits for 0.75 seconds, repeated x number of times.
        print(".....")
        time.sleep(0.75)


def PQP(): #The programme's code is enveloped in a function name, so that it can be called by the user later
    investigation_counter=+1 #Counter is additive on each loop of the PQP function (when the user asks to run another query at the end of the programme). This creates uniquely numbered directoryies for each run of PQP.
    WorkingDir=os.getcwd() #Retrieve the present OS working directory and storing in a variable for python as a string.
    print(WorkingDir) #Show the user the present working directory.
    #Make folders for work and change to main parent folder
    subprocess.call(f"mkdir {WorkingDir}/ProteinInvestigation{investigation_counter}", shell=True) #Create a directory for all files resuling from this work to be held in. 
    print(f"\nChanging working directory to folder: {WorkingDir}/ProteinInvestigation{investigation_counter}\n")
    os.chdir(f"{WorkingDir}/ProteinInvestigation{investigation_counter}") #Changing the OS working directory to the newly-created folder.


    #############################################################################
    #############################################################################
    ################                    P.Q.P                  ##################
    #############################################################################
    #############################################################################
    
    #############################################################################
    ################        Define user's Query terms           #################
    #############################################################################
    queryDetails={} #To create the dictionary for query terms to be held in.
    queryDetailstr=None #To print the query terms ina  nice format for the user.

    def qsf(): #The Query Searching Function


        queryDetails['ProteinFam']=input("Please type the protein family of interest. An example is pyruvate dehydrogenase\n")
        queryDetails['TaxGroup']=input("Please type the taxnomic group of interest. An example is Ascomycete fungi\n")
        queryDetailstr=" ".join(queryDetails.values())
        
        if queryDetails["ProteinFam"] == "" or queryDetails["TaxGroup"] == "":
            print(f"Your query of interest is:\n{queryDetailstr}")
            print("Your query must contain a protein name and a taxonomic group. Please try again")
            qsf()
        else:
            print(f"The query terms you entered are:\n{queryDetailstr}.\n\nContinuing...\n")
    
    dirname=os.getcwd() 
    print(f"The cwd is:{dirname}")#Show the user the new working directory
    ReadingTime(5) #Print 5 lines of dots, and wait about 4 seconds total
    qsf() #Call the Query Searching Function
    queryDetailstr=" ".join(queryDetails.values())

    def cqpf(): #the Change Query Function
    #Give user the option to change query with, also if query is empty then state 'not allwoed' and return to top.
        ChangeQueryPrompt=input(f"Your query is {queryDetailstr}, enter P to Proceed or C to Change your query\n").upper()  
        if ChangeQueryPrompt=="P":
            print(f"Moving onto sequence retrieval with {queryDetailstr}")
            ReadingTime(5)
        elif ChangeQueryPrompt=="C":
            qsf(), cqpf() #call back functionsup to this point to restart from query specification
        else:
            print("That wasn't an option")
            cqpf()
    cqpf() #Calling the Change Query Function

    #############################################################################
    ################        Sequence retrieval section          #################
    #############################################################################

    ##EDirect searching##
    #Creating the variables to make a query string for the Edirect querying commands esearch and efetch.
    edirInBase=f'esearch -db protein -query "{queryDetails["ProteinFam"]} [PROT] AND {queryDetails["TaxGroup"]} [ORGN]"' #query string stem 
    edirInCount="xtract -pattern ENTREZ_DIRECT -element Count" # The createsa small summary of the results in the XML format, including a count of results. 
    edirInSpecies="efetch -format docsum | xtract -pattern DocumentSummary -element Organism| sort | uniq -c" #This extracts the result summary for each entry in the docsum format (a human-readable format). It extracts the species name from each entry, sorts these alphabetically, then removes duplicates, giving the number of duplicates in-line with each name.
    edirInFasta="efetch -format fasta" #This extracts the protein sequences for each protein entry in fasta format
    edirInGbGen="efetch -format gb" #This retrieves the results in genbank format for navigation by the user and for the commands later in PQP
    

    #This command gives the number of entries for this query.... 
    print(f"Counting the number of entries on NCBI Protein for your query: {queryDetailstr}")
    edirInCountString=f"{edirInBase} | {edirInCount}" #assembling the query for the number of entries
    edirInCountNum=(os.popen(edirInCountString).read()).rstrip('\n') #passing the edirect query to the Linux OS
    print(f"There are {edirInCountNum} entries on NCBI Protein for {queryDetails['ProteinFam']} in {queryDetails['TaxGroup']}")
    edirInCountNum=int(edirInCountNum) #turn variable into integer to allow numerical comparison operators to assess the number of results against cuttoffs of 0 and 1000.


    def cpf(): #Module which redirects the user to restart the search or if their initial search query returns either 0 or >1000 results.
        if CountPrompt=="Y":
            print("Restarting")
            ReadingTime(5)
            qsf(), cqpf(),cpf() #call back functions up to this point to restart from query specification
        elif CountPrompt=="N":
            print("That concludes our work in the PQP.\nThank you for using our sevices!\nQuitting PQP")
            ReadingTime(5)
            exit()
        else:
            print("That wasn't an option")
            cpf()

    #....if count>0 say count and list of species,  if count=0, tell them and give option to change query
    if edirInCountNum <1:
        CountPrompt=input("Since there are no results from this query, you can restart PQP to enter new query terms. Restart PQP? (Y/N)\n").upper()
        cpf()
    elif edirInCountNum >=1000:
        CountPrompt=input(f"The search returned {edirInCountNum} results, you can still view the species list, however, alignment and further investgation will not be possible.\n\nPlease refine your search to produce fewer results (a more specific taxonomic group query may be easiest way to do so).").upper()
        cpf()
    else:
        print("Moving onto species retrieval")
        ReadingTime(5)          

    #This command gives the names of species which have an entry for the query protein on the NCBI protein db, and the number of entries per species in the first column. It removes multiple lines for the same species
    edirInSpeciesListString=f"{edirInBase} | {edirInSpecies} | tee speciesList.txt" #assembling the query for the species list, outputting to a file too
    edirInSpeciesList=os.popen(edirInSpeciesListString).read() #passing the edirect query to the Linux OS, saving as variable to display for user.
    UniqueSpeciesCount=edirInSpeciesList.count("\n") #counting how many lines there are to calculate the number of species


    print(f'''The following list shows the species with a recorded {queryDetails['ProteinFam']} protein.
    The numbers in the first column represent the number of entries a species has in NCBI Protein with your query
    ({UniqueSpeciesCount} unique species in total for this search):''')
    ReadingTime(10)
    print(f"{edirInSpeciesList}")

    if edirInCountNum>1000:
        print("Over 1000 results returned, no further processing allowed. Exiting PQP")
        ReadingTime(5)
        exit()
    
    def vslfp(): #the View Species List Finction. Choice to open the .txt file of the names of species or change query.
        try:
            speciesListPrompt=input(f"The list of species has also been saved as file speciesList.txt. \n\tWould you like to view it or change your query? (Y/N/C)\n\t").upper()
            #Evaluating the user's input, must be Y,N OR C
            if speciesListPrompt=="Y":
                print("Opening the file viewer in a new window. This will take some time, and the programme will progress in the meantime.")
                subprocess.call('xdg-open "speciesList.txt"&', shell=True) 
                ReadingTime(5) 
            elif speciesListPrompt=="N":
                print("OK, the image can be found later and is saved as speciesList.txt. Onto the next section")
            elif speciesListPrompt=="C":
                print("Restarting") #reseting variables and recalling modules to return to this point afresh
                ReadingTime(5)
                queryDetails={'ProteinFam':'','TaxGroup':''}
                queryDetailstr=" "
                qsf()
                ChangeQueryPrompt=""
                cqpf()
                CountPrompt=""
                cpf(), vslfp()
            else: #restricting input to a specifed set of responses
                print("Please type Y,N or C")
                speciesListPrompt=""
                vslfp()
                
        except:
            print("Something went wrong with the looping, please re-run and specify a different query")
    vslfp()

    #To allow user to browse information on query hits, create gb files for each:
    edirInGbGenString=f"{edirInBase}|{edirInGbGen} | tee queryGBs.txt" #generating the genbank reports for each entry to a file
    print("Retrieving Genbank reports of the query results, please wait.") 
    edInGb=os.popen(edirInGbGenString).read()
    print("A file with the query results in GenBank format has been saved as queryGBs.txt for your viewing.")

    GBs=edInGb.split("//\n")#split the Gb sequence file into a list object
    GBFileCounter=0
    ReadingTime(3)
    print(f"the first entry in GBs is\n\n {GBs[0]}\n\n")
    for entry in GBs: #for each Gb sequence in the newly-generated list, write the sequence to a Gb file which is numbered using the above counter variable.
        GBFileCounter+=1
        with open(f"GBs{GBFileCounter}.txt","w") as my_file:
            my_file.write(entry)
            print(f"Writing Gb file {GBFileCounter}")

    print("Genbanks files generated for later.\n")
    




    #############################################################################
    ################          Fasta generation section          #################
    #############################################################################

    print("Making separate files of each sequence in Fasta format for use later. The full list of fastas can be found in the file queryseqs.fasta")
    ReadingTime(5)

    #This command makes a file of the query results fastas
    edirInFullFastaString=f"{edirInBase} | {edirInFasta} | tee queryseqs.fasta" #will retrieve fastas and write to a variable and a fileos.popen("chmod 755 queryseqs.fasta")
    edirFastas=os.popen(f"{edirInFullFastaString}").read() #contains the fastas

    # To Extract FASTAs to their own files for use in prosite later:
    # Extract FASTAs to a list
    Fastas=edirFastas.split('\n>') #split the fasta sequence file into a list object
    FilterLengthFastas=[]        
    
    def mmpf(): #the Minimum Maximum Function. Allow user to input min amd max values to filter the length of sequences, use this number directly on fasta seqs file.
        try:
            MinPrompt,MaxPrompt=input("Please enter the desired minimum and maximum length of sequences\n").split() #taking two values in a single prompt
            if MinPrompt>MaxPrompt:
                print("Ensure that your value for Minimum length is smaller than that for Maximum length!\n")
                print(MinPrompt,MaxPrompt)
                mmpf()
            else:
                for seq in Fastas:
                    if int(len(seq)) in range (int(MinPrompt),int(MaxPrompt)):
                        FilterLengthFastas.append(seq)
                
                with open(f"filteredQuerySeqs.fasta","w") as my_file: #Making master list of filtered fastas
                    for items in FilterLengthFastas:
                        my_file.writelines(f">{items}\n\n")
                        

                FastaFileCounter=0
                for entry in Fastas: #for each fasta sequence in the newly-generated list, write the sequence to a fasta file which is numbered using the above counter variable. The if conditional should avoid the addition of an additinoal '>' int he first file.
                    FastaFileCounter+=1
                    with open(f"Fasta{FastaFileCounter}.fasta","w") as my_file:
                        if {FastaFileCounter}==1:
                            my_file.write(entry)#The first fasta sequence retains its '>' character when the file is split, so it does not need to be added when writing the files.
                        else:
                            my_file.write(">"+entry)# fasta file was split using the '>' character ath the start of each sequence, so this needs to be added back
                        print(f"Writing Fasta file {FastaFileCounter}")
                

                print("Fastas generated for later.\n")
                #ReadingTime(5)            

        except ValueError: #If the user enters anything other than a number, python should give an error and the loop will repeat back to the prompt.
            print("Something was wrong with your input, please type two integer numbers using digits, separated by a space.")
            mmpf()

    def Vfp(): #the View Fastas function, used to view the new file of fasta sequences
        ViewfastasPrompt=input("Do you want to view the file of fastas?(Y/N):\n").upper()
        if ViewfastasPrompt=="Y":
            print("Opening the file viewer in a new window. This will take some time, and P.Q.P will progress in the meantime.")
            subprocess.call('xdg-open "filteredQuerySeqs.fasta"&', shell=True)
            print("OK, moving on Moving onto clustering for conservation")
            ReadingTime(3)
        elif ViewfastasPrompt=="N":
            print("OK, moving on Moving onto clustering for conservation")

    mmpf() #Calling the Min/Max function
    Vfp() #Calling the view fastas function

    #############################################################################
    ################            Sequence Clustering              ################
    #############################################################################

    print("Clustering sequences using clustalo with 30 threads.")

    clustering="clustalo --threads 30 -i queryseqs.fasta -o queryalign.fasta --auto --force"
    os.popen(clustering,'w') #generates aligned sequences fasta file
    subprocess.call(clustering,shell=True) #generates aligned sequences fasta file
    ReadingTime(5)
    
    os.chmod(f"{dirname}/queryalign.fasta",775) # Changing permissions to the generated file to allow PQP to alter it. 
    subprocess.call("chmod 775 queryalign.fasta",shell=True) # Changing permissions to the generated file to allow PQP to alter it. 


    subprocess.call("showalign -sequence queryalign.fasta -outfile prettyQueryAlignment.txt -sformat1 fasta -width 100 -auto -warning -die", shell=True) #Using show align to generate a nicely formatted sequence alignment for the user to access
    ReadingTime(5)

    def vaf(): #the View Alignment function
        ViewAlignPrompt=input("A formatted alignment, with a consensus sequence, has been saved as 'prettyQueryAlignment.txt' .\n\tWould you like to view it? (Y/N)\n\t").upper()
        if ViewAlignPrompt=="Y":
            print("Opening the file viewer in a new window. This will take some time, and P.Q.P will progress in the meantime.\n If the file is not created, please re-run the whole programme as this seems to rectify the issue.")
            subprocess.call('xdg-open "prettyQueryAlignment.txt"&', shell=True)
            ReadingTime(5)
        elif ViewAlignPrompt=="N":
            print("OK, the image can be found later and is saved as prettyQueryAlignment.txt. Onto the next section")
        else:
            print("Please type Y or N")
            ViewAlignPrompt=""
            vaf()
    vaf() #Calling the View Alignment function


    print("Moving onto Plotcon conservation alignment.")
    ReadingTime(5)

    
    windowSize=50 #Creating a window size variable to allow the function to write user's input into command.
    def wspf(): #the Window Size Function
        windowSizePrompt=input("The default window size for plotcon conservation alignment has been set to 50. \nIf you would like to use a different window size (lower may be more noisy, higher may lose accuracy), please type the number now in digits\n")
        global windowSize #Using the previously defined variable for windowSize
        try:
            windowSizePrompt=int(windowSizePrompt)
            if windowSizePrompt<150:
                windowSize=windowSizePrompt
            else:
                print("!!Please enter a number below 150!!")
                windowSizePrompt=None
                wspf()
        except ValueError:
            print("!!!Please enter a number below 150 using digits!!!")
            wspf()

    def vcppf(): #the View Conservation Plot Function
        conservationPlotPrompt=input(f"Conservation plot created from your sequences. \n\tDo you want to view the image (sorry it's rotated!)? (Y/N)\n\t").upper()
        if conservationPlotPrompt=="Y":
            print("Attempting to open the plot in new window. To resume investigation: close the image window, return to this terminal and press enter")
            subprocess.call("gs conservationPlot.ps", shell=True)
        elif conservationPlotPrompt=="N":
            print("OK, the image can be found later and is saved as conservationPlot.ps. Onto the next section")
        else:
            print("Please type Y or N")
            conservationPlotPrompt=""
            vcppf()
    

    def plotconUse(): #Prompts user for desired window size, makes a corresponding plot then gives optoin to re-run with an adjusted window size.
        wspf()
        os.system(f"plotcon -sequence queryalign.fasta -winsize {windowSize} -graph ps -goutfile conservationPlot") #Plotcon for conservation across species, uses user input to decide window size, 
        vcppf()

        #Allow the user to repeate the process but change the window size for the generation of a different plot.
        changewindowPrompt=input("If you want to make a new plot with a different window size, please enter C. Otherwise, enter N to move on\n")
        if changewindowPrompt=="C":
            plotconUse()
        elif changewindowPrompt=="N":
            print("OK, moving on")
            ReadingTime(3)
        else:
            print("please enter C or N")
            ReadingTime(2)
        

    plotconUse()#Using the plotcon feature of PQP, defined above
    ReadingTime(5)
    
    #############################################################################
    ################            MOTIF identification           ##################
    #############################################################################
    
    print("Moving onto patmatmotifs search")
    ReadingTime(5)

    # iterating through the fasta files to make motif files
    for seq in os.listdir(dirname):
        if seq.endswith(".fasta"):
            seqName=seq.rstrip(".fasta")
            print(f"Scanning {seqName}")
            os.popen(f"patmatmotifs -sequence {seq} -outfile {seqName}.patmatmotifs -sprotein1 Yes -auto -sformat1 fasta -rformat excel")
        else:
            continue

    #Making the df 
    df=pd.DataFrame(columns=['SeqName','Start','End', 'Score','Strand','Motif'])
    
    print(f"Iterating through detected motif hits and adding them to a table")
    for seq in os.listdir(dirname):
        if seq.endswith(".patmatmotifs"):
            intermediarydf=pd.read_csv(seq, sep="\t")
            df=pd.concat([df,intermediarydf]) #changed to concat from .append()
        else:
            continue

    #Results presentation
    df=df.set_index('SeqName') #to remove the ugly indexing numbers from appended results
    #Changing some of the column names
    df.rename(columns = {'SeqName':'SeqAccession', 'Motif':'MotifName'}, inplace = True)
    df=df.sort_values('SeqName')

    #Exporting results to a file
    df.to_csv(f"{queryDetailstr}MotifResults.tsv",sep="\t", header=True, )
    print(f"{df}\n") #The results table variable is also printed for the user to preview.
    
    #Choice to open excel file of motifs
    def vmf():
        ViewMotifsPrompt=input(f"The detected motifs have been gathered in a dataframe (previewed above) and sorted by SeqAccession alphabetical order. Do you want to open this file for viewing? (This may take a moment) \n\t(Y/N)\n\t").upper()
        if ViewMotifsPrompt=="Y":
            print("Opening the file viewer in a new window. This will take some time, and P.Q.P will progress in the meantime.\n")
            subprocess.call(f'xdg-open "{queryDetailstr}MotifResults.tsv"&', shell=True)
            ReadingTime(5)
        elif ViewMotifsPrompt=="N":
            print(f"OK, the image can be found later and is saved as {queryDetailstr}MotifResults.tsv. Onto the next section")
        else:
            print("Please type Y or N")
            ViewMotifsPrompt=""
            vmf()
    vmf()


    #############################################################################
    ################    Wildcard Section: Amino Acid properties     #############
    #############################################################################
    

    for GbFile in os.listdir(dirname): # iterating through the Gb files to make property files
        if GbFile.startswith("GBs"):
            outputName=GbFile.rstrip(".txt")
            subprocess.call(f"pepinfo {GbFile} -sformat gb -graph ps -goutfile {outputName}graphs -outfile {outputName}pepinfo.ps",shell=True)
        else:
            pass


    def pepinfoView():
        pepinfochoose=input("The amino acid properties of each sequence have been saved to image files.\n To view the image of a particular sequence, type the index number of the sequence in question.\n Otherwise, type N to carry on\n")
        if pepinfochoose!="N":
            try:
                pepinfochoose=int(pepinfochoose) #convert the input into an integer to identify the file of interest
                subprocess.call(f"gs GBs{pepinfochoose}graphs.ps", shell=True)
                anotherpepView=input("If you want to view another pepinfo graph, press Y, if not, type any other character and P.Q.P will return to welcome prompt.").upper()
                if anotherpepView=="Y":
                    pepinfoView()
                else:
                    print("That concludes our work in the PQP.\nThank you for using our sevices!")
            except:
                print(f"The sequence index must be an integer between 0 and {edirInCountNum}")
                pepinfoView()
        else:
            print("That concludes our work in the PQP.\nThank you for using our sevices!")
            ReadingTime(4)

    pepinfoView()
    
    ReadingTime(4)
    PQPactivation()

def PQPactivation(): #Giving the user a choice to run P.Q.P or not.
    StartPrompt=input("Welcome to the Protein Query Programme (P.Q.P.), would you like to run the programme (Y/N) (N exits programme)?\n").upper()
    if StartPrompt=="Y":
        print("\nRunning PQP")
        ReadingTime(2)
        PQP()
    elif StartPrompt=="N":
        print('''Understandable, here's a nice quote then:\n\tKnowledge, like air, is vital to life.\n\tLike air, no one should be denied it\n\t\t-Alan Moore, V for Vendetta''')
        print("Exiting script")
        ReadingTime(5)
        exit()
    else:
        print("\nChoose Y or N, please.\n")
        PQPactivation()

PQPactivation()

response=None
while response!="Continue":
    print(str(range(1,10)))
    print("sequence not chosen")
    response=input("type Continue to move on")

print(f"Sequence chosen")
