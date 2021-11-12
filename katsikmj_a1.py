from sys import argv

program_name=argv[0]
file_name=argv[1]

def readFASTA(file_name): ##Reads in the file
    file_handler=open(file_name)
    seq_file= file_handler.read()
    defLine = seq_file.split('\n')[0][0:]
    seq_name = defLine[0:defLine.index(' ')]
    desc = defLine[len(seq_name):]
    sequence = ''.join(seq_file.split('\n')[1:])
    return seq_name, desc, sequence

def printInFASTA(SeqName, Sequence, SeqDescription): ##Prints Sequence in FASTA format
    FASTA_LINE_NUM=60
    print("{0} {1}".format(SeqName, SeqDescription))
    for i in range(0, len(Sequence), FASTA_LINE_NUM):
        seq_60 = Sequence[i:i+FASTA_LINE_NUM]
        for j in range(0, len(seq_60), 10):
            seq_10 = seq_60[j:j+10]
            print(seq_10, end='')
        print()


def printWithRuler(Sequence, Spacer): #Prints the Sequence in Non-FASTA format
    NON_FASTA_LINE_NUM=100
    counter=0
    counterTwo=0
    print('              ', end='')
    for i in range(0, 100, 10):
        counterTwo = counterTwo + 1
        print(counterTwo, '       ',Spacer, end='')
    print()
    print('Line ', end='')
    for h in range(0,10,1):
        print('1234567890', end=Spacer)
    print()
    for i in range(0, len(Sequence), NON_FASTA_LINE_NUM):
        seq_100 = Sequence[i:i + NON_FASTA_LINE_NUM]
        if counter<10:
            print(counter, '   ', end='')
        else:
            print(counter, '  ',end='')
        counter = counter + 1
        for j in range(0, len(seq_100), 10):
            seq_10 = seq_100[j:j + 10]
            print(seq_10, end=Spacer)
        print()


def nucleotideCounter(Sequence): ##Calculates Nucleotides from a Sequence
    A_Count = 0
    G_Count = 0
    C_Count = 0
    T_Count = 0
    N_Count = 0
    for a in range(0, len(Sequence), 1):
        if Sequence[a]=='A':
            A_Count=A_Count+1
        elif Sequence[a]=='G':
            G_Count=G_Count+1
        elif Sequence[a]=='C':
            C_Count=C_Count+1
        elif Sequence[a]=='T':
            T_Count=T_Count+1
        elif Sequence[a]=='N':
            N_Count=N_Count+1
    counterTuple=(A_Count, G_Count, C_Count, T_Count, N_Count)
    return(counterTuple)


def gcContent(Sequence): ##Calculates GC Content
    G_Count=0
    C_Count=0
    for i in range(0, len(Sequence), 1):
        if Sequence[i]=='G':
            G_Count=G_Count+1
        elif Sequence[i]=='C':
            C_Count=C_Count+1
    GC_Content=(G_Count + C_Count)/len(Sequence)
    float(GC_Content)
    GC_Content = GC_Content * 100 ##Converts to %
    GC_Content = "{:.2f}".format(GC_Content)
    return(GC_Content)


def diCounter(Sequence):
    diNucleotideCounter = {'AA':0,'AT':0, 'AG':0, 'AC':0, ##Establishes dictionary
                           'TA':0,'TT':0,'TG':0,'TC':0,
                           'GA':0,'GT':0,'GG':0,'GC':0,
                           'CA':0,'CT':0,'CG':0,'CC':0}
    AA_counter=0
    AT_counter=0
    AG_counter=0
    AC_counter=0
    TA_counter = 0
    TT_counter = 0
    TG_counter = 0
    TC_counter = 0
    GA_counter=0
    GT_counter=0
    GG_counter=0
    GC_counter=0
    CA_counter=0
    CT_counter=0
    CG_counter=0
    CC_counter=0
    for a in range(0,len(Sequence),1): ##Scans sequence for dinucleotides
        if Sequence[a:a+2] == 'AA':
            AA_counter=AA_counter+1
            diNucleotideCounter['AA'] = AA_counter
        if Sequence[a:a+2] == 'AT':
            AT_counter=AT_counter+1
            diNucleotideCounter['AT'] = AT_counter
        if Sequence[a:a+2] == 'AG':
            AG_counter=AG_counter+1
            diNucleotideCounter['AG'] = AG_counter
        if Sequence[a:a+2] == 'AC':
            AC_counter=AC_counter+1
            diNucleotideCounter['AC'] = AC_counter
        if Sequence[a:a+2] == 'TA':
            TA_counter=TA_counter+1
            diNucleotideCounter['TA'] = TA_counter
        if Sequence[a:a+2] == 'TT':
            TT_counter=TT_counter+1
            diNucleotideCounter['TT'] = TT_counter
        if Sequence[a:a+2] == 'TG':
            TG_counter=TG_counter+1
            diNucleotideCounter['TG'] = TG_counter
        if Sequence[a:a+2] == 'TC':
            TC_counter=TC_counter+1
            diNucleotideCounter['TC'] = TC_counter
        if Sequence[a:a+2] == 'GA':
            GA_counter=GA_counter+1
            diNucleotideCounter['GA'] = GA_counter
        if Sequence[a:a+2] == 'GT':
            GT_counter=GT_counter+1
            diNucleotideCounter['GT'] = GT_counter
        if Sequence[a:a+2] == 'GG':
            GG_counter=GG_counter+1
            diNucleotideCounter['GG'] = GG_counter
        if Sequence[a:a+2] == 'GC':
            GC_counter=GC_counter+1
            diNucleotideCounter['GC'] = GC_counter
        if Sequence[a:a+2] == 'CA':
            CA_counter=CA_counter+1
            diNucleotideCounter['CA'] = CA_counter
        if Sequence[a:a+2] == 'CT':
            CT_counter=CT_counter+1
            diNucleotideCounter['CT'] = CT_counter
        if Sequence[a:a+2] == 'CG':
            CG_counter=CG_counter+1
            diNucleotideCounter['CG'] = CG_counter
        if Sequence[a:a+2] == 'CC':
            CC_counter=CC_counter+1
            diNucleotideCounter['CC'] = CC_counter
    return (diNucleotideCounter)


def inquiry(Sequence):
    a=0
    while (a==0): ##Continues to prompt user until satisfied
        question_three = input('Please enter the start and end positions (e.g., 19::48):')
        try: ##Starts over if invalid format
            start_end = question_three.split("::")
            start=start_end[0]
            end=start_end[1]
            int_start=int(start)
            int_end=int(end)
            difference=int_end - int_start
            if int_start<int_end: ##Checks for correct input
                if int_start>0:
                    adj_seq=Sequence[int_start:int_end]
                    print("The fragment you selected has a length of [{0}] nucleotides:".format(difference))
                    print("<{0}".format(int_start),
                          "-"*(difference-(4+len(start)+len(end))),
                          "{0}>".format(int_end),
                          end="")
                    print()
                    print('|'*difference)
                    print(adj_seq)
                    print()

                    counterTuple = nucleotideCounter(adj_seq) ##Initiates nucleotideCounter function
                    print("A=", counterTuple[0], "G=", counterTuple[1], "C=", counterTuple[2], "T=", counterTuple[3],
                          "N=", counterTuple[4], )
                    print("GC content=", gcContent(adj_seq), "%")

                    print("Di-nucleotide Counts:")
                    diNucleotideCounter = diCounter(adj_seq) ##Initializes diCounter function
                    print("AA=[{0}] AT=[{1}] AG=[{2}] AC=[{3}]".format(diNucleotideCounter['AA'],
                                                                       diNucleotideCounter['AT'],
                                                                       diNucleotideCounter['AG'],
                                                                       diNucleotideCounter['AC'], ))
                    print("TA=[{0}] TT=[{1}] TG=[{2}] TC=[{3}]".format(diNucleotideCounter['TA'],
                                                                       diNucleotideCounter['TT'],
                                                                       diNucleotideCounter['TG'],
                                                                       diNucleotideCounter['TC'], ))
                    print("GA=[{0}] GT=[{1}] GG=[{2}] GC=[{3}]".format(diNucleotideCounter['GA'],
                                                                       diNucleotideCounter['GT'],
                                                                       diNucleotideCounter['GG'],
                                                                       diNucleotideCounter['GC'], ))
                    print("CA=[{0}] CT=[{1}] CG=[{2}] CC=[{3}]".format(diNucleotideCounter['CA'],
                                                                       diNucleotideCounter['CT'],
                                                                       diNucleotideCounter['CG'],
                                                                       diNucleotideCounter['CC'], ))
                    print()
                    question_four=input("Do you want to extract another DNA fragment (Y|N) ? ")
                    if question_four!='Y':
                        a=1

                else:
                    print("Invalid Positions or Format!")
            else:
                print("Invalid Positions or Format!")
        except:
            print("Invalid Positions or Format!")



## Main
print("Welcome to sequence viewer! Programmer: katsikmj")
print("There is", 1, "sequence detected in the file:",file_name)

print("[Part I]: Display mode")

question_one=input("Question 1: Do you want to view the sequence in FASTA format or not (Y|N) ? ")
print()

SeqName, SeqDescription, Sequence = readFASTA(file_name) ##Initiates readFASTA function and reads in variables

if question_one == "Y":
    printInFASTA(SeqName, Sequence, SeqDescription) ##Initiates pirntInFASTA function

else:
    question_two=input("Question 2: Do you need a spacer for viewing nucleotide positions (Y|N) ? ")
    print("{0} {1}".format(SeqName, SeqDescription))
    if question_two == "Y":
        Spacer=' '
    else:
        Spacer=''

    printWithRuler(Sequence, Spacer) ##Initiates printWithRuler function

print()
print("[Part II]: Analysis Mode")
print("Nucleotide Counts:")
counterTuple = nucleotideCounter(Sequence) ##Initiates nucleotideCounter function
print("A=", counterTuple[0],"G=", counterTuple[1],"C=", counterTuple[2],"T=", counterTuple[3],"N=", counterTuple[4],)
print("GC content=", gcContent(Sequence),"%")

print("Di-nucleotide Counts:")
diNucleotideCounter=diCounter(Sequence) ##Initiates diCounter function
print("AA=[{0}] AT=[{1}] AG=[{2}] AC=[{3}]".format(diNucleotideCounter['AA'],
                                                   diNucleotideCounter['AT'],
                                                   diNucleotideCounter['AG'],
                                                   diNucleotideCounter['AC'],))
print("TA=[{0}] TT=[{1}] TG=[{2}] TC=[{3}]".format(diNucleotideCounter['TA'],
                                                   diNucleotideCounter['TT'],
                                                   diNucleotideCounter['TG'],
                                                   diNucleotideCounter['TC'],))
print("GA=[{0}] GT=[{1}] GG=[{2}] GC=[{3}]".format(diNucleotideCounter['GA'],
                                                   diNucleotideCounter['GT'],
                                                   diNucleotideCounter['GG'],
                                                   diNucleotideCounter['GC'],))
print("CA=[{0}] CT=[{1}] CG=[{2}] CC=[{3}]".format(diNucleotideCounter['CA'],
                                                   diNucleotideCounter['CT'],
                                                   diNucleotideCounter['CG'],
                                                   diNucleotideCounter['CC'],))
print()

print("[Part III]: Inquiry Mode")
print("You can extract a DNA fragment from the given sequence.")

inquiry(Sequence) ##initiates inquiry function