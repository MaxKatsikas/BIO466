from sys import argv
import re

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
            print('  ', counter, end=' ')
        else:
            print(' ', counter, end=' ')
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
    counterTuple=(T_Count, C_Count, A_Count, G_Count, N_Count)
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
    diNucleotideCounter = {}
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


def codonProfile(Sequence):
    codon_dict = {}
    for a in range(0, len(Sequence)-3, 3):  ##Scans sequence for dinucleotides
        codon=Sequence[a:a+3]
        codon_dict[codon]= 0 ## Initiates all codon sequences
    for a in range(0, len(Sequence)-3, 3):
        codon=Sequence[a:a+3]
        codon_dict[codon]+=1
    return (codon_dict)


def printCodonProfile(codon_dict):
    codon_order = ['TTT','TCT','TAT','TGT','TTC','TCC','TAC','TGC','TTA','TCA','TAA','TGA','TTG','TCG','TAG','TGG',
                   'CTT','CCT','CAT','CGT','CTC','CCC','CAC','CGC','CTA','CCA','CAA','CGA','CTG','CCG','CAG','CGG',
                   'ATT','ACT','AAT','AGT','ATC','ACC','AAC','AGC','ATA','ACA','AAA','AGA','ATG','ACG','AAG','AGG',
                   'GTT','GCT','GAT','GGT','GTC','GCC','GAC','GGC','GTA','GCA','GAA','GGA','GTG','GCG','GAG','GGG',]
    print('              \t\t\t2nd')
    print('         --------------------------------')
    print("1st       T       C       A       G         3rd")
    print()
    print('T        ',end='')
    for i in range(0, 16 ,1):
        try:
            codon=codon_dict[codon_order[i]]
        except:
            codon=0 ##If there is not a key corresponding to this, the value becomes 0
        print('{0}={1} '.format(codon_order[i], repr(codon).rjust(3)), end='')
        if i==3:
            print('   T\n         ',end='')
        if i==7:
            print('   C\n         ',end='')
        if i==11:
            print('   A\n         ',end='')
        if i==15:
            print('   G\n         ',end='')
    print()
    print('C        ', end='')
    for i in range(16, 32, 1):
        try:
            codon = codon_dict[codon_order[i]]
        except:
            codon = 0
        print('{0}={1} '.format(codon_order[i], repr(codon).rjust(3)), end='')
        if i == 19:
            print('   T\n         ', end='')
        if i == 23:
            print('   C\n         ', end='')
        if i == 27:
            print('   A\n         ', end='')
        if i == 31:
            print('   G\n         ', end='')
    print()
    print('A        ', end='')
    for i in range(32, 48, 1):
        try:
            codon = codon_dict[codon_order[i]]
        except:
            codon = 0
        print('{0}={1} '.format(codon_order[i], repr(codon).rjust(3)), end='')
        if i == 35:
            print('   T\n         ', end='')
        if i == 39:
            print('   C\n         ', end='')
        if i == 43:
            print('   A\n         ', end='')
        if i == 47:
            print('   G\n         ', end='')
    print()
    print('G        ', end='')
    for i in range(48, 64, 1):
        try:
            codon = codon_dict[codon_order[i]]
        except:
            codon = 0
        print('{0}={1} '.format(codon_order[i], repr(codon).rjust(3)), end='')
        if i == 51:
            print('   T\n         ', end='')
        if i == 55:
            print('   C\n         ', end='')
        if i == 59:
            print('   A\n         ', end='')
        if i == 63:
            print('   G\n         ', end='')
    print()


def translation(Sequence): ##Translates a given sequence
    trans_seq={ 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
                'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
                'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
                'NNN':'', 'NNC':'','NCC':'','CCN':'', 'CNN':''}
    dnaseq = Sequence
    rnaseq = dnaseq.replace('T', 'U')
    protein = []
    for i in range(0, len(rnaseq) - 2, 3):
        protein.append(trans_seq[rnaseq[i:i + 3]])
    return ("".join(protein))


def ORF(Sequence):
    seq0=Sequence[0:]
    seq1=Sequence[1:]
    seq2=Sequence[2:]
    rev_seq=Sequence[::-1]
    seq3=rev_seq
    seq4=rev_seq[1:]
    seq5=rev_seq[2:]

    #trans_total = ['NONE'] * 6
    trans1 = translation(seq0)
    #trans_total[0] = trans1
    trans2 = translation(seq1)
    #trans_total[1] = trans2
    trans3 = translation(seq2)
    #trans_total[2] = trans3
    trans4 = translation(seq3)
    #trans_total[3] = trans4
    trans5 = translation(seq4)
    #trans_total[4] = trans5
    trans6 = translation(seq5)
    #trans_total[5] = trans6

    ORF1 = re.findall('M\w{9,}\*', str(trans1)) ##Finds all sequences starting with M and ending with *
    ORF2 = re.findall('M\w{9,}\*', str(trans2))
    ORF3 = re.findall('M\w{9,}\*', str(trans3))
    ORF4 = re.findall('M\w{9,}\*', str(trans4))
    ORF5 = re.findall('M\w{9,}\*', str(trans5))
    ORF6 = re.findall('M\w{9,}\*', str(trans6))

    print("(2.5) 6-Frame Translations")
    print()
    #print(trans1, trans2, trans3)
    #print(ORF1, ORF2, ORF3, ORF4, ORF5, ORF6)

    if ORF1 != []:
        ORF1 = sorted(ORF1, reverse=True)
        for i in ORF1:
            print("ORF (+1)", i)
    else:
        print('ORF (+1) No Valid ORF')
    if ORF2 != []:
        ORF2 = sorted(ORF2, reverse=True)
        for i in ORF2:
            print("ORF (+2)", i)
    else:
        print('ORF (+2) No Valid ORF')
    if ORF3 != []:
        ORF3 = sorted(ORF3, reverse=True)
        for i in ORF3:
            print("ORF (+3)", i)
    else:
        print('ORF (+3) No Valid ORF')
    if ORF4 != []:
        ORF4 = sorted(ORF4, reverse=True)
        for i in ORF4:
            print("ORF (-1)", i)
    else:
        print('ORF (-1) No Valid ORF')
    if ORF5 != []:
        ORF5 = sorted(ORF5, reverse=True)
        for i in ORF5:
            print("ORF (-2)", i)
    else:
        print('ORF (-2) No Valid ORF')
    if ORF6 != []:
        ORF6 = sorted(ORF6, reverse=True)
        for i in ORF6:
            print("ORF (-3)", i)
    else:
        print('ORF (-3) No Valid ORF')


def inquiry(Sequence):
    print(" (3.1) Extract a DNA fragment", "\n", "(3.2) Compare two codon profiles", "\n", "(3.3) Find a motif", "\n")
    x = 0
    while x == 0:
        inquiry_choice = input("Your choice is: ") ##Prompts user for inquiry choice
        print()
        if inquiry_choice == "3.1":
            print("You can extract a DNA fragment from the given sequence.")
            a = 0
            while (a == 0):  ##Continues to prompt user until satisfied
                question_three = input('Please enter the start and end positions (e.g., 19::48):')
                try:  ##Starts over if invalid format
                    start_end = question_three.split("::")
                    start = start_end[0]
                    end = start_end[1]
                    int_start = int(start)
                    int_end = int(end)
                    difference = int_end - int_start
                    if int_start < int_end:  ##Checks for correct input
                        if int_start > 0:
                            adj_seq = Sequence[int_start:int_end]
                            print("The fragment you selected has a length of [{0}] nucleotides:".format(difference))
                            print("<{0}".format(int_start),
                                  "-" * (difference - (4 + len(start) + len(end))),
                                  "{0}>".format(int_end),
                                  end="")
                            print()
                            print('|' * difference)
                            print(adj_seq)
                            print()

                            counterTuple = nucleotideCounter(adj_seq)  ##Initiates nucleotideCounter function
                            print("T=", counterTuple[0], "C=", counterTuple[1], "A=", counterTuple[2], "G=",
                                  counterTuple[3],
                                  "N=", counterTuple[4], )
                            print("GC content=", gcContent(adj_seq), "%")

                            print("Di-nucleotide Counts:")
                            diNucleotideCounter = diCounter(adj_seq)  ##Initializes diCounter function
                            inquiry2(diNucleotideCounter)
                            print()

                            print("(2.4) Codon Profile:")
                            codon_dict = codonProfile(adj_seq)
                            printCodonProfile(codon_dict)

                            print()
                            ORF(adj_seq)
                            print()

                            question_four = input("Do you want to extract another DNA fragment (Y|N) ? ")
                            if question_four != 'Y':
                                a = 1
                        else:
                            print("Invalid Positions or Format!")
                    else:
                        print("Invalid Positions or Format!")
                except:
                    print("Invalid Positions or Format!")
            x = 1

        if inquiry_choice == "3.2":
            print("Please open another sequence file that contains only one sequence.")
            file2 = input("Enter the file name here: ") ##Prompts user for file within directory
            seq1=Sequence
            file1=file_name
            seq2_name, seq2desc, seq2=readFASTA(file2)
            codonProfileCompare(seq1, file1, seq2, file2)
            x = 1

        if inquiry_choice == "3.3":
            print("Please enter the motif that you want to find (e.g., AATTC,GCCTTA).")
            motif_choice = input("Enter the motif sequence here: ") ##Prompts user for target motif
            print()
            motif_dict = motif_finder(Sequence, motif_choice) ##Calls motif_finder function
            locations=motif_dict['motif_locations']
            locations = locations.split(',')
            print("Motif        Frequency    Position(Start-End)")
            try:
                print('{0}'.format(motif_choice).ljust(12), motif_dict['motif_frequencies'], '           ', end='')
                end_point=int(locations[0]) + len(motif_choice)
                print('( {0} - {1} )'.format(locations[0], end_point), end='')
                for m in range(1, len(locations), 1):
                    end_point = int(locations[m]) + len(motif_choice)
                    print('({0} - {1} )'.format(locations[m], end_point), end='')
                print()

            except:
                print("")
            seq_lowercase = re.sub(motif_choice, motif_choice.lower(), Sequence)
            print()
            spacer=' '
            printWithRuler(seq_lowercase, spacer)
            x = 1


def inquiry2(diNucelotideCounter): ##Used just to print dinucleotide counters for the Inquiry function
    print("TT=[{0}] TC=[{1}] TA=[{2}] TG=[{3}]".format(diNucleotideCounter['TT'],
                                                       diNucleotideCounter['TC'],
                                                       diNucleotideCounter['TA'],
                                                       diNucleotideCounter['TG']))
    print("CT=[{0}] CC=[{1}] CA=[{2}] CG=[{3}]".format(diNucleotideCounter['CT'],
                                                       diNucleotideCounter['CC'],
                                                       diNucleotideCounter['CA'],
                                                       diNucleotideCounter['CG']))
    print("AT=[{0}] AC=[{1}] AA=[{2}] AG=[{3}]".format(diNucleotideCounter['AT'],
                                                       diNucleotideCounter['AC'],
                                                       diNucleotideCounter['AA'],
                                                       diNucleotideCounter['AG']))
    print("GT=[{0}] GC=[{1}] GA=[{2}] GG=[{3}]".format(diNucleotideCounter['GT'],
                                                       diNucleotideCounter['GC'],
                                                       diNucleotideCounter['GA'],
                                                       diNucleotideCounter['GG']))


def codonProfileCompare(Seq1, File1, Seq2, File2):
    print("File1: {0} \nFile2: {1}".format(File1, File2))
    profile1 = codonProfile(Seq1)
    profile2 = codonProfile(Seq2)
    profile1_adj={}
    profile2_adj={}
    codon_order = ['TTT', 'TCT', 'TAT', 'TGT', 'TTC', 'TCC', 'TAC', 'TGC', 'TTA', 'TCA', 'TAA', 'TGA', 'TTG', 'TCG','TAG', 'TGG',
                   'CTT', 'CCT', 'CAT', 'CGT', 'CTC', 'CCC', 'CAC', 'CGC', 'CTA', 'CCA', 'CAA', 'CGA', 'CTG', 'CCG','CAG', 'CGG',
                   'ATT', 'ACT', 'AAT', 'AGT', 'ATC', 'ACC', 'AAC', 'AGC', 'ATA', 'ACA', 'AAA', 'AGA', 'ATG', 'ACG','AAG', 'AGG',
                   'GTT', 'GCT', 'GAT', 'GGT', 'GTC', 'GCC', 'GAC', 'GGC', 'GTA', 'GCA', 'GAA', 'GGA', 'GTG', 'GCG','GAG', 'GGG', ]
    for k in range(0, len(codon_order),1): ##Fills secondary profile with all 64 codons
        try:
            profile1_adj[codon_order[k]]=profile1[codon_order[k]]
        except:
            profile1_adj[codon_order[k]]=0
        try:
            profile2_adj[codon_order[k]]=profile2[codon_order[k]]
        except:
            profile2_adj[codon_order[k]]=0
    codon_compare = {}
    codon_compare2 = {}
    for a in range(0, len(codon_order), 1): ##Establishes comparison keys
        codon_compare[codon_order[a]]=profile1_adj[codon_order[a]]-profile2_adj[codon_order[a]]
        if -10<codon_compare[codon_order[a]]<10:
            codon_compare2[codon_order[a]]=''
        if -20<codon_compare[codon_order[a]]<=-10:
            codon_compare2[codon_order[a]] = "<"
        if -20<codon_compare[codon_order[a]]<=-10:
            codon_compare2[codon_order[a]] = "<"
        if -30<codon_compare[codon_order[a]]<=-20:
            codon_compare2[codon_order[a]] = "<<"
        if codon_compare[codon_order[a]]<=-30:
            codon_compare2[codon_order[a]] = "<<<"
        if 10<=codon_compare[codon_order[a]]<20:
            codon_compare2[codon_order[a]] = ">"
        if 20<=codon_compare[codon_order[a]]<30:
            codon_compare2[codon_order[a]] = ">>"
        if codon_compare[codon_order[a]]>=30:
            codon_compare2[codon_order[a]] = ">>>"
    ##print(profile1_adj)
    ##print(profile2_adj)
    ##print(codon_compare)
    print('              \t\t\t2nd')
    print('         --------------------------------')
    print("1st       T       C       A       G         3rd")
    print()
    print('T        ', end='')
    for i in range(0, 16 ,1):
        try:
            codon=codon_compare2[codon_order[i]]
        except:
            codon=0
        print('{0}={1} '.format(codon_order[i], repr(codon).rjust(3)), end='')
        if i==3:
            print('   T\n         ',end='')
        if i==7:
            print('   C\n         ',end='')
        if i==11:
            print('   A\n         ',end='')
        if i==15:
            print('   G\n         ',end='')
    print()
    print('C        ', end='')
    for i in range(16, 32, 1):
        try:
            codon = codon_compare2[codon_order[i]]
        except:
            codon = 0
        print('{0}={1} '.format(codon_order[i], repr(codon).rjust(3)), end='')
        if i == 19:
            print('   T\n         ', end='')
        if i == 23:
            print('   C\n         ', end='')
        if i == 27:
            print('   A\n         ', end='')
        if i == 31:
            print('   G\n         ', end='')
    print()
    print('A        ', end='')
    for i in range(32, 48, 1):
        try:
            codon = codon_compare2[codon_order[i]]
        except:
            codon = 0
        print('{0}={1} '.format(codon_order[i], repr(codon).rjust(3)), end='')
        if i == 35:
            print('   T\n         ', end='')
        if i == 39:
            print('   C\n         ', end='')
        if i == 43:
            print('   A\n         ', end='')
        if i == 47:
            print('   G\n         ', end='')
    print()
    print('G        ', end='')
    for i in range(48, 64, 1):
        try:
            codon = codon_compare2[codon_order[i]]
        except:
            codon = 0
        print('{0}={1} '.format(codon_order[i], repr(codon).rjust(3)), end='')
        if i == 51:
            print('   T\n         ', end='')
        if i == 55:
            print('   C\n         ', end='')
        if i == 59:
            print('   A\n         ', end='')
        if i == 63:
            print('   G\n         ', end='')
    print()


def motif_finder(Sequence, motif_choice): ##Finds locations and frequencies of a motif within a sequence
    motif_dict = {}
    motif_frequencies='motif_frequencies'
    motif_locations='motif_locations'
    locations=[]
    for i in range(0, len(Sequence), 1):
        if Sequence[i:i+len(motif_choice)] == motif_choice:
            locations.append(i)
    motif_dict[motif_frequencies]=Sequence.count(motif_choice)
    str_locations = str(locations)
    motif_dict[motif_locations]=str_locations[1:-1]
    return motif_dict


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
print()
print("(2.1) Nucleotide Counts:")
counterTuple = nucleotideCounter(Sequence) ##Initiates nucleotideCounter function
print("T=", counterTuple[0],"C=", counterTuple[1],"A=", counterTuple[2],"G=", counterTuple[3],"N=", counterTuple[4],)
print()
print("(2.2) GC content: ", gcContent(Sequence),"%")
print()
print("(2.3) Di-nucleotide Counts:")
diNucleotideCounter=diCounter(Sequence) ##Initiates diCounter function
print("TT=[{0}] TC=[{1}] TA=[{2}] TG=[{3}]".format(diNucleotideCounter['TT'],
                                                   diNucleotideCounter['TC'],
                                                   diNucleotideCounter['TA'],
                                                   diNucleotideCounter['TG'],))
print("CT=[{0}] CC=[{1}] CA=[{2}] CG=[{3}]".format(diNucleotideCounter['CT'],
                                                   diNucleotideCounter['CC'],
                                                   diNucleotideCounter['CA'],
                                                   diNucleotideCounter['CG'],))
print("AT=[{0}] AC=[{1}] AA=[{2}] AG=[{3}]".format(diNucleotideCounter['AT'],
                                                   diNucleotideCounter['AC'],
                                                   diNucleotideCounter['AA'],
                                                   diNucleotideCounter['AG'],))
print("GT=[{0}] GC=[{1}] GA=[{2}] GG=[{3}]".format(diNucleotideCounter['GT'],
                                                   diNucleotideCounter['GC'],
                                                   diNucleotideCounter['GA'],
                                                   diNucleotideCounter['GG'],))
print()

print("(2.4) Codon Profile:")
codon_dict = codonProfile(Sequence) ##Calls codonProfile function
printCodonProfile(codon_dict) ##Calls printCodonProfile function

ORF(Sequence) ##Calls ORF function


print()
print("[Part III]: Inquiry Mode")
print()
inquiry(Sequence) ##Calls inquiry function