from sys import argv
import pandas as pd
import re
import SmithWaterman_class as sw

class MyGene(object):

    def __init__(self, File1Name, File2Name=''): ## Initializes class
        self.FastaFileName = File1Name
        self.File2Name = File2Name
        self.SeqDict = {} # A dictionary for gene ID and sequence
        self.SppDict = {} # A dictionary for gene ID and names of gene and species
        self.bs_list = ['TATTGTTTATT',
                        'TATTGTTTATT',
                        'TATTGTTTACT',
                        'AATTGTTTATT',
                        'AATTGTTTATT',
                        'CATTGTTTATT',
                        'TATTGTTTATT',
                        'GATTGTTTACT',
                        'AATTGTTTATT',
                        'AATTGTTTAGT'] # a list of default binding site
        self.Consensus = ''
        self.IUPAC = ''
        self.PFM = pd.DataFrame()
        self.PPM = pd.DataFrame()
        self.__extract_file()
        self.set_motif()
        self.set_PFM()
        self.set_PPM()
        self.set_Consensus()
        self.set_IUPAC()

    def __extract_file(self): ## Handles and extracts file into seq and species dictionaries
        fh = open(self.FastaFileName, 'r')
        file = fh.read()
        file_split = file.split('\n')
        gene_id = ''
        for line in file_split:
            if line.startswith('>'):
                info = line.split(' ')
                gene_id = info[0][1:]
                self.SeqDict[gene_id] = ""
                self.SppDict[gene_id] = (info[1], info[2])
            else:
                self.SeqDict[gene_id] += line

    def set_motif(self): ## Replaces default list with specified list
        try:
            fh = open(self.File2Name, 'r')
            self.bs_list = []
            for line in fh:
                line = line.strip('\n')
                self.bs_list.append(line)
        except:
            self.bs_list

    def set_Consensus(self): ## Finds concensus based on proportions
        for i in range(0, len(self.bs_list[0])):
            a_count, g_count, c_count, t_count = 0, 0, 0, 0
            for j in range(0, len(self.bs_list)):
                if self.bs_list[j][i] == 'A':
                    a_count += 1
                elif self.bs_list[j][i] == 'G':
                    g_count += 1
                elif self.bs_list[j][i] == 'C':
                    c_count += 1
                elif self.bs_list[j][i] == 'T':
                    t_count += 1
            col = {'A': a_count / len(self.bs_list), 'G': g_count / len(self.bs_list), 'C': c_count / len(self.bs_list),
                   'T': t_count / len(self.bs_list)}
            self.Consensus += max(col, key=col.get) ## Adds each max proportion to the consensus seq

    def set_IUPAC(self): ## Finds IUPAC code
        for i in range(0, len(self.bs_list[0])):
            a_count, g_count, c_count, t_count = 0, 0, 0, 0
            for j in range(0, len(self.bs_list)):
                if self.bs_list[j][i] == 'A':
                    a_count += 1
                elif self.bs_list[j][i] == 'G':
                    g_count += 1
                elif self.bs_list[j][i] == 'C':
                    c_count += 1
                elif self.bs_list[j][i] == 'T':
                    t_count += 1
            col = {'A': a_count / len(self.bs_list), 'G': g_count / len(self.bs_list), 'C': c_count / len(self.bs_list),
                   'T': t_count / len(self.bs_list)}
            max_col = max(col, key=col.get) ## Finds max nucleotide per column
            if col[max_col] > 0.5: ## Specifies if another nucleotide is possible
                self.IUPAC += max_col
            elif col['A'] + col['G'] == 1.0:
                self.IUPAC += 'R'
            elif col['C'] + col['T'] == 1.0:
                self.IUPAC += 'Y'
            elif col['G'] + col['C'] == 1.0:
                self.IUPAC += 'S'
            elif col['A'] + col['T'] == 1.0:
                self.IUPAC += 'W'
            elif col['G'] + col['T'] == 1.0:
                self.IUPAC += 'K'
            elif col['A'] + col['C'] == 1.0:
                self.IUPAC += 'M'
            else:
                self.IUPAC == 'N'

    def set_PFM(self): ## Creates frequency matrix
        bs_list_len = len(self.bs_list)
        list = self.bs_list
        for i in range(0, len(self.bs_list[0])):
            a_count, g_count, c_count, t_count = 0, 0, 0, 0
            for j in range(0, len(self.bs_list)):
                if self.bs_list[j][i] == 'A':
                    a_count += 1
                elif self.bs_list[j][i] == 'G':
                    g_count += 1
                elif self.bs_list[j][i] == 'C':
                    c_count += 1
                elif self.bs_list[j][i] == 'T':
                    t_count += 1
            col = (a_count, g_count, c_count, t_count)
            self.PFM["Position " + str(i + 1)] = col
        self.PFM.index = ['A', 'G', 'C', 'T']

    def set_PPM(self): ## Creates proportion matrix
        length = len(self.bs_list)
        for i in range(0, len(self.bs_list[0])):
            a_count, g_count, c_count, t_count = 0, 0, 0, 0
            for j in range(0, len(self.bs_list)):
                if self.bs_list[j][i] == 'A':
                    a_count += 1
                elif self.bs_list[j][i] == 'G':
                    g_count += 1
                elif self.bs_list[j][i] == 'C':
                    c_count += 1
                elif self.bs_list[j][i] == 'T':
                    t_count += 1
            ## Creates proportions
            col = (a_count/len(self.bs_list), g_count/len(self.bs_list), c_count/len(self.bs_list), t_count/len(self.bs_list))
            self.PPM["Position " + str(i + 1)] = col
        self.PPM.index = ['A', 'G', 'C', 'T']

    def find_bs(self): ## Prints info about potential matches for consensus
        for gene in self.SeqDict.keys():
            seq = self.SeqDict[gene]
            gene_name, species = self.SppDict[gene]
            start, end = self.scan_seq(seq)
            print("\nSequence Num:\t{0}".format(species[-2:]))
            print("GeneID:\t\t\t{0}".format(gene))
            print("GeneName:\t\t{0}".format(gene_name))
            print("SpeciesName:\t{0}".format(species))
            if start >= 2:
                ext_start = start - 2
            else:
                ext_start = start
            if end != len(seq):
                ext_end = end + 2
            else:
                ext_end = end
            ext_seq = seq[ext_start:ext_end]
            ## Calls align_seq method and passes out SW results
            seq1_aligned, seq2_aligned, ident_perc, gap_perc = self.align_seq(ext_seq, self.Consensus)
            print("Identity:\t\t{0:.1%}".format(ident_perc))
            print("Gaps:\t\t\t{0:.1%}\n".format(gap_perc))
            ## Passes data into get_alignment method to print alignment
            self.get_alignment(start, start+len(seq1_aligned), gene, seq1_aligned, seq2_aligned)
            seq_print = seq[:start] + seq[start:end].lower() + seq[:end]
            ## Prints spaced results
            printWithRuler(seq_print)

    def scan_seq(self, seq): ## Returns locations of most likely binding site
        prob_dict = {}
        for i in range(0, len(seq)-len(self.bs_list[0])):
            seq_frag = seq[i:i+len(self.bs_list[0])]
            frag_prob = 0
            for j in range(0, len(self.bs_list[0])):
                nuc = seq_frag[j]
                ppm = self.PPM.iloc[:, j]
                nuc_prob = ppm[nuc]
                frag_prob += nuc_prob
            prob_dict[i] = frag_prob
        max_start = max(prob_dict, key=prob_dict.get)
        return (max_start, max_start + len(self.bs_list[0]))

    def align_seq(self, ext_seq, motif): ## Calls Smith Waterman class to align sequences
        seq1_aligned, seq2_aligned, ident_perc, gap_perc = sw.Smith_Waterman(ext_seq, motif).give_final_result()
        return(seq1_aligned, seq2_aligned, ident_perc, gap_perc)

    def get_alignment(self, start, end, gene, seq, motif): ## Prints alignment based on SW returned sequences
        alignment = ''
        for i in range(0, len(motif)):
            if seq[i] == motif[i]:
                alignment += '|'
            elif seq[i] == '-' or motif[i] == '-':
                alignment += ' '
            else:
                alignment += ':'
        print("{0}\t\t{1}\t{2}\t{3}".format(gene, start, seq, end))
        print("\t\t\t\t{0}".format(alignment))
        print("BS_Site\t\t1\t{0}\t{1}\n".format(motif, end-start))

def printWithRuler(Sequence, Spacer=' '): #Prints the Sequence in Non-FASTA format
    NON_FASTA_LINE_NUM=100
    counter=1
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

if __name__== '__main__': ## Main function
    print("\n********************************************************************************************")
    print("Use this program to detect potential protein binding sites for a set of sequences.")
    print("********************************************************************************************\n")
    print("Program: katsikmj@a4.py")
    print("Developer: Maxwell Katsikas")
    gene_file = MyGene(argv[1])
    ## Prints amount of sequences in input file
    print("Input file: {0}\t(There are a total of {1} sequences)\n".format(argv[1], len(gene_file.SeqDict)))
    ## Prompts user to enter bs_list file
    input_file = input("Please enter the name of file that contains binding-site sequences (e.g., bs.txt)\n")
    ## Assigns class to my_gene
    my_gene = MyGene(argv[1], input_file)
    print()
    print("[1]. The protein binding sequences are:\n")
    ## Prints bs_list
    for i in range (1, len(my_gene.bs_list)):
        print(i, my_gene.bs_list[i])
    print('\nThe consensus sequence : {0}'.format(my_gene.Consensus))
    print('The IUPAC sequence : {0}'.format(my_gene.IUPAC))
    print('\n[2]. The Position Frequency Matrix for the protein binding sequences is:')
    print(my_gene.PFM)
    print('\n[3]. The Position Probability Matrix for the protein binding sequences is:')
    print(my_gene.PPM)
    print('\n[4]. Results of sequence scan')
    ## Calls find.bs method in MyGene class
    my_gene.find_bs()