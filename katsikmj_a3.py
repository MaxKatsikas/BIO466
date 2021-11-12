from sys import argv
import pandas as pd

def get_seq_fasta(File): ## Reads in FASTA file
    file_dict = {}
    fh = open(File, 'r')
    name, desc, seq = ('', '', '')
    for line in fh:
        line = line.strip()
        if line.startswith('>'):
            if seq != '':
                name, desc, seq = ('', '', '')
            name = line[1:line.index(' ')]
            desc = line[line.index(' ') + 1:]
            file_dict[name] = ''
        else:
            seq = seq+line
            file_dict[name] = seq
    return(file_dict)

def get_expression(lung_exp, prot_exp): ## Fills dictionaries for lung and prostate cancer
    dict_lung, dict_prot = {}, {}
    lung_df = pd.read_table(lung_exp)
    prot_df = pd.read_table(prot_exp)

    for i in lung_df.index: ## Fills dict_lung
        row = lung_df.iloc[i]
        numL = row[1]
        numC = row[2]
        if numL/numC >= 1.5:
            value = '+'
        elif numL/numC <= 0.667:
            value = '-'
        else:
            value = '.'
        ratio = '{0}:{1}:{2}'.format(numL, numC, value)
        dict_lung[row[0]] = ratio

    for i in prot_df.index: ## Fills dict_prot
        row = prot_df.iloc[i]
        numP = row[1]
        numC = row[2]
        if numP/numC >= 1.5:
            value = '+'
        elif numP/numC <= 0.667:
            value = '-'
        else:
            value = '.'
        ratio = '{0}:{1}:{2}'.format(numP, numC, value)
        dict_prot[row[0]] = ratio

    return (dict_lung, dict_prot)

def common_genes(dict_lung, dict_prot): ## Creates common set between cancer genes
    common_set = set()
    for i in dict_lung.keys():
        if i in dict_prot.keys():
            common_set.add(i)
    return(common_set)

def unique_genes(dict_lung, dict_prot): ## Creates sets of unique genes
    UniqueLung = set()
    UniqueProt = set()
    for i in dict_lung.keys():
        if i not in dict_prot.keys():
            UniqueLung.add(i)
    for i in dict_prot.keys():
        if i not in dict_lung.keys():
            UniqueProt.add(i)
    return(UniqueLung, UniqueProt)

def compare(dict_lung, dict_prot, common_set): ## Compares up and down regulation of genes
    print("4. Expression Comparison (lung_vs_control) vs (prostate_vs_control):")
    for i in common_set:
        if dict_lung[i][-1] == dict_prot[i][-1]:
            congruence = 'Congruence'
        else:
            congruence = 'Disparity'
        print("{0} NumL:NumC:Exp {1} NumP:NumC:Exp {2} {3}".format(i,
                                                                   dict_lung[i],
                                                                   dict_prot[i],
                                                                   congruence))

def get_nuc_pos_profile(dict, gene_set, start=0, end=140):
    dict_adj = {}
    for gene in gene_set: ## Creates dictionary with only values in selected gene set
        dict_adj[gene] = dict[gene]

    dict_a = {}
    dict_t = {}
    dict_g = {}
    dict_c = {}

    for i in range(start, end+1):
        dict_a[i] = 0
        dict_g[i] = 0
        dict_t[i] = 0
        dict_c[i] = 0

    for seq in dict_adj.values():
        for j in range(start, end+1):
            nucleotide = seq[j]
            if nucleotide == 'A':
                dict_a[j] += 1
            if nucleotide == 'G':
                dict_g[j] += 1
            if nucleotide == 'T':
                dict_t[j] += 1
            if nucleotide == 'C':
                dict_c[j] += 1

    for k in range(start, end, 10): ## Creates printout
        for l in range(0,10):
            position = k + l
            if l == 0:
                print("\nPosition\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}"
                      .format(position, position+1, position+2, position+3,
                              position+4, position+5, position+6, position+7,
                              position+8, position+9))
                print("--------\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--")
                print("A\t\t\t", end='')
            try:
                print(dict_a[position], end='\t')
            except:
                print("--\t", end='')
        print("\nG\t\t\t", end='')
        for l in range(0,10):
            position = k + l
            try:
                print(dict_g[position], end='\t')
            except:
                print("--\t", end='')
        print("\nT\t\t\t", end='')
        for l in range(0,10):
            position = k + l
            try:
                print(dict_t[position], end='\t')
            except:
                print("--\t", end='')
        print("\nC\t\t\t", end='')
        for l in range(0,10):
            position = k + l
            try:
                print(dict_c[position], end='\t')
            except:
                print("--\t", end='')
        print()

## Main
FASTA_name = argv[1]
lung_exp_name = argv[2]
prot_exp_name = argv[3]
file_dict = get_seq_fasta(FASTA_name) ## Passes FASTA_name into get_seq_fasta function
print("1: Files: FASTA ({0}) Cancer Expression Documents ({1}),({2})\n".format(FASTA_name, lung_exp_name, prot_exp_name))

## Passes in cancer expression files into get_expression function and filling dictionaries
dict_lung, dict_prot = get_expression(lung_exp_name, prot_exp_name)

print("2. Differentially expressed genes detected for lung cancer-control comparison:")
print("Gene ID \t NumL:NumC:Exp")
for i in dict_lung.keys():
    print(i, '\t', dict_lung[i])
print()
print("3. Differentially expressed genes detected for prostate cancer-control comparison:")
print("Gene ID \t NumP:NumC:Exp")
for i in dict_prot.keys():
    print(i, '\t', dict_prot[i])
print()

common_set = common_genes(dict_lung, dict_prot) ## Finds common set between cancer genes
UniqueLung, UniqueProt = unique_genes(dict_lung, dict_prot) ## Finds unique genes in each cancer type
compare(dict_lung, dict_prot, common_set) ## Passes in dictionaries and common set into compare function

print("\n5. Common DEGs (expressed genes) for both lung and prostate cancer:")
print(common_set)

print("\n6. DEGs uniquely expressed in lung cancer:")
print(UniqueLung)

print("\n7. DEGs uniquely expressed in prostate cancer:")
print(UniqueProt)
print()

set_choice = input("8. Which gene set do you want to examine (5|6|7) ")
print("\nYour answer is {0}".format(set_choice))
print("We will extract nucleotide position profile for the following sequences:")
if set_choice == '5':
    gene_set = common_set
    print(common_set)
elif set_choice == '6':
    gene_set = UniqueLung
    print(UniqueLung)
elif set_choice == '7':
    gene_set = UniqueProt
    print(UniqueProt)
print()

k = 0
while k == 0:
    try:
        set_range = input(
        "Please enter the start and end positions (e.g., 1-9; positions should be shorter than the gene length:140) ")
        range_split = set_range.split("-")
        start = range_split[0]
        end = range_split[1]
        start = int(start)
        end = int(end)
        if end - start > 0:
            if end - start < 140:
                get_nuc_pos_profile(file_dict, gene_set, start, end) ## Passes in data to get_nuc_pos_profile function
                k = 1
            else:
                print("Invalid entry: length must be shorter than 140")
        else:
            print("Invalid entry: end position must be higher than start position")
    except:
        print("Invalid Entry")