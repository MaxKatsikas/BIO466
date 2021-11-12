def readFile(file_name):
    name_list=[]
    desc_list=[]
    seq_list=[]
    fh = open(file_name, 'r')
    name, desc, seq=('','','')
    for line in fh:
        line = line.strip()
        if line.startswith('>'):
            if seq !='':
                print("name=[{0}] desc=[{1}] seqbase=[{2}]".format(name, desc, seq))
                name = line[1:line.index(' ')]
                desc = line[line.index(' ')+1:]
        else:
            seq=seq+line
        return(name_list, desc_list, seq_list)