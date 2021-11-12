from sys import argv

def readFASTQ(filename):
    print("INside readFASTA{0}".format(filename))
    fh = open(filename, 'r')
    lines=fh.read()
    line=lines.split('\n')
    print("{0}".format(line[-1]))
    print(len(line))
#    line_number=0
#    for line in fh:
#        line_number=line_number+1
#        line=line.strip()
#        print("Inside readFASTA{0}".format(line))
#
#    file_list=fh.split()
#    for i in range(0, line_number,4):
#        print("line number {0}".format(i))

file=argv[1]
print("{0".format(file))
readFASTQ(file)