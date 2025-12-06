#This script will extract information alpha helices and beta sheets from the raw file from chimera
def get_file(name, out):
    f = open(name, 'r')
    data = f.readlines()
    f.close()
    data = ','.join(list(set([x.split()[2][:-2] for x in data])))
    newn = name.split('/')[-1].replace('.txt','')+'_list.txt'
    newf = open(out + newn,'w')
    newf.write(data)
    newf.close()

def main():
    infolder = '../data/input/'
    ofolder = '../data/output/'
    get_file(infolder + 'helices.txt', ofolder)
    get_file(infolder + 'sheets.txt', ofolder)

if __name__ == "__main__":
    main()
