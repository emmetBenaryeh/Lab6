from SequenceAnalysis import FastAreader
import sys 
class TRNA():
    """
    Class to represent a tRNA sequence
    """
    list = [] #global list variable will contain all instances of the class 
    def powerset(self,s):
        """
        Return the powerset of a given string
        """
        result = set()
        n = len(s)
        for i in range(n):
            for j in range(i, n):
                result.add(s[i:j+1])
        return result
    def __init__(self,header,seq):
        """
        Initializes an instance of the class with a name seqeunce, amino acid name, and a powerset and adds itself to the list
        """
        self.name = header.replace(" ","")
        self.aa = self.getAAname()
        self.seq = seq
        self.powerSet = self.powerset(seq)
        self.list.append(self)
    def getseq(self):
        """
        return seqeunce 
        """
        return self.seq
    def uniqueSeq(self):
        """
        Returns essential unique seqeunces
        """
        localList = self.list.copy()
        localList.remove(self) #list of all tRNA objects not including itself
        set1 = set()
        for index,val in enumerate(localList):
            set1 = set1.union(val.powerSet) #union of all powersets 
        set2 = list(self.powerSet - set1) #list of unique substrings
        set3 = set2.copy()
        for index, val in enumerate(set2): #loop compares each value in the list to each other value and deletes values containing other values until only essentials remain
            for i in range(0,len(set2)):
                if val != set2[i] and val in set2[i]:
                    try:
                        set3.remove(set2[i])
                    except ValueError:
                        continue
        self.clearList() #clear list to prevent screw ups down the line 
        return set3
    def clearList(self):
        """
        Clears global list variable
        """
        self.list = []
    def find_all_positions(self,string, substring):
        """
        Returns all positions where a given substring is found in a parent string
        """
        start = 0
        positions = []
        while True:
            index = string.find(substring, start)
            if index == -1:
                break
            positions.append(index)
            start = index + 1
        return positions
    def getAAname(self):
        """
        Returns the name of the amino acid associated with a tRNA
        """
        return1 = self.name.split('|')[1]
        return return1
    def printer(self):
        """
        Prints all essential unique sequences alligned with their position in the sequence
        """
        returnStr = ""
        out = self.uniqueSeq()
        printer = []
        sequence = self.getseq()
        for item in out: 
            pos = self.find_all_positions(sequence,item) #find all positions assoicated with an essential 
            for position in pos:
                printer.append([str(item),position])
        printer = sorted(printer, key=lambda x: x[1]) #sort by position in the sequence
        sys.stdout.reconfigure(encoding='utf-8') #handles special characters
        print(self.name)
        print(sequence)
        for item in printer:
            print('.'*item[1]+str(item[0])) #print with number of periods that matches position to allign with the sequence

def main(infile=None):
    #if infile is None:
        #infile = sys.argv[1]
    trna_list = []
    #with open(infile, 'r', encoding='utf-8') as f:
    reader = FastAreader(infile)
    for head,seq in reader.readFasta():
        #stripping allignment characters and white space
        seq = seq.replace("-","")
        seq = seq.replace("_","")
        seq = seq.replace(".","")
        trna_list.append([head,seq])
    """
    I deserve whatever I get for this monstrous list of tRNA instantiation but it seemed like the easiest way at the time and now i just feel if it aint broke... 
    """
    Trna1 = TRNA(trna_list[0][0], trna_list[0][1])
    Trna2 = TRNA(trna_list[1][0], trna_list[1][1])
    Trna3 = TRNA(trna_list[2][0], trna_list[2][1])
    Trna4 = TRNA(trna_list[3][0], trna_list[3][1])
    Trna5 = TRNA(trna_list[4][0], trna_list[4][1])
    Trna6 = TRNA(trna_list[5][0], trna_list[5][1])
    Trna7 = TRNA(trna_list[6][0], trna_list[6][1])
    Trna8 = TRNA(trna_list[7][0], trna_list[7][1])
    Trna9 = TRNA(trna_list[8][0], trna_list[8][1])
    Trna10 = TRNA(trna_list[9][0], trna_list[9][1])
    Trna11 = TRNA(trna_list[10][0], trna_list[10][1])
    Trna12 = TRNA(trna_list[11][0], trna_list[11][1])
    Trna13 = TRNA(trna_list[12][0], trna_list[12][1])
    Trna14 = TRNA(trna_list[13][0], trna_list[13][1])
    Trna15 = TRNA(trna_list[14][0], trna_list[14][1])
    Trna16 = TRNA(trna_list[15][0], trna_list[15][1])
    Trna17 = TRNA(trna_list[16][0], trna_list[16][1])
    Trna18 = TRNA(trna_list[17][0], trna_list[17][1])
    Trna19 = TRNA(trna_list[18][0], trna_list[18][1])
    Trna20 = TRNA(trna_list[19][0], trna_list[19][1])
    Trna21 = TRNA(trna_list[20][0], trna_list[20][1])
    Trna22 = TRNA(trna_list[21][0], trna_list[21][1])
    trnas = [Trna1, Trna2, Trna3, Trna4, Trna5, Trna6, Trna7, Trna8, Trna9, Trna10, Trna11, Trna12, Trna13, Trna14, Trna15, Trna16, Trna17, Trna18, Trna19, Trna20, Trna21, Trna22]
    trnas = sorted(trnas, key=lambda x: x.aa) #sorting tRNAs alpha by amino acid 
    sys.stdout.reconfigure(encoding='utf-8')#handles special characters
    for item in trnas: #printing all essentials for each tRNA in alpha order
        item.printer()
    Trna1.clearList
if __name__ == "__main__":
    sys.stdin.reconfigure(encoding='utf-8')#handles special characters
    main()