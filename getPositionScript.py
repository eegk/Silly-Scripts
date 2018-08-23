# -*- coding: utf-8 -*-
import sys 

my_positions = open(sys.argv[1], 'r')
my_fasta = open(sys.argv[2], 'r')

with open(sys.argv[1], 'r') as f:    
    lines = len(f.readlines())

names= []
seqs = []

def read_positions(my_positions):
            
        line = my_positions.read()            
        values = line.split()         
        return values
            
def read_fasta(my_fasta):

        name, seq = None, []   
        
        for line in my_fasta:                        

            line = line.rstrip()
            
            if line.startswith(">"):

                if name: yield (name, ''.join(seq))
                    
                name, seq = line, []                            

            else:
                
                seq.append(line)
        
        if name: yield (name, ''.join(seq))                

for name, seq in read_fasta(my_fasta):
    
        seqs.append(seq)
        names.append(name)


if __name__ == '__main__':
    output_file=open(sys.argv[3],"w+")    
    a=0
    b=a+1
    positions = read_positions(my_positions)

    for i in range(lines):    

            single_seqs = seqs [i]        
            #print names[i]+'_my_gene_name', '\n', single_seqs[int(positions[a])-1:int(positions[b])]
            print >> output_file, names[i]+'_my_gene_name', '\n', single_seqs[int(positions[a])-1:int(positions[b])]
    
            a = b+1       
            b = a+1    
            
    output_file.close()
    


