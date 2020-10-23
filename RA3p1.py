'''

create fasta file for mtDNA and for Y chromosome

Do the MSA, maybe save it into another file?
For MSA we need to choose gap penalty model along with gap penalty values and 
subs matrix
0. Find the two seqs that give the highest alignment score
1. look at each two seqs (e.g. a and b), calculate pairwise alignment score
2.total alignment score is then a sum of all pairwise alig.scores
- so e.g. if we have seqs a,b,c,d, then the sum S= Sab + Sac + Sad +Sbc + Sbd + 
Scd 
'''
#parsing the files into mtDNA fasta and Ychr fasta. Atm Im parsing them into 
#text files as I cant open fastas on my laptop.
with open('GeneticData.txt', 'r') as genefile, open('mtDNA.txt','w') as outputdna, open('Ychr.txt','w') as outputY:
    #genefile=genefile.readlines()
    mtDNA_dict={}
    seq_listDNA=[]
    Ychr_dict={}
    seq_listY=[]
    
   # name=next(genefile)
    
    seqmtDNA=''
    seqY=''
    for lines in genefile:
        #print(lines)
        if lines.strip():
            #print(lines)    #strip the empty lines (\n)
            
            if lines.startswith("mtDNA"):
                seqmtDNA=next(genefile).rstrip()
                
            elif lines.startswith("Y"):
                seqY=next(genefile).rstrip()
            elif "hemophilia" in lines:
                continue
            else:           #by this point we are on the last line of the indicidual
                name=lines
                if seqmtDNA:    #true as lomg as the line isnt empty line
                    outputdna.write('>' + name + seqmtDNA +'\n')
                   # name='' #name emptied
                if seqY:
                    outputY.write('>' + name + seqY + '\n')

#get seq1 and compare its every NT to the ones in the equivalent pos in seq2,3,4...
#do scoring based on this
#we have the dictionary with the names as keys and values as seqs, no align. make single dics with one key, one value

with open("mtDNA.txt", 'r') as DNAfile, open('Ychr.txt','r') as Ychrfile, open("scores.txt","w"):
    
    DNA_dict={}
    lista1 = []
    idheader = None
    for lines in DNAfile:
        if lines.startswith('>'):
            if idheader:
                DNA_dict[idheader]=''.join(lista1) # adding the id to the dic, making the whole 3 lines into a one line, then deleting the lista content so we can start again
                del lista1[:]

            idheader = lines.strip().replace('>','')
        else:
            lista1.append(lines.strip())
            #print(lista1)
    DNA_dict[idheader]=''.join(lista1)       # need to add this once more
                                            #as the last id+seq ends with a seq
                                            #instead of the next id
    del lista1[:]
    #print(DNA_dict)

    Ychrom_dict={}
    lista2 = []
    idheader2 = None
    for lines in Ychrfile:
        if lines.startswith('>'):
            if idheader2:
                Ychrom_dict[idheader2]=''.join(lista2) # adding the id to the dic, making the whole 3 lines into a one line, then deleting the lista content so we can start again
                del lista2[:]

            idheader2 = lines.strip().replace('>','')
        else:
            lista2.append(lines.strip())
            #print(lista)
    Ychrom_dict[idheader2]=''.join(lista2)       # need to add this once more
                                            #as the last id+seq ends with a seq
                                            #instead of the next id
    del lista2[:]
    #print(DNA_dict)
    
    def Scores_mtDNA(DNA_dict):
        transition = ['AG', 'TC', 'GA', 'CT'] #the transition scores
        NTscoreTotal=0
        Seqscore_list=[]  #list for the total scores of NTs when comparing two seqs
        Identical_NTs=0   #for counting the identical ones once we start off the loops
       
        checked_pair=[]
        
        for key1, seq1 in DNA_dict.items():
            for key2, seq2 in DNA_dict.items():
                if key1==key2:
                    continue
                if [key1,key2] not in checked_pair:     #so we only add the key1-key2 comb if its not already in the set, avoiding repeats
                    checked_pair.append([key1,key2])
                Identical_NTs=0
    
                for NTa,NTb in zip(seq1, seq2): #parallel iteration-in one loop calculates the overall score for the alignment
    
                    if NTa==NTb:
                        if '-' in NTa and '-' in NTb:
                            score=0
                            #NTscoreTotal.append(score)
                            NTscoreTotal +=score
                            
                        else:
                            score=1
                            NTscoreTotal +=score
                            Identical_NTs+=1   #keep count of identical NTs, identical gaps dont count
                    
                    else: #if the nucleotides are not the same
                        if '-' in NTa or '-' in NTb: #if one sequence hNTas NTa gNTap when aligned to the other one the score is -1 
                            score = -1
                            NTscoreTotal +=score
                        
                        else: #if one of them is not a -, it creates a variable nt_sum that is the sum of the two nucleotides and compares them to the transition list
                            NT_both = NTa + NTb
                            if NT_both in transition: # transition:  score is -1
                                score = -1
                                NTscoreTotal +=score
                            else: #if not transition, then transversion, which has a score of -2 that is added to NTscoreTotal
                                score = -2
                                NTscoreTotal +=score
                    perc_identity=Identical_NTs/len(seq1)*100
                    
                    Header_Seq=str(key1) +'\t'+str(key2) + '\t' + str(perc_identity) + '\t' +str(NTscoreTotal)   #a string, works as the header id, i.e. the two person's seqs that are compared
                    Seqscore_list.append(Header_Seq)
                    NTscoreTotal=0
        #print(Seqscore_list)
        return Seqscore_list

    def Scores_Ychr(Ychrom_dict):
        transition = ['AG', 'TC', 'GA', 'CT'] #the transition scores
        NTscoreTotal=0
        Seqscore_list=[]  #list for the total scores of NTs when comparing two seqs
        Identical_NTs=0   #for counting the identical ones once we start off the loops
       
        checked_pair=[]
        
        for key1, seq1 in DNA_dict.items():
            for key2, seq2 in DNA_dict.items():
                if key1==key2:
                    continue
                if [key1,key2] not in checked_pair:     #so we only add the key1-key2 comb if its not already in the set, avoiding repeats
                    checked_pair.append([key1,key2])
                Identical_NTs=0
    
                for NTa,NTb in zip(seq1, seq2): #parallel iteration-in one loop calculates the overall score for the alignment
    
                    if NTa==NTb:
                        if '-' in NTa and '-' in NTb:
                            score=0
                            #NTscoreTotal.append(score)
                            NTscoreTotal +=score
                            
                        else:
                            score=1
                            NTscoreTotal +=score
                            Identical_NTs+=1   #keep count of identical NTs, identical gaps dont count
                    
                    else: #if the nucleotides are not the same
                        if '-' in NTa or '-' in NTb: #if one sequence hNTas NTa gNTap when aligned to the other one the score is -1 
                            score = -1
                            NTscoreTotal +=score
                        
                        else: #if one of them is not a -, it creates a variable nt_sum that is the sum of the two nucleotides and compares them to the transition list
                            NT_both = NTa + NTb
                            if NT_both in transition: # transition:  score is -1
                                score = -1
                                NTscoreTotal +=score
                            else: #if not transition, then transversion, which has a score of -2 that is added to NTscoreTotal
                                score = -2
                                NTscoreTotal +=score
                    perc_identity=Identical_NTs/len(seq1)*100
                    
                    Header_Seq=str(key1) +'\t'+str(key2) + '\t' + str(perc_identity) + '\t' +str(NTscoreTotal)   #a string, works as the header id, i.e. the two person's seqs that are compared
                    Seqscore_list.append(Header_Seq)
                    NTscoreTotal=0
        #print(Seqscore_list)
        return Seqscore_list
    
    mtDNA_results=Scores_mtDNA(DNA_dict)
    print(mtDNA_results)
    Ychrom_results=Scores_Ychr(Ychrom_dict)
    #print(mtDNA_results)
    #print(Ychrom_results)

    def output_dict(my_list, output_file):
        with open(output_file, "w") as out: #open a given file to write to
            label_set=set()
            for line in my_list: # unpack dictionary
                    names, the_rest=line.split(':')
                    label_set.add(names)
                    
                    if names not in label_set:
                        bothnames = names.replace('-','\t')
                        print(bothnames, the_rest.replace('','\t'), sep="\t")
                    
                    else:
                        continue
               # print(names)

                #print(label_set)
                #scores, A=the_rest.split()
               # print(A)
                #Header="SampleA \t SampleB \t IdentityScore \t Score"
               # perc_ident, A.split()
                #Names=line.split(':')[0]
               # label_sets.add(Names)
                #Score=line.split('')[0]
                #print(label_sets)
                
                #print(f"{}", sep="\t", file=out)

    
    output_dict(mtDNA_results, "output_file_mtDNA.txt")
    output_dict(Ychrom_results, "output_file_Ychr.txt")
    
