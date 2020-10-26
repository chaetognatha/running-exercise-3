#parsing the files into mtDNA fasta and Ychr fasta. Atm Im parsing them into
#text files as I cant open fastas on my laptop.
with open('GeneticData.txt', 'r') as genefile, open('mtDNA.txt','w') as outputdna, open('Ychr.txt','w') as outputY:
    #genefile=genefile.readlines()
  #  mtDNA_dict={}
   # seq_listDNA=[]
    #Ychr_dict={}
   # seq_listY=[]
    seqmtDNA=''   #define string for mtDNA seq
    seqY=''     ##define string for mtDNA seq
#create separate files for mtDNA and Ychr with the name included
    for lines in genefile:
        #print(lines)
        if lines.strip():
            #print(lines)    #strip the empty lines (\n)

            if lines.startswith("mtDNA"):     #identify mtDNA, the seq is in the next line
                seqmtDNA=next(genefile).rstrip()    #save the seq line into the variable, remove line changes

            elif lines.startswith("Y"): #identify Y chromosome and repeat the process
                seqY=next(genefile).rstrip()
            elif "hemophilia" in lines:     #if the line contains this word, then skip it and its seq
                continue
            else:           #by this point we are on the last line of the individual (e.g. Rasputin)
                name=lines
                if seqmtDNA:    #true as lomg as the line isnt empty line
                    outputdna.write('>' + name + seqmtDNA +'\n')
                   # name='' #name emptied
                    seqmtDNA=''
                if seqY:
                    outputY.write('>' + name + seqY + '\n')
                    seqY=''
#get seq1 and compare its every NT to the ones in the equivalent pos in seq2,3,4...
#do scoring based on this
#we have the dictionary with the names as keys and values as seqs, no align. make single dics with one key, one value

with open("mtDNA.txt", 'r') as DNAfile, open('Ychr.txt','r') as Ychrfile:
    #open both files and make them into dictionary
    DNA_dict={}
    lista1 = []
    idheader = None
    for lines in DNAfile:
        if lines.startswith('>'):
            if idheader:
                DNA_dict[idheader]=''.join(lista1) # adding the id to the dic
                del lista1[:]   #empty the list for the next iteration

            idheader = lines.strip().replace('>','')    #remove line changes, get rid of the >, save as the id header
        else:
            lista1.append(lines.strip())    #otherwise we're in the sequence line, append to the list and add into the dict as value
            #print(lista1)
    DNA_dict[idheader]=''.join(lista1)       # need to add this once more
                                            #as the last id+seq ends with a seq
                                            #instead of the next id
    del lista1[:]

#repeat the same procedure as for mtDNA dict to create a dict for Y chrom
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
    #print(Ychrom_dict)

    def Scores(scoring):
        transition = ['AG', 'TC', 'GA', 'CT'] #the transition scores
        NTscoreTotal=0     #a counter to count identical NTs. Done later on in the parallel iteration
        Seqscore_list=[]  #list for the total scores of NTs when comparing two seqs
        Identical_NTs=0   #for counting the identical ones once we start off the loops

        checked_pair=[]

        for key1, seq1 in scoring.items():
            checked_pair.append(key1) #NEW: Add the id to the list of compared sequences
            for key2, seq2 in scoring.items():
                if key2 not in checked_pair:     #so we only add the key1-key2 comb if its not already in the set, avoiding repeats
                    # Everything after the if is indented
                    Identical_NTs=0

                    for NTa,NTb in zip(seq1, seq2): #parallel iteration-in one loop calculates the overall score for the alignment

                        if NTa==NTb:
                            if '-' in NTa and '-' in NTb:    #if empty positions in both NTs, then score is 0
                                score=0
                                #NTscoreTotal.append(score)
                                NTscoreTotal +=score

                            else:      #else the nucleotides pair up with content in them
                                score=1   #score is 1 when the NTs match
                                NTscoreTotal +=score   #add this to the total score
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
                        perc_identity=Identical_NTs/len(seq1)*100   #calculate perc.identity. Len
                        perc_identity=round(perc_identity,2)

                        Header_Seq = str(key1) + '\t' + str(key2) + '\t' + str(perc_identity) + '%\t' + str(NTscoreTotal) + '\n'   #a string, works as the header id, i.e. the two person's seqs that are compared
                    Seqscore_list.append(Header_Seq)
                    NTscoreTotal=0
        #print(Seqscore_list)
        return Seqscore_list



    mtDNA_results=Scores(DNA_dict)
    Ychrom_results=Scores(Ychrom_dict)
    #print(mtDNA_results)
    #print(Ychrom_results)   #encoding had to be performed to be able to run the outputs in part 2
with open('output_mtDNA.txt','w', encoding='utf-8') as mtDNAoutput, open('output_ychr.txt','w', encoding='utf-8') as Ychroutput:
     Header="sample1 \t sample2 \t Perc_identity \t Score \n"
     mtDNAoutput.write(Header)
     Ychroutput.write(Header)
     for line in mtDNA_results:
         #print(line)
         mtDNAoutput.write(line)
     for line in Ychrom_results:
         #print(line)
         Ychroutput.write(line)

