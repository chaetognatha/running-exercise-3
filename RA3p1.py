
def gen_data_parser(fh):
    with open(fh) as f:
        mt_dict = {}
        y_dict = {}
        my_count = 0
        for line in f:
            if my_count == 0:
                name = line.strip()
                my_count+=1
                continue
            if line.startswith("mtDNA"):
                seq = next(f)
                edit = seq.strip()
                mt_dict[name] = edit
                seq = ''
                edit = ''
                continue
            if line.startswith("Y chromosome"):
                seq = next(f)
                edit = seq.strip()
                y_dict[name] = edit
                seq = ''
                edit = ''
                continue
            if "hemophilia" in line:
                continue
            if line.strip() != '':
                name = line.strip().strip(">")
                continue
        return mt_dict, y_dict


def scoring(seq1, seq2, header1, header2): #to check the alignment score
    transition = ['AG', 'TC', 'GA', 'CT'] #the transition scores
    score_list = [] #a list to add all the scores when they are calculated
    identical_nt_list = []
    for a, b in zip(seq1, seq2): #this for loop goes through the nucleotides in seq1 and seq2 (a and b respectively) and compares them
        a = a.upper() #in case the sequence is lower case
        b = b.upper() #in case the sequence is upper case
        if a == b: #if the nucleotide in a is the same as in b
            if '-' in a and '-' in b: #if both are - the score is 0 and is added to the list score_list
                score = 0
                score_list.append(score)
            if '?' in a and '?' in b:
                score = -1
                score_list.append(score)
            else: #if they are not - they are nucleotides and the score is 1, which is added to the score_list
                score = 1
                score_list.append(score)
                identical_nt_list.append(score)
        else: #if the nucleotides are not the same
            if '-' in a or '-' in b: #if one sequence has a gap when aligned to the other one the score is -1 and is added to score_list
                score = -1
                score_list.append(score)
            if '?' in a or '?' in b:
                score = -1
                score_list.append(score)
            else:
                nt_sum = a + b
                if nt_sum in transition: #if they are in transition the score is -1 and it is appended to score_list
                    score = -1
                    score_list.append(score)
                else: #if it is not a transition it means that it is a transversion, which has a score of -2 that is added to score_list
                    score = -2
                    score_list.append(score)
    align_score = str(sum(score_list)) #then i do a sum to sum all the values in the score_list and assign it to a variable score_sum
    id_score = str(sum(identical_nt_list)/len(seq1)*100)
    return id_score, align_score #id and alignment must be returned


def dict_score_loop(dict):
    my_list = []
    for k1, v1 in dict.items():
        for k2, v2 in dict.items():
            id_score, align_score = scoring(v1, v2, k1, k2)
            row = [k1, k2, id_score, align_score]
            my_list.append(list(row))
    return my_list
mt_dict, y_dict = gen_data_parser("GeneticData.txt")
mt_score_out = dict_score_loop(mt_dict)
y_score_out = dict_score_loop(y_dict)
#print(mt_score_out)
#print(y_score_out)

def un_wrapper(my_list):
    check_list = []
    new_list = []
    for row in my_list:
        print(row[0])
        if row[0] != row[1]:
            if row[1] not in check_list:
                new_list.append(row)
        check_list.append(row[0])


    return new_list

mt_unwrapped = un_wrapper(mt_score_out)
y_unwrapped = un_wrapper(y_score_out)

def printer(my_list_of_lists, output_file):
    header = "SampleA\tSampleB\tIdentityScore\tScore"
    with open(output_file, 'w') as f:
        print(f"{header}", file=f)
        for row in my_list_of_lists:
            print(f"{row[0]}\t{row[1]}\t{row[2]}%\t{row[3]}\t", file=f)
printer(mt_unwrapped, "mt_test.txt")
printer(y_unwrapped, "y_test.txt")