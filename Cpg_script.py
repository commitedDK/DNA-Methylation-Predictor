#!/usr/bin/python3

import sys
import re

def main():# this function call all the functions used
    verifyInput()
    input_file = sys.argv[1]
    pattern = 'cg'
    ids, sequen = making_seq_nd_ids_lists(input_file)
    methylatedC,unmethylatedC,window_name = sliding_window(sequen,pattern)
    converted_seq = bisulfite_conversion(unmethylatedC,sequen)
    digests = restriction_analysis(converted_seq)
    percent_methylation(digests)


def verifyInput(): # for user to run the script correctly
    if len(sys.argv) != 2:
        print("Usage: {} {}".format(sys.argv[0], 'upstream1000.fa'))


def making_seq_nd_ids_lists(person_file) :
    with open(person_file,'r') as file_fna:
        seq_list = []
        ids = []
        seq = ''
        for line in file_fna:
            line = line.rstrip()
            if line.startswith('>'):
                if seq != '':
                    seq_list.append(seq)
                    seq = ''
                ids.append(line)
            else:
                seq += line
        seq_list.append(seq)
    return ids,seq_list

def sliding_window(seq,pat) :
    my_seq = seq[0]
    window_len = 200
    methylatedC, unmethylatedC, window_name = [],[],[]
    
    for index,value in enumerate(range(0, len(my_seq), 1)):
        if len(my_seq[value:value + window_len]) == window_len:
            seq_low = my_seq[value:value + window_len].lower()
            gc_count = seq_low.count('g') + seq_low.count('c')
            ta_count = seq_low.count('a') + seq_low.count('t')
            gc_expected = seq_low.count('g') * seq_low.count('c')/window_len
            gc_content = gc_count / (gc_count + ta_count)
            gc_observed = seq_low.count(pat)
            ratio = gc_observed/gc_expected
            ratio = round(ratio,2)
            
            if (gc_content) > 0.5 and ratio > 0.6:
                upperlimit=str(index+window_len)
                windowname= "CpG island found in:"+str(index)+':'+upperlimit +'  Ratio:'+str(ratio)+'   GC Content'+str(gc_content)
                methC_pos = my_seq.find('c',index,index+window_len)
                window_name.append(windowname)
                methylatedC.append(methC_pos)
            else :
                unmethC_pos = my_seq.find('c',index,index+window_len)
                unmethylatedC.append(unmethC_pos)             
    methylatedC = list(set(methylatedC))
    methylatedC = sorted(methylatedC)
    unmethylatedC = list(set(unmethylatedC))
    unmethylatedC = sorted(unmethylatedC)
    #print(methylatedC,unmethylatedC)
    print(window_name)
                   
    return methylatedC, unmethylatedC, window_name

def bisulfite_conversion(unmethylatedC,seq) :# converts the cytosine not present in cpg islands only to thymine 
    my_seq = seq[0]
    for i in unmethylatedC :
        my_seq = my_seq[0:i] + 't' + my_seq[i+1:]
        
    return my_seq

def restriction_analysis(converted_seq):
    con_seq_list,digests = [],[]
    con_seq_list.append(0)# it adds last position of dna into dna_list
    for res_seq in re.finditer('ccgg', converted_seq):#find all restriction sites
        con_seq_list.append(res_seq.start())#find first position of every restriction site
    con_seq_list.append(len(converted_seq) - 1)# it adds last position of dna into con_seq_list
    print(con_seq_list)
   
    for idx,value in enumerate(con_seq_list):# it iterated through con_seq_list and for every element that it hasfound it returns its index in idx and element 
        nextelem = idx + 1# it finds the next element for the found element
        if nextelem < len(con_seq_list):
            if value == con_seq_list[0]:# then it checks if it is the first element of the con_seq_list
                digests.append(converted_seq[con_seq_list[idx]:con_seq_list[nextelem]+2])#in this particular case only it doesnt need to find the seond position of the restriction site, instead it shouldfind the starting point of converted_seq, and adds 0 position to 2nd position of first time it found restriction site
                con_seq_list.pop(0)#and then remove the 0th position of con_seq_list

            str_idx = con_seq_list[idx] + 2 # all other cases add one to the starting and ending point in converted_seq
            digests.append(converted_seq[str_idx:con_seq_list[nextelem] + 2]) ## add the seq to digests

        else :
            pass
    print('Digested segments found:',digests)
    return digests

def percent_methylation(digests):
    no_digested = len(digests)
    no_undigested = 1
    percent_methylation = 100 * (no_digested / (no_digested + no_undigested))
    
    print('Percent methylation:',percent_methylation,'%')

        

main()
