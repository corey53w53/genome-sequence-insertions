#This program reads a file with test sequences and scans the dictionary json file for matches, and a dictionary containing the possible insertions both within and between sequences is exported to a provided json file.

#Command line processing
import argparse
import valet
import json

parser = argparse.ArgumentParser(description="Finds k-mers of a given length in a fasta file, prints out a dictionary with the first 5 k-mers in the DNA sequence")
parser.add_argument('--fastq_file', '-f', required=True, metavar='fastq_file',
                    help="Input the fastq file")  #fastq file with the code samples
parser.add_argument('--json_file', '-j', required=True, metavar='json_file',
                    help="The json dictionary created by kmer_dict.py")  #JSON file with dictionary
parser.add_argument('--kmer_length', '-l', required=True, metavar='kmer_length',
                    help="Length of the k-mers, should be the same as kmer_dict.py")  #K-mer length
parser.add_argument('--out_file', '-o', required=True, metavar='out_file',
                    help="Output json file")  #Output file

args = parser.parse_args()

def reverse_complement(testcase):
    reverse=""
    for letter in testcase[::-1]:
        if letter=='A': reverse+='T'
        elif letter=='C': reverse+='G'
        elif letter=='G': reverse+='C'
        elif letter=='T': reverse+='A'
    return reverse
    """
    NAME: reverse_complement()

    PURPOSE:
        To take a DNA sequence and return the reverse complement to help account for DNA being read in opposite directions when using the lexographically correct sequence

    :param testcase: The DNA sequence
    :type testcase: string
    :return: The reverse complement of the DNA sequence
    :rtype: string
    """

def to_number(testcase):
    key=""
    for letter in testcase:
        if letter=='A': key+='00'
        elif letter=='C': key+='01'
        elif letter=='G': key+='10'
        elif letter=='T': key+='11'
    return int(key,base=2)
    """
    NAME: to_number()

    PURPOSE:
        To convert a DNA sequence into an integer for efficiency purposes when storing and retrieving information in the dictionary.

    :param testcase: The DNA sequence
    :type testcase: string
    :return: The DNA sequence in integer form
    :rtype: integer
    """

#The dictionary is read from the json file and set to the variable dictionary.
with open(args.json_file) as json_file:
    dictionary = json.load(json_file)
f = open(args.fastq_file, "r")
#The first line is read as the important information starts on the 2nd line.
f.readline()
#The length of a kmer is converted to an integer from the user and set to a variable.
kmer_length=int(args.kmer_length)
#The extrema list is created so that 
extrema_list=[]
#The intersequence dict stores the values of extremas that are from different sequences of the fasta file.
intersequence_dict={}
selective_intersequence_dict={}
#The 2nd line is read and set to variable line.
line=f.readline().strip()
#The max match list is the highest number that can be achieved if every k-mer in a single fastq line matched a value from the dictionary.
max_match_list_len=len(line)-kmer_length
extrema_dict={}
#This variable tracks the splits that span multiple sequences on the fasta file.
intersequence_list=[]
#This loop runs to the end of the file
while line!="":
    #The match list is reset for each reference sequence.
    match_list=[]
    no_match_counter=0
    for number in range(len(line)-kmer_length+1):
        #The k_mer variable is created and set to start at number with length kmer_length
        k_mer=line[number:number+kmer_length]
        #The reverse complement of k_mer is taken if it is lexographically correct.
        if k_mer>reverse_complement(k_mer):
            k_mer=reverse_complement(k_mer)
        #The k_mer is converted to a number
        k_mer=to_number(k_mer)
        if str(k_mer) in dictionary:
            for list in dictionary[str(k_mer)]:
                #The sequence number is set as sequence_counter and each value of the list is appended as a tuple to match_list
                sequence_counter=list[0]
                for number in list[1:]:
                    match_list.append((sequence_counter,number))
        else:
            #If there are no matches, then the no_match_counter is incremented.
            no_match_counter+=1
            if no_match_counter>(max_match_list_len*0.25):
                break
    #This condition only runs if there are enough matches in the match_list.
    if len(match_list)>(max_match_list_len*0.75):
        #The boolean same_sequence tests if all the tuples in match_list are from the same sequence, or if they are from different sequences.
        same_sequence=True
        sequence=match_list[0][0]
        for tuple in match_list:
            if tuple[0]!=sequence:
                same_sequence=False
        #If the tuples are from the same sequence, then the match_list is set to just be the 2nd value in the tuple (at position 1).
        if same_sequence:
            for num in range(len(match_list)):
                match_list[num]=match_list[num][1]
            match_list.sort()
            difference=match_list[-1]-match_list[0]
            predicted_split_length=1000
            #The difference is the change from the smallest value to the largest value, the filter checks to see if it is greater than 500 but less than 2000.
            if difference<(predicted_split_length*2) and difference>(predicted_split_length/2):
                #A list containing the extrema is made, and any value who is not consecutive on both sides of the value will be added.
                counter=0
                #This loop runs until the end of the match_list. A while loop is used so that the difference between a position and the next position can be found
                while counter<len(match_list)-2:
                    if match_list[counter+1]-match_list[counter]!=1:
                        if extrema_dict.get(sequence)==None:
                            extrema_dict[sequence]=[match_list[counter]+kmer_length,match_list[counter+1]]
                        else:
                            extrema_dict[sequence].append(match_list[counter]+kmer_length)
                            extrema_dict[sequence].append(match_list[counter+1])
                    counter+=1
        #The following code runs if the match_list contains tuples from different sequences from the fasta file.
        else:
            match_list.sort()
            for value in range(len(match_list)-1):
                #Each value is tested to find the tuples where the list switches from one sequence to another.
                if match_list[value][0]!=match_list[value+1][0]:
                    #The intersequence tuple is set to a tuple of tuples, each having the sequence number as the first value and the position as the second value.
                    intersequence_tuple=((match_list[value],match_list[value+1]))
                    #The intersequence tuple is added to the dictionary if it is not already there.
                    if intersequence_dict.get(intersequence_tuple)==None:
                        intersequence_dict[intersequence_tuple]=1
                    else:
                        intersequence_dict[intersequence_tuple]+=1
                        intersequence_list=sorted(intersequence_dict.items(), key=lambda kv: kv[1])
                        intersequence_list.reverse()
                    #The selective dictionary is sorted and reversed so that the most common values are shown at the beginning of the list.
    #3 lines are read so that the next iteration will start with a relevant line.
    f.readline()
    f.readline()
    f.readline()
    line=f.readline().strip()
#The extrema_list is sorted so it can be used by the poisswin function.

#The final dictionary contianing all the important variables is made.
final_dict={}
if intersequence_list:
    final_dict['intersequence_list'] = intersequence_list
for sequence in extrema_dict:
    extrema_list=extrema_dict[sequence]
    #extrema_list is sorted so the poisswin function can be used.
    extrema_list.sort()
    final_dict["extrema_list"+str(sequence)]=extrema_list
    poisswin_list=valet.poisswin(extrema_list,extrema_list[-1])
    final_dict["poisswin_list"+str(sequence)]=poisswin_list
    for poisswin_dict in poisswin_list:
        #The best end must be less than the length of the extrema list to avoid an index error.
        if poisswin_dict['be']<len(extrema_list):
            #The best start and best end are found by plugging the values of bs and be from the dictionary back into the extrema list.
            best_start=extrema_list[poisswin_dict['bs']]
            best_end=extrema_list[poisswin_dict['be']]
        elif poisswin_dict['bs']>0:
            #The best start and best end are found by plugging the values of bs and be from the dictionary back into the extrema list.
            best_start=extrema_list[poisswin_dict['bs']-1]
            best_end=extrema_list[poisswin_dict['be']-1]
        if final_dict.get("best_split"+str(sequence))==None:
            final_dict["best_split"+str(sequence)]=[(best_start,best_end)]
        else:
            final_dict["best_split"+str(sequence)].append((best_start,best_end))

#The final dictionary is printed and also dumped into the given json out_file.
print(final_dict)
with open(args.out_file, "w") as outfile:
    json.dump(final_dict, outfile)