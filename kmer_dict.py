#This program reads multiple strands from a FASTA file and outputs a dictionary into a json file with a k-mer in number form as the key and its positions as the values.

#Command line processing
import argparse

parser = argparse.ArgumentParser(description="Finds k-mers of a given length in a fasta file, prints out a dictionary with the first 5 k-mers in the DNA strand")
parser.add_argument('--infile', '-i', required=True, metavar='fasta_file',
                    help="Input the fasta file that this program should find k-mers")  # Input file
parser.add_argument('--kmer_length', '-l', required=True, metavar='kmer_length',
                    help="Length of the k-mers")  #K-mer length
parser.add_argument('--outfile', '-o', required=True, metavar='output_file',
                    help="Output file for the json dictionary")  # Output file

args = parser.parse_args()

import json

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

def add(list,addition):
    final=()
    for counter in range(len(list)-1):
        final+=(list[counter],)
    last=list[-1]+(addition,)
    final+=(last,)
    return final
    """
    NAME: add()

    PURPOSE:
        To append an integer to the end of the last tuple in a list

    :param list: The original list of tuples
    :type list: list of tuples
    :param addition: The number being added
    :type addition: integer
    :return: The list of tuples with the number appended
    :rtype: list of tuples
    """
    
#The infile is opened, and an error message is printed if there is an error.
try:
    f = open(args.infile, "r") # open input file
except OSError as err:
    print("Could not open input: {e}".format(e=err))
    exit(1)

dictionary={}
#The length of the kmer can be changed by changing the kmer_length variable.
kmer_length=int(args.kmer_length)
position=0
#The first line is read and set to testcase.
testcase=f.readline()
#The second line is read and set to next_testcase, which is necessary for reading k-mers that span more than 1 line.
next_testcase=f.readline()
#The line counter is made to determine the position of the k-mer with respect to the entire file.
line=0
#The strand number is counted so the output can include strand numbers.
strandcounter=0
#This loop runs until the end of the file.
while testcase!="":
    #To get rid of the first line, which does not contain relevant information, the next_testcase is set to testcase.
    testcase=next_testcase.strip()
    #next_testcase is set to the next line.
    next_testcase=f.readline()
    #If the length of next_testcase is greater than 0 and it starts with >, this means that a strand has ended. Thus, the strand counter will be incremented and the line will be reset to 0.
    if testcase:
        if testcase[0]==">":
            line=0
            strandcounter+=1
            #The next testcase is read
            testcase=next_testcase.strip()
            next_testcase=f.readline()
    #The actual length of the testcase is necessary to calculate the position of the genome.
    actual_length=len(testcase)
    if next_testcase:
        # Testcase adds on the next (kmer_length-1) characters from the next line so that k-mers that require 2 lines are included.
        for char in next_testcase[0:kmer_length-1].strip():
            #The character "N" is not added
            if char!="N":
                testcase+=char
    #Extended_length is necessary to calculate the range of the for-loop.
    extended_length=len(testcase)
    #This loop repeats so that all the positions are included so that position+kmer_length does not exceed the length of the testcase, which is extended_length.
    for position in range(0,extended_length-kmer_length+1):
        #A string starting at position with a length of kmer_length is stored in new_segment.
        new_segment=testcase[position:kmer_length+position]
        #This if-statement checks to see whether the new segment is already lexigraphically correct, or if a reverse complement would be necessary.
        if new_segment>reverse_complement(new_segment):
            new_segment=reverse_complement(new_segment)
        #If the string new_segment is not empty, it is converted into a number.
        if len(new_segment)!=0:key=to_number(new_segment)
        #The program tries to add the position to the new segment if the last term is not -1 and if the key is already created, if the key is not created, the except statement does so.
        try:
            #If the first number of the last list is not equal to the strand counter, a new list is made with the strand counter being the first value of the list.
            if strandcounter != dictionary[key][-1][0]:
                dictionary[key]+=((strandcounter,),)
                #The position of that testcase is added right after the strand counter
                dictionary[key]=add(dictionary[key],(position+(actual_length*line)))
            else:
                #If the last number of the last list is not equal to -1, the position is added to the last list.
                if -1 != dictionary[key][-1][-1]:
                    dictionary[key]=add(dictionary[key],(position+(actual_length*line)))
                    #If there are greater than 5 occurences of the segment (6 numbers in the list due to the strand counter), a -1 will be added. 
                    if len(dictionary[key][-1])==6:
                        dictionary[key]=add(dictionary[key],-1)
        except KeyError:
            #If there is not already a key in the dictionary and there are no "N"s in the string and the length of the new_segment is equal to the kmer_length, the dictionary adds the new_segment.
            if "N" not in new_segment and len(new_segment)==kmer_length:
                #A 2d array is made with the value strandcounter at the first position
                dictionary[key]=((strandcounter,),)
                #The position of the testcase is added
                dictionary[key]=add(dictionary[key],position+(actual_length*line))
    #The line is incremented to help calculate the position of the k-mer in the genome.
    line+=1
#The dictionary is dumped to the given json file
with open(args.outfile, "w") as outfile:
    json.dump(dictionary, outfile)