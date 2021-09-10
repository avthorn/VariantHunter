# By avthorn

# SampleFinder is a companion to VariantHunter

# First argument: A file containing the target gene sequence with or without fasta header

# Second Argument: The sequence_sample_metadata.tsv for the VariantHunter Run with the lowest treashold.

# Third Argument: The sequence_sample_metadata.tsv for the VariantHunter Run with the highest treashold.

# Non essential output lines have # so they can be edited out with | grep -v '#


import os
import sys


target_seq_file = open(str(sys.argv[1]),"r")

target_seq = ""

for line in target_seq_file:
    line = line.strip()
    if not line.startswith(">"):
        target_seq = target_seq + line.upper()

target_seq_file.close()

lower_treashold_file = open(str(sys.argv[2]),"r")
upper_treashold_file = open(str(sys.argv[3]),"r")

print("# Target Sequence file: " + str(sys.argv[1]))

for line in lower_treashold_file:
    line = line.strip()
    line_list = line.split("\t")
    line_seq = line_list[2]
    if line_seq == target_seq:
        seq_name = line_list[0]  
        sample_str = line_list[3]
        sample_l_set = set(sample_str.split(","))
        print("\n## Lower ###############################")
        string = "# Sequence name in lower threshold file " + str(sys.argv[2]) + " is " + seq_name + "."
        print(string)
        string = "\n# Samples are:"
        print("# " + str(sample_l_set))

for line in upper_treashold_file:
    line = line.strip()
    line_list = line.split("\t")
    line_seq = line_list[2]
    if line_seq == target_seq:
        seq_name = line_list[0]
        sample_str = line_list[3]
        sample_u_set = set(sample_str.split(","))
        print("\n## Upper ###############################")
        string = "# Sequence name in upper threshold file " + str(sys.argv[3]) + " is " + seq_name + "."
        print(string)
        string = "\n# Samples are:"
        print("# " + str(sample_u_set))

sample_set = sample_l_set.difference(sample_u_set)

string = "\n\n## Samples to use: \n"
print(string)

for sample in sample_set:
    print(sample)

upper_treashold_file.close()
lower_treashold_file.close()
