import sys

with open(sys.argv[1], "r") as ref_dict_file:
    sequence_tuple_list = []
    longest_sequence = 0
    for line in ref_dict_file:
        if line.startswith("@SQ"):
            line_split = line.split("\t")
            # (Sequence_Name, Sequence_Length)
            sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
    longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    print(sequence_tuple_list[0])
    print(longest_sequence)
# We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
# the last element after a :, so we add this as a sacrificial element.
hg38_protection_tag = ":1+"
# initialize the tsv string with the first sequence
tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
temp_size = sequence_tuple_list[0][1]
max_line = 0
for sequence_tuple in sequence_tuple_list[1:]:
    if temp_size + sequence_tuple[1] <= longest_sequence and max_line <= 1600:
        temp_size += sequence_tuple[1]
        tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        max_line = max_line + 1
    else:
        tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
        temp_size = sequence_tuple[1]
        max_line = 0
# add the unmapped sequences as a separate line to ensure that they are recalibrated as well
sequence_groups = tsv_string.split("\n")
print(sequence_groups[0])
for i in range(0, len(sequence_groups)):
    with open("sequence_grouping_{0}.txt".format(i), "w") as tsv_file:
        tsv_file.write(sequence_groups[i])
        tsv_file.close()

tsv_string += '\n' + "unmapped"
sequence_groups_unmapped = tsv_string.split("\n")
for i in range(0, len(sequence_groups_unmapped)):
    with open("sequence_grouping_with_unmapped_{0}.txt".format(i), "w") as tsv_file_with_unmapped:
        tsv_file_with_unmapped.write(sequence_groups_unmapped[i])
        tsv_file_with_unmapped.close()
