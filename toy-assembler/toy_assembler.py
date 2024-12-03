import argparse
from itertools import groupby

parser = argparse.ArgumentParser(description='toy assembler')
parser.add_argument('-i', type = str, required = True)
args = parser.parse_args()
# # # # # # # # # # # # 
class IdSequence:
    def __init__(self, id, sequence):
        self.id = id
        self.sequence = sequence
# # # # # # # # # # # 
def fasta_iter(filename):
    fin = open(filename, 'rb')
    faiter = (x[1] for x in groupby(fin, lambda line: str(line, 'utf-8')[0] == ">"))
    for header in faiter:
        headerstr = str(header.__next__(), 'utf-8')
        long_name = headerstr.strip().replace('>', '')
        name = long_name.split()[0]
        seq = "".join(str(s, 'utf-8').strip() for s in faiter.__next__())
        yield name, seq

sequences = []
for i in fasta_iter(args.i):
    n = IdSequence(i[0], i[1])
    sequences.append(n.sequence)
sequences = sorted(sequences, key = len, reverse = True)

# # # # # # # # # # # # 
def score(sequence1, sequence2, offset):
    start_of_overlap = max(0 - offset, 0)
    end_of_overlap = min([len(sequence2) - offset, len(sequence2), len(sequence1) - offset])
    total_score = 0
    for position in range(start_of_overlap, end_of_overlap):
        if sequence2[position] == sequence1[position + offset]:
            total_score += 1
    return total_score
    
def find_best_offset(sequence1, sequence2):
    lowest_offset = 1 - len(sequence2)
    highest_offset = len(sequence1)
    all_offsets = []
    for offset in range(lowest_offset, highest_offset):
        # add the 4-tuple for this offset
        all_offsets.append((score(sequence1, sequence2, offset), offset, sequence2, sequence1))
    return max(all_offsets)
    
def find_best_match(sequence1, others):
    all_matches = []
    for sequence2 in others:
        if sequence2 != sequence1:
            all_matches.append(find_best_offset(sequence1, sequence2))
    return max(all_matches)

def consensus(score, offset, sequence1, sequence2):
    sequence2_left_overhang = sequence2[0:max(0, offset)]
    #sequence2_right_overhang = sequence2[len(sequence1) + offset:len(sequence2)]
    if len(sequence1) + offset > len(sequence2):
        sequence2_right_overhang = ""
    else:
        sequence2_right_overhang = sequence2[len(sequence1) + offset:]
    return sequence2_left_overhang + sequence1 + sequence2_right_overhang
    
def assemble(sequence, others):
    # remember, best_matching_other is a 4-tuple
    best_matching_other = find_best_match(sequence, others)
    # the * expands the elements of the tuple so we can use them as arguments to consensus()
    consensus_sequence = consensus(*best_matching_other)
    if len(others) == 1:
        return consensus_sequence
    else:
        # get the second element of the best_matching_other tuple, which is the sequence
        best_matching_sequence = best_matching_other[2]
        others.remove(best_matching_sequence)
        return assemble(consensus_sequence, others)

def assemble_helper(dnas):
    return assemble(dnas[0], dnas[1:])
# # # # # # # # # # # # 
f = open("assembly.txt", "w")
f.write(assemble_helper(sequences))
f.close()
