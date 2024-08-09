#!/usr/bin/env python

import sys
import pysam
from Levenshtein import distance as levenshtein_distance

class gene:
    def reverse_complement(seq):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement.get(base, base) for base in reversed(seq))
    def encode_sequence(sequence):
        """Encodes a DNA sequence using 2-bit encoding."""
        encoding = {'A': 0b00, 'C': 0b11, 'G': 0b01, 'T': 0b10}
        encoded = 0
        for nucleotide in sequence:
            encoded = (encoded << 2) | encoding[nucleotide]
        return encoded


def get_duo_soft(read, seq):
    softclip_seqL = []
    cigar = read.cigartuples
    if cigar:
        # Check for left soft-clipping
        if cigar[0][0] == 4:
            left_soft_clip_length = cigar[0][1]
            left_soft_clip_seq = seq[:left_soft_clip_length]
            softclip_seqL.append(left_soft_clip_seq)  
        
        # Check for right soft-clipping
        if cigar[-1][0] == 4:
            right_soft_clip_length = cigar[-1][1]
            right_soft_clip_seq = seq[-right_soft_clip_length:]
            softclip_seqL.append(right_soft_clip_seq)  
    
    return softclip_seqL

def comp(jackpot):
    if jackpot == "3":
        return "5"
    elif jackpot == "5":
        return "3"
    elif jackpot == "r3":
        return "r5"
    elif jackpot == "r5":
        return "r3"
    else:
        return None  # Return None or raise an error if the input is invalid


def main():
    inputbam = sys.argv[1]
    directory = sys.argv[2]
    target = sys.argv[3]
    # target = "CCTGTACTTCGTTCAGTTACGTATTGCTAATGATACGGCGACCACCGAGATCTACACTATAGCCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
    rtarget = reverse_complement(target)
    max_distance = 5
    soft_dict = {}
    whole_num = 0
    match_num = 0 
    bam = pysam.AlignmentFile(inputbam, "rb")
    with open(f"{directory}/c_soft_{inputbam}.txt", 'w', encoding="utf-8") as correspondent:
        with open(f"{directory}/match_soft_{inputbam}.txt", 'w', encoding='utf-8') as output:
            # output.write(f"Primer : {target}\nReverse Primer : {rtarget}\n>\n")
            for read in bam: ###--
                whole_num += 1
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                
                curr_seq = read.query_sequence
                try:
                    duo_soft = get_duo_soft(read,curr_seq)
                    soft_dict['5'] = duo_soft[0]
                    soft_dict['3'] = duo_soft[1]
                    soft_dict['r5'] = reverse_complement(duo_soft[0])
                    soft_dict['r3'] = reverse_complement(duo_soft[1])


                    jackpot = ''
                    for key in soft_dict:
                        if len(soft_dict[key]) > len(target):
                            for scan_num in range(len(soft_dict[key]) - len(target)):
                                subseq = soft_dict[key][scan_num:scan_num + len(target)]
                                if levenshtein_distance(subseq, target) <= max_distance:
                                    jackpot = key
                                    break
                            if jackpot:
                                output.write(f">{jackpot}\' | {read.query_name} | f:{read.flag}\n{soft_dict[jackpot]}\n")
                                correspondent.write(f">{comp(jackpot)}\' | {read.query_name} | f:{read.flag}\n{soft_dict[comp(jackpot)]}\n")
                                #  {len(soft_dict[jackpot])}
                                # output.write(f"{subseq}\n>\n{curr_seq} \n{jackpot}\' | {read.cigarstring} | {read.query_name} | flag : {read.flag}\n@\n")
                                match_num += 1
                                jackpot =''
                                break

                        else:
                            pass        
                except IndexError:
                    pass                

    bam.close()
    per = round(match_num/whole_num * 100 , 2)
    with open(f"{directory}/stat_{inputbam}.txt","w", encoding='utf-8') as stat:
        stat.write(f"match percentage : {per}")


if __name__ == "__main__":
    main()
            
        

