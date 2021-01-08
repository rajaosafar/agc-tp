#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
from operator import itemgetter
import itertools
import nwalign3 as nw


__author__ = "Sarah Rajaosafara"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Sarah Rajaosafara"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Sarah Rajaosafara"
__email__ = "rajaosafar@eisti.eu"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default=400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default=10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default=100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default=8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    with gzip.open(amplicon_file, "rt") as myfile:
        sequence = ""
        for line in myfile:
            # Update the sequence if it is written on multiple lines
            if not line.startswith(">"):
                sequence += line.replace('\n', '')
            # Yield the sequence
            else:
                if len(sequence) >= minseqlen:
                    yield sequence
                sequence = ""
        # Yield the last sequence
        if len(sequence) >= minseqlen:
            yield sequence


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    # Create a dictionnary with the occurence for each sequence
    seq_dict = dict()
    for seq in read_fasta(amplicon_file, minseqlen):
        if seq in seq_dict.keys():
            seq_dict[seq] += 1
        else:
            seq_dict[seq] = 1

    # Create a list with all the sequences having an occurence >= mincount
    seq_list = []
    for seq, occ in seq_dict.items():
        if occ >= mincount:
            seq_list.append((seq, occ))

    # Sort the list
    seq_list.sort(key=itemgetter(1), reverse=True)

    # Yield
    for item in seq_list:
        yield item


def get_chunks(sequence, chunk_size):
    segments = []
    for k in range(4):
        # Create 4 segments of size chunk_size
        segments.append(sequence[k*chunk_size:(k+1)*chunk_size])
    return segments


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    return list(set(lst1) & set(lst2))


def cut_kmer(sequence, kmer_size):
    for i, _ in enumerate(sequence):
        # Only take kmer of size kmer_size
        if i <= len(sequence)-kmer_size:
            yield sequence[i:i+kmer_size]


def get_identity(alignment_list):
    nb_similar = 0
    for i in range(len(alignment_list[0])):
        # Count the number of similar letters
        if alignment_list[0][i] == alignment_list[1][i]:
            nb_similar += 1
    # Divide by the length
    return nb_similar / len(alignment_list[0]) * 100


def search_mates(kmer_dict, sequence, kmer_size):
    list_ids = []
    # Get all the ids
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer in kmer_dict.keys():
            for i in kmer_dict[kmer]:
                list_ids.append(i)

    # Keep the 8 most common ids
    return [id for id, _ in Counter(list_ids).most_common(8)]


def detect_chimera(perc_identity_matrix):
    std_list = []
    parent0 = 0
    parent1 = 0
    for line in perc_identity_matrix:
        # Create a list of all the stds to compute the mean
        std_list.append(statistics.stdev(line))
        # If this segment is similar to the parent 0,
        # increment the number of segments similar to the parent 0
        if line[0] > line[1]:
            parent0 += 1
        # If this segment is similar to the parent 1,
        # increment the number of segments similar to the parent 1
        else:
            parent1 += 1

    # Check if this is a chimera
    return statistics.mean(std_list) > 5 and parent0 > 0 and parent1 > 0


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    # Get all the sequences
    generator = dereplication_fulllength(amplicon_file, minseqlen, mincount)

    # Add the 2 first sequences to the list of references (non chimeric sequences)
    references = []
    i = 0
    for line in generator:
        references.append(line)
        if i == 2:
            break

    # Yield the 2 first sequences
    for ref in references:
        yield ref

    # For all the other sequences
    for candidate, occurence in generator:

        is_chimera = False

        # Get the segments of all the references
        segments_reference_list = []
        for reference, _ in references:
            segments_reference_list.append(get_chunks(reference, chunk_size))

        # Get the kmers of the segments of the references
        kmers_segments_reference_list = []
        for segments_reference in segments_reference_list:
            kmers_segments_reference = []
            for segment in segments_reference:
                kmers_reference = list(cut_kmer(segment, kmer_size))
                kmers_segments_reference.append(kmers_reference)
            kmers_segments_reference_list.append(kmers_segments_reference)

        # 1. Divide each candidate in 4 segments of size L= chunk_size
        segments_candidate = get_chunks(candidate, chunk_size)

        # 2. For each segment, identify 8 sequences with similar kmers
        # Get the kmers of the segments of the candidate
        kmers_segments_candidate = []
        for segment in segments_candidate:
            kmers_candidate = list(cut_kmer(segment, kmer_size))
            kmers_segments_candidate.append(kmers_candidate)

        list_mates = []

        # For each list of kmers for each segment
        for k, kmer_list in enumerate(kmers_segments_candidate):
            kmer_dict = dict()
            # For each kmer of segment k
            for kmer in kmer_list:
                # For each segment n°k of each reference
                for i, seg in enumerate([a[k] for a in kmers_segments_reference_list]):
                    if kmer in seg:
                        if kmer in kmer_dict.keys():
                            if i not in kmer_dict[kmer]:
                                kmer_dict[kmer].append(i)
                        else:
                            kmer_dict[kmer] = [i]
            # For each segment, identify 8 sequences with similar kmers
            list_mates.append(set(search_mates(kmer_dict, candidate, kmer_size)))

        # 3. Find two parents
        parents_ids = set.intersection(*list_mates)

        # 4. Compute the similarities
        if len(parents_ids) > 1:
            # If we have more than 2 parents, test with all the combinations of parents
            for parent1, parent2 in itertools.combinations(parents_ids, 2):

                # Create the matrix
                perc_identity_matrix = [[] for k in range(4)]
                for k in range(4):
                    for id_parent in [parent1, parent2]:
                        segments_parent = get_chunks(references[id_parent][0], chunk_size)

                        alignement = nw.global_align(segments_candidate[k], segments_parent[k],
                                                     gap_open=-1, gap_extend=-1,
                                                     matrix=os.path.abspath(os.path.join(os.path.dirname(__file__), "MATCH")))

                        perc_identity_matrix[k].append(get_identity(alignement))

                # Detect a chimera
                is_chimera = detect_chimera(perc_identity_matrix)
                if is_chimera:
                    break

        # If this is not a chimera, add it to the references and yield
        if not is_chimera:
            references.append((candidate, occurence))
            yield candidate, occurence


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    # Get all the non chimeric sequences
    generator = chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)

    # Add the first sequence to the list of OTUs
    otu_final = []
    otu_final.append(next(generator))

    # For each other sequence
    for sequence, occ in generator:
        is_otu = True
        # Compare it to each sequence in the final list
        # (because the other sequences have a lower occurence)
        for sequence2, occ2 in otu_final:
            alignement = nw.global_align(sequence, sequence2,
                                         gap_open=-1, gap_extend=-1,
                                         matrix=os.path.abspath(os.path.join(os.path.dirname(__file__), "MATCH")))
            if occ2 > occ and get_identity(alignement) > 97:
                is_otu = False
                break
        if is_otu:
            otu_final.append((sequence, occ))

    return otu_final


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_OTU(OTU_list, output_file):
    # Write each OTU to the file
    with open(output_file, "wt") as myfile:
        for i, (otu, occ) in enumerate(OTU_list):
            myfile.write(">OTU_"+str(i+1)+ " occurrence:"+str(occ)+"\n")
            myfile.write(fill(otu, width=80)+"\n")


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    amplicon_file = args.amplicon_file
    minseqlen = args.minseqlen
    mincount = args.mincount
    chunk_size = args.chunk_size
    kmer_size = args.kmer_size
    output_file = args.output_file

    OTU_list = abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)
    write_OTU(OTU_list, output_file)


if __name__ == '__main__':
    main()
