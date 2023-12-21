"""The purpose of this program is to analyse DNA sequences and look for the presence of genetic markers. This program
could form the basis for medical research, identifying whether there are connections between different diseases"""

from Bio import Entrez
from Bio import SeqIO
import os


def get_sequence(accession_number):
    """Retrieve DNA sequence from GenBank"""
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record.seq


def find_genetic_markers(sequence, markers):
    found_markers = []
    for marker in markers:
        if marker in sequence:
            found_markers.append(marker)
    return found_markers


email = os.environ.get("MY_EMAIL")

genetic_markers = ["ATCGA", "CAGT", "GGCT", "TACGA"]


# Samples with INS gene
record1 = get_sequence("OM489474")
record2 = get_sequence("KR710185")
record3 = get_sequence("L15440")
record4 = get_sequence("NG_007114")
record5 = get_sequence("MT335692")

ins_gene_records = [record1, record2, record3, record4, record5]
print("Results from DNA samples with INS gene:")
for item in ins_gene_records:
    markers_found = find_genetic_markers(item, genetic_markers)
    print("Genetic Markers Found:", markers_found)


# Samples with APOE gene
record6 = get_sequence("KJ905144")
record7 = get_sequence("NM_001302690")
record8 = get_sequence("NM_000041")
record9 = get_sequence("KU177913")
record10 = get_sequence("DQ286969")

apoe_gene_records = [record6, record7, record8, record9, record10]
print("\nResults from DNA samples with APOE gene:")
for item in apoe_gene_records:
    markers_found = find_genetic_markers(item, genetic_markers)
    print("Genetic Markers Found:", markers_found)


# Samples with MYH7 gene
record11 = get_sequence("MF399812")
record12 = get_sequence("KJ905833")
record13 = get_sequence("NM_001407004")
record14 = get_sequence("NM_000257")
record15 = get_sequence("EF560725")

myh7_gene_records = [record11, record12, record13, record14, record15]
print("\nResults from DNA samples with MYH7 gene:")
for item in myh7_gene_records:
    markers_found = find_genetic_markers(item, genetic_markers)
    print("Genetic Markers Found:", markers_found)
