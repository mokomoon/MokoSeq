from Bio import SeqIO

def parse_fasta(file_path):
    """Return the first sequence from a FASTA file."""
    for record in SeqIO.parse(file_path, "fasta"):
        return record.seq
    return ""

