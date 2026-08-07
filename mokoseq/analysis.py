def preview_sequence(text):
    return text[:500] + "\n... \n"

from Bio import SeqIO
import csv #creates csv file
import argparse #lets the script accept command-line arguments 
import numpy as np


#Main analysis 
def analyze_sequence(file_path, output_csv):
    with open(output_csv, "w", newline="") as f: #open CSV file for writing
        writer = csv.writer(f)
        writer.writerow(["ID", "Length", "GC_Content", "ORF_Count", "Longest_Peptide_Length", "Longest_Peptide", "Longest_ORF (start, end, length)", "Adv_ORF_Count", "Longest_Adv_ORF",
                         "A_to_I_like", "C_to_U_like", "uORF_Count", "MicroORF_Count", "ncRNA_Type", "Most_Used_Codon", "Most_Used_Count", "Codon_Bias", 
                         "ORF6_Count", "Longest_ORF6", "Motif_Hits", "TATA_Count", "Kozak_Count", "PolyA_Count", "CpG_Count", "ShineDalgarno_Count"]) #header row
        
       # Loop through each sequence in the FASTA file`
        for record in SeqIO.parse(file_path, "fasta"):
            seq = record.seq
            gc = gc_content(seq) #calculate GC content
            orfs = find_orfs(seq) #Find ORFs
            translated_orfs = []
            for start, end, length in orfs: 
                dna = str(seq)[start:end]
                aa = translate_sequence(dna)
                translated_orfs.append((start, end, length, aa))
            longest_peptide = ""
            if translated_orfs:
                longest_peptide = max(translated_orfs, key=lambda x: len(x[3]))[3]
            longest_orf = max(orfs) if orfs else 0 
            advanced_orfs = find_advanced_orfs(seq)
            advanced_orf_count = len(advanced_orfs)
            longest_advanced = max([orf[4] for orf in advanced_orfs]) if advanced_orfs else 0
            editing = detect_rna_editing(seq)
            a_to_i = len(editing["A_to_I_like"])
            c_to_u = len(editing["C_to_U_like"])
            uorf_data = find_uorfs_and_microorfs(seq)
            uorf_count = len([o for o in uorf_data if o[0] == "uORF"])
            micro_count = len([o for o in uorf_data if o[0] == "microORF"])
            ncrna_type = classify_ncrna(seq)
            codons = codon_usage(seq)
            most_used, most_used_count, bias = codon_usage_summary(codons)
            orfs6 = find_orfs_6frame(seq)
            sixframe_count = len(orfs6)
            longest6 = max([o[3] for o in orfs6]) if orfs6 else 0
            motifs = find_motifs(seq)
            motif_hits = sum(len(v) for v in motifs.values())
            tata = len(motifs["TATA_box"])
            kozak = len(motifs["Kozak_consensus"])
            polya = len(motifs["PolyA_signal"])
            cpg = len(motifs["CpG_island"])
            shine = len(motifs["Shine_Dalgarno"])
            writer.writerow([record.id, len(seq), f"{gc:.2f}", len(orfs), len(longest_peptide), longest_peptide, longest_orf, advanced_orf_count, longest_advanced,
                             a_to_i, c_to_u, uorf_count, micro_count, ncrna_type, most_used, most_used_count, bias, sixframe_count, longest6,
                             motif_hits, tata, kozak, polya, cpg, shine]) #Write 1 row per sequence with 2 decimal points 


    print(f"CSV file created: {output_csv}")

#GC content 
def gc_content(seq):
    g = seq.count("G") #count G bases 
    c = seq.count("C") #count C bases 
    return (g + c) / len(seq) * 100 #add the bases and divide by the length of sequence * 100

#ORFs Finder in 3 forward frames 
def find_orfs(seq):
    """
    Returns ORFs in the forward 3 frames.
    Each ORF is returned as (start, end, length)
    """
    seq_str = str(seq)
    results = []
    stop_codons = {"TAA", "TAG", "TGA"}

    for frame in range(3):
        i = frame
        while i < len(seq_str) - 2:
            codon = seq_str[i:i+3]

            if codon == "ATG":
                start = i
                j = i

                # extend ORF
                while j < len(seq_str) - 2:
                    codon_j = seq_str[j:j+3]
                    if codon_j in stop_codons:
                        end = j + 3
                        length = end - start
                        results.append((start, end, length))
                        break
                    j += 3
            i += 3

    return results


#extending the ORF until stop codon
def extend_orf(seq_str, start_index): 
    for i in range(start_index, len(seq_str), 3):
        codon = seq_str[i:i+3]
        if codon in ["TAA", "TAG", "TGA"]: #stop codons
            return (i + 3) - start_index
    return 0 #if no stop codon is found

#advanced ORFs - overlapping, nexted, +1, and -1 frameshifts
def find_advanced_orfs(seq):
    seq_str = str(seq)
    orf_info = []

    stop_codons = {"TAA", "TAG", "TGA"}

    # Precompute codons for speed
    codons = [seq_str[i:i+3] for i in range(0, len(seq_str)-2, 1)]

    # Precompute all start and stop positions
    start_positions = [i for i in range(len(codons)) if codons[i] == "ATG"]
    stop_positions = {i for i in range(len(codons)) if codons[i] in stop_codons}

    # Convert codon index → nucleotide index
    def nt(i): 
        return i * 3

    # Loop through each start codon
    for start_codon in start_positions:
        start_nt = nt(start_codon)

        # NORMAL IN-FRAME ORF
        for stop_codon in range(start_codon + 1, len(codons)):
            if stop_codon in stop_positions:
                end_nt = nt(stop_codon) + 3
                length = end_nt - start_nt
                orf_info.append(("normal", start_codon % 3, start_nt, end_nt, length))
                break

        # +1 FRAMESHIFT
        for stop_codon in stop_positions:
            if stop_codon > start_codon:
                end_nt = nt(stop_codon) + 3
                length = end_nt - start_nt
                orf_info.append(("+1_frameshift", start_codon % 3, start_nt, end_nt, length))
                break

        # -1 FRAMESHIFT
        for stop_codon in stop_positions:
            if stop_codon > start_codon:
                end_nt = nt(stop_codon) + 3
                length = end_nt - start_nt
                orf_info.append(("-1_frameshift", start_codon % 3, start_nt, end_nt, length))
                break

    return orf_info

#RNA editing detection -- inspired by the recent challenge of the Central Dogma
def detect_rna_editing(seq): 
    """
    Detects likely RNA editing events (A->I and C->U) based on sequence patterns. Returns counts and positions."""

    seq_str = str(seq)
    edits = {
        "A_to_I_like": [], #A->G pattern
        "C_to_U_like": [] #C->T pattern
    }

    for i in range(len(seq_str) - 1):
        dinuc = seq_str[i:i+2]

        #A->I editing often occurs in double stranded regions but we detect A->G
        if seq_str[i] == "A" and seq_str[i+1] == "G":
            edits["A_to_I_like"].append(i)

        #C->U editing appears as C->T
        if seq_str[i] == "C" and seq_str[i+1] == "T":
            edits["C_to_U_like"].append(i)
    
    return edits 

#uORFs and Micropeptide Detection 
def find_uorfs_and_microorfs(seq, micro_max_len=300): 
    """
    Detects: 
    - uORFs (upstream ORFS)
    - micro-ORFs (< micro_max_len nucleotides)
    Returns a list of ORFs with type labels."""

    seq_str = str(seq)
    results = []

    stop_codons = {"TAA", "TAG", "TGA"}

    #Scan all 3 forward frames 
    for frame in range(3):
        i = frame 
        while i < len(seq_str) - 2: 
            codon = seq_str[i:i+3]

            if codon == "ATG": 
                start = i
                j = i

                # Extend ORF
                while j < len(seq_str) - 2: 
                    codon_j = seq_str[j:j+3]
                    if codon_j in stop_codons:
                        end = j + 3 
                        length = end - start 

                        # Classify ORF 
                        if start < 100 : # arbitrary UTR threshold 
                            results.append(("uORF", frame, start, end, length))
                        if length <= micro_max_len: 
                            results.append(("microORF", frame, start, end, length))

                        break
                    j += 3 
            i += 3 

    return results 

# Classify non-coding RNA such as tRNA, rRNA, snRNA, miRNA, lncRNA
def classify_ncrna(seq): 
    """
    Classifies non-coding RNA types using simple biological rules.
    Returns one of: 
    - tRNA
    - rRNA
    - miRNA_precursor
    - lncRNA
    - snRNA
    - unknown
    """

    seq_str = str(seq)
    length = len(seq_str)
    gc = (seq_str.count("G") + seq_str.count("C")) / length * 100

    # tRNA: 70-90 nt, high GC, conserved loops
    if 65 <= length <= 95 and gc > 45: 
        if "TTC" in seq_str or "GGG" in seq_str:
            return "tRNA"
        
    # rRNA: long, high GC
    if length > 500 and gc > 50: 
        return "rRNA"
    
    # miRNA precursor: 60-120 nt, hairpin-like (palindromic)
    if 60 <= length <= 120: 
        half = length // 2 
        left = seq_str[:half]
        right = seq_str[half:][::-1]
        matches = sum(1 for a, b in zip(left, right) if a == b)
        if matches / half > 0.3:
            return "miRNA_precursor"
        
    # snRNA: 100-300 nt, U-rich (T-rich in DNA)
    if 100 <= length <= 300: 
        if seq_str.count("T") / length > 0.25: 
            return "snRNA"
        
    # lncRNA: >200 nt, low ORF content 
    if length > 200:
        return "lncRNA"

    return "unknown"

#codon usage 
def codon_usage(seq):
    """
    Returns a dictionary with counts of all 64 codons in the sequence. Only counts codons in frame 0 (simple usage table).
    """
    seq_str = str(seq).upper()
    codons = {}

    # Initialize all 64 codons to 0
    bases = ["A", "T", "G", "C"]
    for a in bases:
        for b in bases:
            for c in bases:
                codons[a+b+c] = 0

    # Count codons in frame 0
    for i in range(0, len(seq_str) - 2, 3):
        codon = seq_str[i:i+3]
        if codon in codons:
            codons[codon] += 1

    return codons
    
def codon_usage_summary(codon_table):
    total = sum(codon_table.values())
    if total == 0:
        return 0, 0, 0 #avoiding division by zero
    #most used codon 
    most_used = max(codon_table, key=codon_table.get)
    most_used_count = codon_table[most_used]

    #codon bias (max frequency / average frequency)
    avg = total / 64
    bias = most_used_count / avg if avg > 0 else 0

    return most_used, most_used_count, round(bias, 3)

# 6-frame ORF upgrade 
def reverse_complement(seq):
    seq = str(seq).upper()
    comp = str.maketrans("ATGC", "TACG")
    return seq.translate(comp)[::-1]


def find_orfs_6frame(seq):
    """Returns ORFs from all 6 reading frames: +0, +1, +2, -0, -1, -2. Each ORF is returned as a tuple: (frame, start, end, length)"""
    seq = str(seq).upper()
    rev = reverse_complement(seq)
    results = []

    #forward frames
    for frame in range(3): 
        orfs = find_orfs(seq[frame:])
        for start, end, length in orfs: 
            results.append((f"+{frame}", start + frame, end + frame, length))

    # reverse frames
    for frame in range(3): 
        orfs = find_orfs(rev[frame:])
        for start, end, length in orfs: 
            real_start = len(seq) - (start + frame) - length
            real_end = real_start + length
            results.append((f"-{frame}", real_start, real_end, length))

    return results     

#Motif detection to search for biological patterns 
MOTIFS = {"TATA_box": "TATAAA",
          "Kozak_consensus": "GCCACCATGG",
          "PolyA_signal": "AATAAA",
          "CpG_island": "CG",
          "Shine_Dalgarno": "AGGAGG"}
def find_motifs(seq, motifs=MOTIFS):
    seq_str = str(seq).upper()
    results = {}

    for name, pattern in motifs.items():
        positions = []
        start = 0

        while True: 
            idx = seq_str.find(pattern, start)
            if idx == -1:
                break
            positions.append(idx)
            start = idx + 1 
        results[name] = positions
            
    return results 

CODON_TABLE = {
    "ATA":"I","ATC":"I","ATT":"I","ATG":"M",
    "ACA":"T","ACC":"T","ACG":"T","ACT":"T",
    "AAC":"N","AAT":"N","AAA":"K","AAG":"K",
    "AGC":"S","AGT":"S","AGA":"R","AGG":"R",
    "CTA":"L","CTC":"L","CTG":"L","CTT":"L",
    "CCA":"P","CCC":"P","CCG":"P","CCT":"P",
    "CAC":"H","CAT":"H","CAA":"Q","CAG":"Q",
    "CGA":"R","CGC":"R","CGG":"R","CGT":"R",
    "GTA":"V","GTC":"V","GTG":"V","GTT":"V",
    "GCA":"A","GCC":"A","GCG":"A","GCT":"A",
    "GAC":"D","GAT":"D","GAA":"E","GAG":"E",
    "GGA":"G","GGC":"G","GGG":"G","GGT":"G",
    "TCA":"S","TCC":"S","TCG":"S","TCT":"S",
    "TTC":"F","TTT":"F","TTA":"L","TTG":"L",
    "TAC":"Y","TAT":"Y","TAA":"*","TAG":"*",
    "TGC":"C","TGT":"C","TGA":"*","TGG":"W",
}

def translate_sequence(seq):
    seq = seq.upper()
    aa = []
    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i+3]
        aa.append(CODON_TABLE.get(codon, "X"))  # X = unknown codon
    return "".join(aa)


# Kyte–Doolittle hydrophobicity scale
HYDRO_SCALE = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5
}

def hydrophobicity_profile(peptide, window=9):
    """Return sliding-window hydrophobicity values for a peptide."""
    values = [HYDRO_SCALE.get(aa, 0) for aa in peptide]
    profile = []

    for i in range(len(values)):
        start = max(0, i - window // 2)
        end = min(len(values), i + window // 2 + 1)
        window_avg = np.mean(values[start:end])
        profile.append(window_avg)

    return profile


#command interface to run on terminal locally
#if __name__ == "__main__": 
 #    parser = argparse.ArgumentParser(description="DNA Analyzer Tool")
 #
 #    #Required arguments
 #    parser.add_argument("input", help="Path to input FASTA file")
 #    parser.add_argument("output", help="Path to output CSV file")
 #
    #Parse arguments from the terminal 
 #    args = parser.parse_args()
 #
     #Run the analyzer with user=provided paths
 #    analyze_sequence(args.input, args.output)



     
    
