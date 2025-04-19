import os
import ssl
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import Bio.SeqUtils
from Bio.SeqUtils import MeltingTemp as mt

# For SSL issues (if needed)
if not os.environ.get("PYTHONHTTPSVERIFY", "") and getattr(ssl, "_create_unverified_context", None):
    ssl._create_default_https_context = ssl._create_unverified_context

# ---------------------------------------------------------------------
# 1. fetch_gene_by_gb_id
# ---------------------------------------------------------------------
def fetch_gene_by_gb_id(genbank_id, gene_name):
    """
    Given a GenBank accession (e.g. 'NC_000913.3') and a gene name (e.g. 'metB'),
    retrieve the full GenBank record from NCBI and extract:
      - The entire genomic sequence
      - The exact gene sequence (5'->3')
      - The protein sequence (for the matching CDS)
    Returns a dict with:
      {
        "genomic_accession": str,
        "Whole genome": str,
        "gene_sequence": str,
        "gene_start": int,
        "gene_end": int,
        "protein_sequence": str
      }
    """
    Entrez.email = "cypete.cp@gmail.com"
    # Entrez.api_key = "17da3*****************" 

    handle = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="gb", retmode="text")
    gb_record = SeqIO.read(handle, "gb")
    handle.close()

    # Full genome
    full_seq = gb_record.seq
    seq_len = len(full_seq)

    # Find the named gene feature
    gene_feature = None
    for feat in gb_record.features:
        if feat.type == "gene":
            if "gene" in feat.qualifiers and gene_name.lower() in [g.lower() for g in feat.qualifiers["gene"]]:
                gene_feature = feat
                break

    if not gene_feature:
        raise ValueError(f"Gene '{gene_name}' not found in GenBank record '{genbank_id}'.")

    strand = gene_feature.location.strand
    gene_start = int(gene_feature.location.start)
    gene_end = int(gene_feature.location.end)

    if gene_start < 0 or gene_end > seq_len or gene_start >= gene_end:
        raise ValueError(f"Invalid gene range {gene_start}-{gene_end} for genome length {seq_len}.")

    gene_seq = full_seq[gene_start:gene_end]
    if strand == -1:
        gene_seq = gene_seq.reverse_complement()

    # Find the corresponding protein sequence
    protein_seq = ""
    for feat in gb_record.features:
        if feat.type == "CDS":
            if "gene" in feat.qualifiers and gene_name.lower() in [g.lower() for g in feat.qualifiers["gene"]]:
                protein_seq = feat.qualifiers.get("translation", [""])[0]
                break

    return {
        "genomic_accession": gb_record.id,
        "Whole genome": str(full_seq),       
        "gene_sequence": str(gene_seq),
        "gene_start": gene_start,
        "gene_end": gene_end,
        "protein_sequence": protein_seq
    }

# ---------------------------------------------------------------------
# 2. Helper functions
# ---------------------------------------------------------------------
def check_sequence(protein_sequence, position, residue):
    """
    Checks if the protein has the correct residue at a given position (1-based).
    Raises ValueError if there's a mismatch.
    """
    p_index = position - 1
    if len(protein_sequence) < position:
        raise ValueError(
            f"Protein is only {len(protein_sequence)} residues long; position {position} is out of range."
        )
    actual_aa = protein_sequence[p_index]
    if actual_aa != residue:
        raise ValueError(
            f"Residue mismatch at position {position}: expected '{residue}', found '{actual_aa}'."
        )

def reverse_complement(seq_str):
    """
    Returns the reverse complement of a DNA sequence string.
    """
    return str(Seq(seq_str).reverse_complement())
# Helper to pick a codon for a given amino acid.
def pick_codon_for_aa(aa):
    aa_to_codons = {
        'A': ['GCG','GCC','GCA','GCT'],
        'C': ['TGC','TGT'],
        'D': ['GAT','GAC'],
        'E': ['GAA','GAG'],
        'F': ['TTT','TTC'],
        'G': ['GGC','GGT','GGG','GGA'],
        'H': ['CAT','CAC'],
        'I': ['ATT','ATC','ATA'],
        'K': ['AAA','AAG'],
        'L': ['CTG','TTA','TTG','CTT','CTC','CTA'],
        'M': ['ATG'],
        'N': ['AAT','AAC'],
        'P': ['CCG','CCA','CCT','CCC'],
        'Q': ['CAG','CAA'],
        'R': ['CGT','CGC','CGA','CGG','AGA','AGG'],
        'S': ['AGC','TCT','AGT','TCC','TCA','TCG'],
        'T': ['ACC','ACG','ACT','ACA'],
        'V': ['GTG','GTT','GTC','GTA'],
        'W': ['TGG'],
        'Y': ['TAT','TAC'],
    }
    aa = aa.upper()
    possible = aa_to_codons.get(aa, [])
    if not possible:
        raise ValueError(f"No known codon for amino acid '{aa}'.")
    return possible[0]  # pick first for simplicity
# ---------------------------------------------------------------------
# 3. Design Partner Primers
# ---------------------------------------------------------------------
def design_partner_primers(genbank_id, gene_name,upstream_rev, downstream_fwd,
                           flank_distance=500, tolerance=100,
                           primer_length=25):
    """
    Designs partner primers ~500 ± 100 bp away from gene boundaries, oriented toward the gene.
    Returns upstream_partner_primer, downstream_partner_primer.
    """
    fetch_output = fetch_gene_by_gb_id(genbank_id, gene_name)
    genome_seq = fetch_output['Whole genome']
    gene_start = fetch_output['gene_start']
    gene_end = fetch_output['gene_end']

    genome_len = len(genome_seq)
    tm_up = mt.Tm_NN(upstream_rev)
    tm_down = mt.Tm_NN(downstream_fwd)

    # Upstream window
    up_min = gene_start - (flank_distance + tolerance)
    up_max = gene_start - (flank_distance - tolerance)
    if up_min < 0: up_min = 0
    if up_max < 0:
        raise ValueError("No valid upstream window; gene is near start of genome.")

    def attempt_forward_primer_3prime(pos):
        start_idx = pos - primer_length
        end_idx = pos
        if start_idx < 0:
            return None
        p_seq = genome_seq[start_idx:end_idx]
        tm_val = mt.Tm_NN(p_seq)
        gc_val = Bio.SeqUtils.gc_fraction(p_seq)
        if (tm_up - 20) <= tm_val <= (tm_up + 20) and gc_val <= 0.60:
            return str(p_seq).upper(), (start_idx, end_idx)
        return None

    upstream_partner = None
    for pos in range(up_max, up_min - 1, -1):
        if pos < primer_length:
            continue
        candidate = attempt_forward_primer_3prime(pos)
        if candidate:
            upstream_partner = candidate
            break

    if not upstream_partner:
        raise ValueError("No suitable upstream partner primer found.")
    upstream_partner_primer, _ = upstream_partner

    # Downstream window
    down_min = gene_end + (flank_distance - tolerance)
    down_max = gene_end + (flank_distance + tolerance)
    if down_max > genome_len:
        down_max = genome_len

    def attempt_reverse_primer_3prime(pos):
        end_idx = pos + primer_length
        if end_idx > genome_len:
            return None
        p_seq = genome_seq[pos:end_idx]
        rc_p_seq = reverse_complement(p_seq)
        tm_val = mt.Tm_NN(rc_p_seq)
        gc_val = Bio.SeqUtils.gc_fraction(rc_p_seq)
        if (tm_down - 20) <= tm_val <= (tm_down + 20) and gc_val <= 0.60:
            return rc_p_seq.upper(), (pos, end_idx)
        return None

    downstream_partner = None
    for pos in range(down_min, down_max + 1):
        candidate = attempt_reverse_primer_3prime(pos)
        if candidate:
            downstream_partner = candidate
            break

    if not downstream_partner:
        raise ValueError("No suitable downstream partner primer found.")
    downstream_partner_primer, _ = downstream_partner

    return {
        "upstream_partner_primer": upstream_partner_primer,
        "downstream_partner_primer": downstream_partner_primer
    }


# ---------------------------------------------------------------------
# 4. Design Knockout Primers
# ---------------------------------------------------------------------
def design_knockout_primers(genbank_id, gene_name,
                            upstream_anneal_len=25,
                            downstream_anneal_len=25,
                            overlap_len=20, flank_distance = 500 , tolerance =100,
                           partner_primer_length =25):
    """
    Designs two-part primers for a gene knockout.
    Returns a dict with upstream_rev_homolog_primer, downstream_fwd_homolog_primer, etc.
    """
    fetch_output = fetch_gene_by_gb_id(genbank_id, gene_name)
    genome_seq = fetch_output['Whole genome']
    gene_start = fetch_output['gene_start']
    gene_end = fetch_output['gene_end']
    # Upstream annealing
    upstream_anneal_region = str(genome_seq[gene_start - upstream_anneal_len: gene_start])
    # Downstream annealing
    downstream_anneal_region = str(genome_seq[gene_end: gene_end + downstream_anneal_len])

    if overlap_len > len(downstream_anneal_region) or overlap_len > len(upstream_anneal_region):
        raise ValueError("Overlap length cannot exceed the annealing region lengths.")

    # Overlaps
    upstream_overlap_region = reverse_complement(downstream_anneal_region[:overlap_len])
    downstream_overlap_region = upstream_anneal_region[-overlap_len:]

    # Construct primers
    upstream_anneal_rc = reverse_complement(upstream_anneal_region)
    downstream_anneal_fwd = downstream_anneal_region

    upstream_primer = (upstream_overlap_region + upstream_anneal_rc).upper()
    downstream_primer = (downstream_overlap_region + downstream_anneal_fwd).upper()

    output = design_partner_primers(genbank_id, gene_name,upstream_primer, downstream_primer,
                           flank_distance, tolerance,
                           partner_primer_length)
    upstream_fwd_primer = output['upstream_partner_primer']
    downstream_rv_primer = output ['downstream_partner_primer']

    return {
        "upstream_rev_homolog_primer": upstream_primer,
        "downstream_fwd_homolog_primer": downstream_primer,
        'upstream_fwd_primer' : upstream_fwd_primer,
        'downstream_rv_primer': downstream_rv_primer,
        "upstream_anneal_region": upstream_anneal_region.upper(),
        "downstream_anneal_region": downstream_anneal_region.upper(),
        "upstream_overlap_region": upstream_overlap_region.upper(),
        "downstream_overlap_region": downstream_overlap_region.upper()
    }


# ---------------------------------------------------------------------
# 5. Multiple mutations: design_single_mutation_primers, design_combined_primers, design_multiple_knockin_primers
# ---------------------------------------------------------------------


def design_multiple_knockin_primers(genbank_id, gene_name, mutations,
                                    left_flank_length=10, right_flank_length=10,
                                    merge_distance=48):
    """
    1) Creates individual primers for each mutation.
    2) Groups mutations if their intervals are <= merge_distance apart.
    3) Produces combined primers for grouped sets of mutations.
    """

    fetch_output = fetch_gene_by_gb_id(genbank_id, gene_name)
    dna_sequence = fetch_output['gene_sequence']
    protein_sequence = fetch_output['protein_sequence']

    def design_single_mutation_primers(genbank_id, gene_name,
                                   position, original_residue, desired_residue,
                                   left_flank_length, right_flank_length):
        """
        Designs forward/reverse primers for one site-directed mutation.
        Also returns region_start/region_end for interval usage.
        """

        check_sequence(protein_sequence, position, original_residue)
        codon_index_start = (position - 1) * 3
        codon_index_end = codon_index_start + 3
        new_codon = pick_codon_for_aa(desired_residue)

        left_flank_start = max(0, codon_index_start - left_flank_length)
        right_flank_end = min(len(dna_sequence), codon_index_end + right_flank_length)

        left_flank_seq = dna_sequence[left_flank_start:codon_index_start]
        right_flank_seq = dna_sequence[codon_index_end:right_flank_end]

        mutated_segment = left_flank_seq + new_codon + right_flank_seq
        forward_primer = mutated_segment
        reverse_primer = reverse_complement(mutated_segment)

        return {
            "mutation_position": position,
            "original_aa": original_residue,
            "desired_aa": desired_residue,
            "forward_primer": forward_primer,
            "reverse_primer": reverse_primer,
            "region_start": left_flank_start,
            "region_end": right_flank_end 
        }
    def design_combined_primers(genbank_id, gene_name, mutations,
                            left_flank_length=10, right_flank_length=10):
        """
        Given multiple mutations in close proximity, designs a single pair of primers
        that encodes all changes in one region.
        """

        # Confirm each original residue
        for m in mutations:
            check_sequence(protein_sequence, m["position"], m["original_aa"])

        positions = [m["position"] for m in mutations]
        left_most_pos = min(positions)
        right_most_pos = max(positions)

        codon_index_start = (left_most_pos - 1) * 3
        codon_index_end = (right_most_pos - 1) * 3 + 3

        left_flank_start = max(0, codon_index_start - left_flank_length)
        right_flank_end = min(len(dna_sequence), codon_index_end + right_flank_length)

        # Rebuild the region codon by codon, applying new codons
        pos_to_codon = {}
        for m in mutations:
            cstart = (m["position"] - 1) * 3
            new_c = pick_codon_for_aa(m["desired_aa"])
            pos_to_codon[cstart] = new_c

        mutated_region = []
        for i in range(codon_index_start, codon_index_end, 3):
            orig_triplet = dna_sequence[i:i+3]
            if i in pos_to_codon:
                mutated_region.append(pos_to_codon[i])
            else:
                mutated_region.append(orig_triplet)
        mutated_region_str = "".join(mutated_region)

        left_flank_seq = dna_sequence[left_flank_start:codon_index_start]
        right_flank_seq = dna_sequence[codon_index_end:right_flank_end]

        mutated_segment = left_flank_seq + mutated_region_str + right_flank_seq
        forward_primer = mutated_segment
        reverse_primer = reverse_complement(mutated_segment)

        return {
            "mutations_covered": [(m["position"], m["original_aa"], m["desired_aa"]) for m in mutations],
            "forward_primer": forward_primer,
            "reverse_primer": reverse_primer,
            "region_start": left_flank_start,
            "region_end": right_flank_end
        }
     
    # Make individual designs
    individual = []
    for mut in mutations:
        single = design_single_mutation_primers(
            dna_sequence, protein_sequence,
            position=mut["position"],
            original_residue=mut["original_aa"],
            desired_residue=mut["desired_aa"],
            left_flank_length=left_flank_length,
            right_flank_length=right_flank_length
        )
        individual.append(single)

    # Sort by region_start
    individual.sort(key=lambda x: x["region_start"])

    # We'll track intervals plus the actual mutation data
    intervals_with_mutdata = []
    for i, mut in enumerate(mutations):
        start_ = individual[i]["region_start"]
        end_ = individual[i]["region_end"]
        intervals_with_mutdata.append((start_, end_, [mut]))

    intervals_with_mutdata.sort(key=lambda x: x[0])

    merged_intervals = []
    current_start, current_end, current_muts = intervals_with_mutdata[0]

    for i in range(1, len(intervals_with_mutdata)):
        next_start, next_end, next_muts = intervals_with_mutdata[i]
        # If next is within merge_distance, merge
        if next_start - current_end <= merge_distance:
            current_end = max(current_end, next_end)
            current_muts.extend(next_muts)
        else:
            # finalize the current group
            merged_intervals.append((current_start, current_end, current_muts))
            current_start, current_end, current_muts = next_start, next_end, next_muts

    merged_intervals.append((current_start, current_end, current_muts))

    combined = []
    for start, end, mut_list in merged_intervals:
        if len(mut_list) > 1:
            # Design combined primer for these multiple changes
            combined_primers = design_combined_primers(
                dna_sequence, protein_sequence, mut_list,
                left_flank_length=left_flank_length,
                right_flank_length=right_flank_length
            )
            combined.append(combined_primers)

    return {
        "individual_primers": individual,
        "combined_primers": combined
    }