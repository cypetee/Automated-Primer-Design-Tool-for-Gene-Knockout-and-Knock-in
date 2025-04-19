import streamlit as st
from Bio import Entrez
from primer_design_ import (
    design_knockout_primers,
    design_multiple_knockin_primers
)
api_key = st.secrets["api"]["key"]
st.set_page_config(page_title="Primer Design Tool", layout="wide")
st.title("🧬 Primer Design Web Tool")


mode = st.selectbox("Select Primer Design Mode", [
    "Knock-out Primers",
    "Knock-in Primers"])

with st.sidebar:
    if mode == "Knock-out Primers":
         st.markdown("## Knock-out Primer Design")

         st.markdown(
            """
Use this mode to design primers for gene knockouts in *E. coli* using Overlap PCR.  
The tool generates:

- Two **homologous primers** flanking the target gene  
- Two **partner primers** for amplifying upstream and downstream regions  

---

**Steps:**

1. **Enter a valid GenBank Accession ID** and **gene name**  
2. **Adjust the following parameters as needed:**

   - **Upstream/Downstream Anneal Length**: Region adjacent to the gene for primer binding  
   - **Overlap Length**: Homologous region shared between primers  
   - **Flank Distance**: Distance from gene to partner primer site (default: 500 bp)  
   - **Tolerance**: Flexibility range for partner primer position (±100 bp)  
   - **Partner Primer Length**: Length of primers used to amplify flanks  

--- 
Partner primers are complementary to the homologous primers and positioned ~500 bp from the gene, with a ±100 bp tolerance for optimal placement.
        
        """)
    elif mode == "Knock-in Primers":
        st.sidebar.markdown("## Knock-in Primer Design")

        st.sidebar.markdown(
            """
Use this mode to design primers for introducing single or multiple point mutations into *E. coli* genes using Overlap PCR.  
The tool checks the protein sequence to ensure that the specified positions are valid and designs forward and reverse primers of the mutation region.

The tool automatically generates:

- Individual **mutation-specific primers**  
- **Combined primers** if multiple mutations are within 48 bp of each other (to enable a single primer carrying multiple changes)  

---

**Steps:**

1. **Enter a valid GenBank Accession ID** and **target gene name**  
2. **Input desired mutations** in the format:  
   `Q58K, T75S` (i.e., OriginalResidue + Position + NewResidue)  
3. **Adjust parameters as needed:**

   - **Primer Length**: Total length of the designed mutation primer (default: 60 bp max for combined mutations)  
   - **Minimum Spacing**: Distance threshold (e.g., 48 bp) to merge nearby mutations into a single primer  

--- 
The tool validates that the original amino acid matches the reference protein sequence before introducing the mutation(s). If multiple mutations are close together, a combined primer is designed to carry all mutations in a single oligo.
            """
        )


st.markdown("---")

if mode == "Knock-out Primers":
    st.subheader("Design Knock-out Primers")
    genbank_id = st.text_input("GenBank ID", value="U00096.3", key="ko_gb")
    gene_name = st.text_input("Gene Name", value="metB", key="ko_gene")
    up_len = st.number_input("Upstream Anneal Length", value=25)
    down_len = st.number_input("Downstream Anneal Length", value=25)
    overlap = st.number_input("Overlap Length", value=20)
    flank_distance = st.number_input("Flank Length", value=500)
    tolerance = st.number_input("Tolerance", value=100)
    partner_primer_length = st.number_input("Partner Primer Length", value=20)

    if st.button("Design KO Primers"):
        try:
            primers = design_knockout_primers(genbank_id, gene_name, up_len, down_len, overlap,flank_distance, tolerance, partner_primer_length )
            st.json(primers)
        except Exception as e:
            st.error(str(e))



elif mode == "Knock-in Primers":
    st.subheader("Design Knock-in Primers")
    genbank_id = st.text_input("GenBank ID", value="U00096.3", key="multi_gb")
    gene_name = st.text_input("Gene Name", value="metB", key="multi_gene")
    mut_str = st.text_area("Enter Mutations (e.g. Q58K, T75S)", value="Q58K\nT75S")
    merge_dist = st.number_input("Merge Distance (bp)", value=48)

    def parse_mutations(mut_str):
        mut_lines = mut_str.strip().splitlines()
        muts = []
        for m in mut_lines:
            orig, pos, new = m[0], int(m[1:-1]), m[-1]
            muts.append({"position": pos, "original_aa": orig, "desired_aa": new})
        return muts

    if st.button("Design KI Primers"):
        try:
            mutations = parse_mutations(mut_str)
            primers = design_multiple_knockin_primers(genbank_id, gene_name, mutations, merge_distance=merge_dist)
            st.subheader("Individual Primers")
            st.json(primers["individual_primers"])
            st.subheader("Combined Primers")
            st.json(primers["combined_primers"])
        except Exception as e:
            st.error(str(e))
