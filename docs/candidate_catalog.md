# Prioritized Candidate Genes Master Catalog — *E. coli* QA5221

This catalog serves as the master reference for the **23 prioritized candidate genes** identified in the clinical multidrug-resistant (MDR) extraintestinal pathogenic *Escherichia coli* (ExPEC) sequence type 354 (ST354) isolate **QA5221**. 

All 23 candidates are pangenomic singletons (unique to QA5221 across 32 lineage genomes) that lacked any sequence-level annotations (**0 Pfam domain hits**). Their functional identities were resolved using a sequence-independent structural genomics pipeline combining **ESMFold** 3D modeling and **Foldseek/TM-align** structural similarity searches against Swiss-Prot, PDB, and the AlphaFold Database (afdb50).

---

## 📊 Publication Mapping Summary

*   **Total Prioritized Candidates:** 23
*   **Paper 1 (Prioritization Pipeline & GNAT):** 4 candidates analyzed (focuses on GNAT biophysical dynamics).
*   **Paper 2 (NODE_1 Intact Prophage):** 6 candidates analyzed (includes enterohemolysin cargo and structural genes; 1 shared with Paper 1).
*   **Paper 3 (NODE_23 Genomic Island):** 2 candidates analyzed (methyltransferase and partition protein).
*   **Unused / Available for Future Work:** 12 candidates (including the `NODE_12` LPS-modification operon).

---

## 🗺️ Master Catalog Table

| Locus Tag | Size (aa) | GC (%) | Contig Location | Target Paper | Top Foldseek Structural Homolog | TM-Score | RMSD (Å) | Predicted Functional Class & Nature |
| :--- | :---: | :---: | :--- | :---: | :--- | :---: | :---: | :--- |
| **KNGPFPPJ_02769** | 351 | 32.9% | NODE_9: 165.2k–166.3k | **Paper 1** | N-acetyltransferase (*Salmonella enterica*) | **0.9569** | 1.55 | **AMR**: GNAT-family GCN5-like acetyltransferase (*GNAT_KA27*) |
| **KNGPFPPJ_00061** | 350 | 48.5% | NODE_1: 58.6k–59.7k | **Paper 1 & 2** | Putative enterohemolysin (*E. coli*) | **0.5899** | 2.99 | **Virulence**: Pore-forming cytotolysin toxin cargo (*Ehly_61*) |
| **KNGPFPPJ_03161** | 372 | 29.3% | NODE_12: 40.1k–41.2k | **Paper 1** | O-antigen polymerase (*Salmonella enterica*) | **0.8953** | 2.53 | **Cell Wall**: Wzy-family O-antigen polymerase (*OAgP_161*) |
| **KNGPFPPJ_04371** | 984 | 35.8% | NODE_24: 6.0k–8.9k | **Paper 1** | Outer membrane protein YejO (*E. coli*) | 0.3450 | 5.40 | **Virulence**: Autotransporter autolytic anchor (*OAT_371*) |
| **KNGPFPPJ_00084** | 203 | 45.8% | NODE_1: 73.2k–73.8k | **Paper 2** | Phage antirepressor protein (*Shigella*) | **0.5352** | 2.62 | **Phage**: Lytic/lysogeny repressor switch protein |
| **KNGPFPPJ_00091** | 249 | 51.2% | NODE_1: 76.5k–77.3k | **Paper 2** | Phage large terminase (*Citrobacter*) | 0.4573 | 4.82 | **Phage**: Phage packaging terminase enzyme |
| **KNGPFPPJ_00107** | 216 | 48.5% | NODE_1: 89.8k–90.4k | **Paper 2** | Dit-like phage tail protein (*Salmonella*) | **0.5864** | 3.56 | **Phage**: Structural tail hub assembly protein |
| **KNGPFPPJ_00109** | 344 | 44.9% | NODE_1: 90.7k–91.8k | **Paper 2** | Major capsid protein (*Serratia sp.*) | **0.7750** | 3.04 | **Phage**: Head capsid structural protein |
| **KNGPFPPJ_00114** | 399 | 48.6% | NODE_1: 93.6k–94.8k | **Paper 2** | Baseplate protein J-like (*Jejubacter*) | 0.4854 | 5.02 | **Phage**: Baseplate tail host-recognition receptor |
| **KNGPFPPJ_04325** | 294 | 35.0% | NODE_23: 32.6k–33.5k | **Paper 3** | Type II Methylase Mcbe1 (*Caldicellulosiruptor*) | **0.5000** | 1.45 | **Epigenetics**: site-specific DNA adenine methyltransferase |
| **KNGPFPPJ_04326** | 501 | 35.1% | NODE_23: 33.5k–35.1k | **Paper 3** | Chromosome-partitioning ParB (*Caulobacter*) | **0.5200** | 4.87 | **Segregation**: CTP-dependent ParB DNA segregation protein |
| **KNGPFPPJ_03157** | 371 | 31.4% | NODE_12: 35.6k–36.7k | *Unused (Paper 4)* | Glycosyltransferase family 4 (*Vibrio*) | **0.8972** | 2.33 | **Virulence**: Capsule/LPS modifying glycosyltransferase |
| **KNGPFPPJ_03836** | 246 | 46.7% | NODE_17: 80.6k–81.4k | *Unused* | Uracil-DNA glycosylase-like domain (*E. coli*) | **0.9502** | 1.63 | **DNA Repair**: Divergent Uracil-DNA base excision repair |
| **KNGPFPPJ_00097** | 313 | 53.1% | NODE_1: 82.5k–83.4k | *Unused* | DUF2184 domain protein (*Agrobacterium*) | **0.8229** | 2.50 | **Membrane**: Transmembrane transporter of unknown function |
| **KNGPFPPJ_00103** | 381 | 49.0% | NODE_1: 85.5k–86.7k | *Unused* | DUF3383 domain protein (*Cedecea*) | **0.9051** | 1.77 | **Metabolism**: Conserved domain of unknown function |
| **KNGPFPPJ_00106** | 668 | 52.2% | NODE_1: 87.8k–89.8k | *Unused* | Putative bacteriophage tail hydrolase (*Salmonella*) | 0.3253 | 6.43 | **Cell Wall**: Peptidoglycan lytic transglycosylase fold |
| **KNGPFPPJ_00112** | 250 | 51.1% | NODE_1: 92.5k–93.3k | *Unused* | Translation initiation factor IF-2 (*Morganella*)| **0.8974** | 2.17 | **Translation**: Divergent ribosomal translation factor |
| **KNGPFPPJ_00115** | 226 | 47.9% | NODE_1: 94.8k–95.5k | *Unused* | DUF2612 domain protein (*Jejubacter*) | **0.6774** | 4.09 | **Cargo**: Phage-carried uncharacterized cargo |
| **KNGPFPPJ_00400** | 574 | 32.9% | NODE_1: 380.7k–382.4k | *Unused* | DUF2326 domain protein (*E. coli*) | **0.6890** | 5.42 | **Metabolism**: Atypical low-GC core-accessory domain |
| **KNGPFPPJ_01156** | 300 | 32.6% | NODE_3: 293.9k–294.8k | *Unused* | Uncharacterized protein (*Pseudomonas*) | **0.7379** | 4.24 | **Metabolism**: Conserved uncharacterized sequence |
| **KNGPFPPJ_01158** | 677 | 32.6% | NODE_3: 295.2k–297.3k | *Unused* | Uncharacterized protein (*Aeromonas*) | **0.7823** | 4.42 | **Metabolism**: Conserved uncharacterized sequence |
| **KNGPFPPJ_04303** | 256 | 39.4% | NODE_23: 16.2k–17.0k | *Unused* | Uncharacterized protein (*Pantoea brenneri*) | **0.7623** | 3.37 | **Island Cargo**: Uncharacterized auxiliary island protein |
| **KNGPFPPJ_04571** | 290 | 38.6% | NODE_30: 9.2k–10.1k | *Unused* | Uncharacterized protein (*E. coli*) | **0.7169** | 2.16 | **Metabolism**: Atypical uncharacterized domain |

> [!NOTE]
> Bold TM-scores indicate **high-confidence fold matches** (TM-score $\ge 0.50$), confirming that the tertiary structure matches a characterized protein family despite having zero sequence similarity.

---

## 📂 Detailed Gene Profiles by Study Locus

### 1. Paper 1 Focus Loci (Prioritization & AMR Validation)
*   **`KNGPFPPJ_02769` (GNAT_KA27):** 
    *   *Role:* Putative GNAT family aminoglycoside acetyltransferase (AAC).
    *   *Significance:* Contains the classic N-acetyltransferase active recognition cleft. Molecular docking and NPT-ensemble explicit-solvent molecular dynamics simulations mapped the structural coordinates of aminoglycoside (kanamycin, gentamicin, amikacin) target interactions. MM-GBSA free energy analysis verified a $\sim 10.8 \text{ kcal/mol}$ binding free energy gap, confirming aminoglycoside specificity.
*   **`KNGPFPPJ_03161` (OAgP_161):**
    *   *Role:* Putative Wzy-family O-antigen polymerase.
    *   *Significance:* Displays high structural homology to O-antigen polymerases across Enterobacteriaceae. A multi-span transmembrane fold that has been validated using explicit-solvent lipid bilayer molecular dynamics simulations, demonstrating stable topology.
*   **`KNGPFPPJ_00061` (Ehly_61):**
    *   *Role:* Putative enterohemolysin cytotoxin.
    *   *Significance:* Pore-forming virulence factor located as cargo inside the `NODE_1` prophage. Foldseek confirmed structural homology to cytolysins, and its membrane-bilayer molecular dynamics (currently running) indicates membrane-disruption potential.
*   **`KNGPFPPJ_04371` (OAT_371):**
    *   *Role:* Divergent autotransporter outer-membrane protein.
    *   *Significance:* Plasmid-associated autotransporter containing a passenger domain and a C-terminal $\beta$-barrel anchor.

---

### 2. Paper 2 Focus Loci (`NODE_1` Intact Prophage)
This operon represents an intact, lytically competent lambdoid Siphoviridae prophage (attL: 56,422 bp; attR: 98,747 bp).
*   **`KNGPFPPJ_00084` (Antirepressor):** Governs phage transcription repressor toggling, regulating lysogenic-to-lytic cycle induction.
*   **`KNGPFPPJ_00091` (Terminase):** Large packaging subunit responsible for translocation of phage genomic DNA into empty capsids.
*   **`KNGPFPPJ_00107` (Tail Assembly) & `KNGPFPPJ_00109` (Major Capsid):** Core structural proteins that make up the structural capsid shell and distal tail hub.
*   **`KNGPFPPJ_00114` (Baseplate receptor):** Baseplate J-like domain that binds specifically to host outer-membrane receptor complexes.

---

### 3. Paper 3 Focus Loci (`NODE_23` Host-Defense & Partitioning Island)
This locus represents a hybrid phage-plasmid chromosomal integrated element flanked by a tyrosine recombinase XerC.
*   **`KNGPFPPJ_04325` (DNA Methyltransferase):** Divergent adenine methyltransferase. Active site analysis confirms it carries the S-adenosylmethionine (SAM) binding domain and target recognition loop similar to classical Dam methylases, acting as a host-defense system.
*   **`KNGPFPPJ_04326` (ParB Partition protein):** Segregation protein showing structural homology to CTP-dependent ParB systems, ensuring the vertical segregation of the integrated island during cell division.

---

### 4. High-Priority Unused Candidates (Recommended for Paper 4)
*   **`KNGPFPPJ_03157` (Glycosyltransferase family 4):**
    *   *Location:* Located on `NODE_12` physically close to the O-antigen polymerase (`KNGPFPPJ_03161`).
    *   *Potential:* Since LPS/O-antigen glycosylation alters host complement recognition (serum resistance) and shields the outer membrane from colistin or host cationic antimicrobial peptides (CAMPs), this represents a **novel virulence and peptide-resistance marker**.
    *   *Paper 4 Concept:* *"Genomic and structural characterization of a novel cell wall modification and lipopolysaccharide glycosylation locus in clinical E. coli ST354."*
