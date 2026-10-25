# phiX174 reference dataset

The bundled data/phiX174_NC_001422.1.fasta contains the real NCBI reference sequence for Escherichia phage phiX174.

- Accession and version: NC_001422.1
- Sequence length: 5,386 bases
- Alphabet: A, C, G, T
- Retrieved: 25 October 2026
- NCBI record: https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1
- Download URL: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text
- Downloaded FASTA SHA-256: caa155c83967cc080f4f188f14f1f28c18d77448e87c3d00d3579ed348f4de2a
- Sequence-only SHA-256 (uppercase bases, no header or whitespace): 97038c7e1edea2297667d7f0426ba942b322c74cb30e072ec66ba47f9c0448d0

## Interpretation

The demo loads the reference without modifying its bases or planting repeats. The comparison copy and short reads are simulated: the copy uses a 1% substitution probability per base, and the reads use a 2% probability per base. These probabilities do not guarantee an exact percentage of changes. Seeds 9 and 10 make the simulations reproducible.

The implementation treats the stored sequence as linear. It does not find matches crossing from the end of the sequence back to the beginning. Read alignment supports substitutions and both orientations, but not insertions or deletions. Positions are zero-based.

## Example results

With the bundled reference and query ACGT: 5,386 bases, 44.76% GC content, 19 exact query matches, and a longest repeated substring of length 12 at positions 263 and 2759. These are sequence-analysis results, not evidence that a repeat has a particular biological function.

Run the demo with --fasta to select this reference. Without --fasta, the demo continues to generate the synthetic 3,000-base genome.
