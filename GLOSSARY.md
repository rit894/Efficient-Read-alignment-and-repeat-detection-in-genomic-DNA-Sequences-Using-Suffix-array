# Glossary

| Term | Plain-language meaning |
|---|---|
| **Genome** | The full DNA sequence of an organism, written as a string over A, C, G, T. |
| **Base** | One letter of DNA (A, C, G or T). Genome length is counted in bases. |
| **Read** | A short fragment (50-300 bases) produced by a sequencing machine. |
| **Reference genome** | The known genome that reads are located in. |
| **Read alignment / mapping** | Finding where in the genome a read comes from. |
| **Pattern** | The string being searched for (a read, a gene, a motif). |
| **Motif** | A short, biologically meaningful pattern, e.g. `GAATTC` (EcoRI cut site). |
| **SNP / substitution** | A single base that differs between two sequences. |
| **Mismatch** | A position where the read and the genome disagree. |
| **Reverse complement** | The opposite DNA strand: reverse the string and swap A<->T, C<->G. |
| **Strand (+/-)** | Reads can come from either strand, so both orientations are searched. |
| **Repeat** | A substring that occurs more than once in a genome. |
| **Tandem repeat** | A short unit repeated back to back, e.g. `ATATATAT`. |
| **GC content** | Fraction of bases that are G or C. |
| **FASTA** | Text file format: a `>name` header line followed by the sequence. |
| **Suffix** | A string from position `i` to the end of the text. |
| **Sentinel `$`** | A unique end marker, smaller than every letter. |
| **Suffix array (SA)** | The start positions of all suffixes, sorted alphabetically by suffix. |
| **LCP** | Longest Common Prefix of two strings; the LCP array stores it for neighbouring suffixes. |
| **BWT** | Burrows-Wheeler Transform: the character before each sorted suffix. |
| **LF-mapping** | The rule that links a row's last-column character to its row in the first column; lets us walk backwards through the text. |
| **FM-index** | An index built on the BWT that supports pattern counting without the full text or SA. |
| **Occ(c, i)** | How many times character `c` appears in the first `i` characters of the BWT. |
| **C table** | For each character, how many characters in the text are smaller. |
| **Backward search** | Matching a pattern from its last character to its first by shrinking a row interval. |
| **Pigeonhole principle** | If a pattern is cut into k+1 pieces and has at most k errors, at least one piece is error-free. |
| **ADT** | Abstract Data Type: a description of operations and rules, not of the implementation. |
