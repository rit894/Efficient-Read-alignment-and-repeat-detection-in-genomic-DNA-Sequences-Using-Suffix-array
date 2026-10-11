# Presentation Outline and Demo Script

Suggested length: 12 to 14 slides, 10 to 15 minutes + 3 minutes live demo.

## Slide plan

| # | Slide | Content | Visual |
|---|---|---|---|
| 1 | Title | Project name, team, course | DNA helix image |
| 2 | The problem | Millions of reads, genomes of billions of bases; scanning is too slow | "Read -> where in the genome?" sketch |
| 3 | Why a naive scan fails | O(n) per query x millions of queries | `results/benchmark.png` (left plot) |
| 4 | Idea: index the text once | Build once, query fast | Diagram: text -> index -> queries |
| 5 | Suffix array | Sorted suffixes of `banana$`, binary search | Table from `docs/ADT_SPECIFICATION.md` |
| 6 | LCP array | Neighbouring suffix overlap; longest repeat = max LCP | `banana$` table with LCP column |
| 7 | BWT and FM-index | Preceding-character column, backward search, sampled SA | Step trace from `python main.py demo` |
| 8 | ADTs | Operations and invariants of each structure | Dependency diagram (SA -> LCP, BWT -> FM) |
| 9 | Applications | Repeats, genome comparison, approximate search, read alignment (both strands) | Pigeonhole picture (read cut in k+1 pieces) |
| 10 | Complexity | Build, query, space tables | Tables from `docs/COMPLEXITY_ANALYSIS.md` |
| 11 | Experiments | Synthetic genomes 1k to 200k bases, 100 queries, correctness check | `results/benchmark.png` |
| 12 | Results | Index search about 1,500x (SA) and 650x (FM) faster than scanning at 200k bases; FM stores 6.25 % of SA | Results table |
| 13 | Limitations and future work | Python speed, substitutions only, SA-IS, wavelet trees | Bullet list |
| 14 | Conclusion and Q&A | Take-aways | |

## Live demo script (about 3 minutes)

```bash
python main.py demo                      # banana table, FM backward-search trace, DNA results
python main.py search --fasta data/sample_genome.fasta --pattern GAATTC --method kmp
python main.py search --fasta data/sample_genome.fasta --pattern GAATTC --method sa
python main.py align  --fasta data/sample_genome.fasta --reads 50 --mismatch 2
python main.py test                      # 17 tests pass
```

Talking points while it runs:
* "Same answers from every method, only the speed differs."
* "The FM-index trace shows the row interval shrinking one character at a time."
* "Reads have errors and can come from either strand, so we search both orientations
  and allow mismatches using the pigeonhole principle."

## Questions you may be asked (and short answers)

| Question | Answer |
|---|---|
| Why not just use KMP? | KMP is O(n) per query. With millions of reads the genome is rescanned millions of times; an index is built once. |
| Why a sentinel `$`? | It makes every suffix unique and gives the BWT a well-defined start/end. |
| Why is the FM-index slower than the SA here? | It does extra counting work per query to save memory; its value shows at genome scale where the SA does not fit in RAM. |
| What does Kasai's trick save? | LCP of neighbouring suffixes is computed in O(n) instead of O(n^2) because the LCP drops by at most 1 between consecutive text positions. |
| How is the SA built? | Prefix doubling: sort by first 1, 2, 4, 8... characters using the previous ranks. O(n log^2 n) here. |
| How are mismatches handled? | Cut the read into k+1 pieces; one piece must match exactly; find it with the index, then verify the whole read. |
| What would you improve? | SA-IS construction in a compiled language, indel support, wavelet-tree BWT. |
