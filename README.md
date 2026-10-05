# RiboseQC, repaired for current R

A copy of [RiboseQC](https://github.com/lcalviell/Ribo-seQC) (commit `152fa66`,
the last one, June 2020) with the changes needed to run it on a current R, to
run several samples at the same time, and to get the P-site offsets right when
ribosomes do not pile up on start codons.

## What differs from the original

| change | where | why |
|---|---|---|
| `disjointExons(annotation, aggregateGenes = T)` replaced by `exonicParts(annotation, linked.to.single.gene.only = FALSE)` | `prepare_annotation_files` | Bioconductor removed `disjointExons`, so the annotation step stopped with an error. The value is only stored in the `_Rannot` file; nothing reads it. |
| Closing code-fence markers no longer indented (six lines per template) | `rmd/riboseqc_template.Rmd`, `rmd/riboseqc_template_full.Rmd` | Newer rmarkdown/knitr does not recognise the indented form. |
| The report is rendered from a private copy of the template; new argument `tmp_dir` | `create_html_report`, `RiboseQC_analysis` | See below. From version 0.99.1. |
| P-site offsets ("cutoffs") from where the 5'ends switch on instead of where they peak, and chosen for all read lengths together; new arguments `cutoff_method` and `max_cutoff_diff` | `calc_cutoffs_from_profiles`, new `constrain_cutoffs_readlengths`, `RiboseQC_analysis` | See below. From version 0.99.2. **Changes results.** |

The version number differs from the original (0.99.0) so that this copy can be
told apart from it: 0.99.1 from the third change on, 0.99.2 from the fourth.

## The report step and `tmp_dir`

rmarkdown writes an intermediate file next to the file it renders. The original
rendered the template where it sits in the installed package, so every R
session wrote the same intermediate file. With several samples analysed at the
same time, one report could end up with another sample's content, or fail.

`create_html_report` now copies the template into a new subfolder of `tmp_dir`,
renders that copy, and removes the subfolder when it is done, also after an
error. Nothing else in `tmp_dir` is touched, and `tmp_dir` is created if it
does not exist.

`tmp_dir` defaults to the temporary folder of the R session, which is private
to that session. `RiboseQC_analysis` has the same argument and passes it on. A
pipeline that gives each job its own temporary folder can pass that folder; the
default is safe without it.

## P-site offsets

The offset of a read length is the distance from the 5'end of its reads to the
P-site. RiboseQC estimates it per read length from the 5'ends in the 24 nt
before start codons, over the most covered transcripts.

**What the original does.** It takes the position with the most 5'ends in that
window (the peak), found in several ways that then vote. That assumes
ribosomes on the start codon outnumber those on the next codons. It holds when
cells are pre-treated with cycloheximide or harringtonine. It does not hold in
general, for instance for frozen tissue lysed in the presence of cycloheximide:
there codons 2 to 5 carry as many reads as the start codon, the highest
position can be any of them, and the offset comes out too small by a whole
number of codons (0 instead of 12 for the most abundant read length in the data
this was found in). Periodicity stays perfect, so it is easy to miss.

**What this version does** (`cutoff_method = "edge"`, the default):

- *The signal is the rise, not the height.* At each position: the coverage
  minus the mean coverage of the two positions 3 and 6 nt before it, divided by
  the square root of their sum. The 5'ends switch on at the offset whatever
  happens after it, so the largest rise is at the offset. The design of the
  original is kept: the signal is read from the aggregate profile and from the
  transcripts one by one, with and without k-means, all positions or in-frame
  only, and the in-frame readings vote. The k-means now starts from fixed
  centers, so the result no longer varies between runs.
- *Read lengths are decided together* (`max_cutoff_diff = 1`). Read lengths
  that differ by 1 nt cannot have offsets that differ by more than about 1 nt.
  The possible offsets of one read length are 3 nt apart (its frame), so a
  neighbour leaves only one of them possible. The chain of offsets with the most
  support over all read lengths is chosen; read lengths with many reads and a
  clear frame weigh the most. This is not in the original.

`cutoff_method = "peak"` gives the original behaviour, unchanged.
`max_cutoff_diff = NA` decides each read length on its own. If a library really
has neighbouring read lengths whose offsets differ by more than `max_cutoff_diff`,
the joint choice forces them together; `cutoff_chain_margin` in the results is
then low.

The results keep the fields of the original (`final_cutoff`, `cutoffs`,
`frames`, `profiles_start`) and add, per read length, `final_cutoff_independent`,
`cutoff_support`, `cutoff_constraint` and `cutoff_chain_margin`.

**Checking existing results.** `recalculate_cutoffs.R` compares the offsets
stored in finished `*_results_RiboseQC` files with what the edge method gives,
without touching the BAM files or the result files:

```
Rscript recalculate_cutoffs.R /path/to/*/*_results_RiboseQC > cutoffs_old_vs_new.tsv
```

It reads the functions from `source/`, so it works before the package is rebuilt.

## Layout

This folder is the **installed** package: it can sit in an R library as
`RiboseQC` and be loaded as it is. R loads the built files, not the source.

```
DESCRIPTION, NAMESPACE, R/, Meta/, help/, html/, INDEX    the built package
rmd/                                                      the report templates, read as they are
source/RiboseQC/                                          the source the built package is made from
rebuild_from_source.sh                                    rebuilds the built package from the source
recalculate_cutoffs.R                                     compares stored offsets with those of the edge method
```

## Changing the package

- **A template:** edit it in `source/RiboseQC/inst/rmd/` and rebuild. (Editing
  `rmd/` directly also works, but the next rebuild overwrites it.)
- **A function:** edit `source/RiboseQC/R/riboseqc.R` and rebuild. Editing the
  source alone changes nothing, because R loads `R/RiboseQC.rdb`.

Rebuild with the R that will use the package, and with no jobs running that
use RiboseQC:

```
./rebuild_from_source.sh
```

The script builds in a temporary folder and replaces the built files here only
if the build succeeds. Then commit the source and the rebuilt files together.

The built files carry the R version they were built with (`Built:` in
`DESCRIPTION`). Build them on the machine that runs the pipeline.
