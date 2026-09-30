# RiboseQC, repaired for current R

A copy of [RiboseQC](https://github.com/lcalviell/Ribo-seQC) (commit `152fa66`,
the last one, June 2020) with the changes needed to run it on a current R and
to run several samples at the same time.

## What differs from the original

| change | where | why |
|---|---|---|
| `disjointExons(annotation, aggregateGenes = T)` replaced by `exonicParts(annotation, linked.to.single.gene.only = FALSE)` | `prepare_annotation_files` | Bioconductor removed `disjointExons`, so the annotation step stopped with an error. The value is only stored in the `_Rannot` file; nothing reads it. |
| Closing code-fence markers no longer indented (six lines per template) | `rmd/riboseqc_template.Rmd`, `rmd/riboseqc_template_full.Rmd` | Newer rmarkdown/knitr does not recognise the indented form. |
| The report is rendered from a private copy of the template; new argument `tmp_dir` | `create_html_report`, `RiboseQC_analysis` | See below. From version 0.99.1. |

The version number is 0.99.1 from the third change on, so that this copy can be
told apart from the original, which is 0.99.0.

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

## Layout

This folder is the **installed** package: it can sit in an R library as
`RiboseQC` and be loaded as it is. R loads the built files, not the source.

```
DESCRIPTION, NAMESPACE, R/, Meta/, help/, html/, INDEX    the built package
rmd/                                                      the report templates, read as they are
source/RiboseQC/                                          the source the built package is made from
rebuild_from_source.sh                                    rebuilds the built package from the source
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
