#!/usr/bin/env Rscript
# Compare the P-site cutoffs stored in existing RiboseQC results with the cutoffs the "edge" method gives.
#
#   Rscript recalculate_cutoffs.R <sample>_results_RiboseQC [more result files ...] > cutoffs_old_vs_new.tsv
#
# Nothing is rerun on the BAM files and no result file is changed. For every read length the stored 5'end profiles
# around start codons (25 nt, one row per transcript) and the stored frame of each transcript are put back into the
# layout calc_cutoffs_from_profiles() expects, and the functions are taken from source/RiboseQC/R/riboseqc.R next to
# this script, so the installed package does not have to be rebuilt first.
#
# One difference with a real run: the stored profiles start 24 nt before the start codon, a real run also has
# position -25 as reference. This can only matter for a cutoff of 23 or 24.
#
# Columns: sample, compartment, read_length, reads_start (5'ends in the 25 nt), frame_preference, selected (used for
# P-sites in the stored run), cutoff_stored, cutoff_edge (edge method, read length on its own), cutoff_edge_joint
# (edge method, read lengths together), constraint (chain/imputed), chain_margin, changed (stored != joint).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("give one or more *_results_RiboseQC files")
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
fns <- new.env()
for (x in parse(file.path(here, "source", "RiboseQC", "R", "riboseqc.R"), keep.source = FALSE))
    if (is.call(x) && as.character(x[[1]]) %in% c("<-", "=") &&
        as.character(x[[2]]) %in% c("edge_from_profiles", "calc_cutoffs_from_profiles", "constrain_cutoffs_readlengths")) eval(x, fns)

cols <- c(paste("5_UTR", 1:50, sep = "_"), paste("CDS", 1:99, sep = "_"), paste("3_UTR", 1:50, sep = "_"))

# stored results of one read length -> profile in the layout of RiboseQC_analysis (5_UTR_1..50, CDS_1..99, 3_UTR_1..50)
rebuild_profile <- function(res) {
    win <- as.matrix(res$profiles_start)
    prof <- matrix(0, nrow = nrow(win), ncol = length(cols), dimnames = list(rownames(win), cols))
    prof[, 51 + as.numeric(colnames(win))] <- win                       # position 0 is CDS_1, column 51
    # frame of the 5'ends over the CDS: large values on CDS_4..6, so that the counts on CDS_1 do not change it
    fr <- as.matrix(res$frames)
    if (!is.null(rownames(fr)) && !is.null(rownames(win)) && any(rownames(win) %in% rownames(fr))) {
        tx <- intersect(rownames(win), rownames(fr))
        prof[tx, c("CDS_4", "CDS_5", "CDS_6")] <- fr[tx, , drop = FALSE] * 1e9
    } else {
        prof[, c("CDS_4", "CDS_5", "CDS_6")] <- matrix(colMeans(fr) * 1e9, nrow = nrow(prof), ncol = 3, byrow = TRUE)
    }
    as.data.frame(prof, check.names = FALSE)
}

out <- NULL
for (f in args) {
    sample <- sub("_results_RiboseQC(_all)?$", "", basename(f))
    res_all <- get(load(f))
    for (comp in names(res_all$selection_cutoffs$analysis_frame_cutoff)) {
        stored <- res_all$selection_cutoffs$analysis_frame_cutoff[[comp]]
        rls <- names(stored)[names(stored) != "all"]
        if (length(rls) == 0) next
        new <- list()
        for (rl in rls) {
            res <- tryCatch(fns$calc_cutoffs_from_profiles(rebuild_profile(stored[[rl]]), length_max = as.numeric(rl), method = "edge"),
                            error = function(e) NULL)
            if (!is.null(res)) new[[rl]] <- res
        }
        joint <- fns$constrain_cutoffs_readlengths(new, max_diff = 1)
        choice <- as.data.frame(res_all$selection_cutoffs$results_choice[[comp]]$final_choice)
        for (rl in rls) {
            one <- function(x, field) if (is.null(x[[rl]][[field]])) NA else x[[rl]][[field]]
            out <- rbind(out, data.frame(sample = sample, compartment = comp, read_length = as.numeric(rl),
                reads_start = sum(stored[[rl]]$profiles_start),
                frame_preference = round(as.numeric(stored[[rl]]$frames_res[2]), 1),
                selected = rl %in% as.character(choice$read_length),
                cutoff_stored = stored[[rl]]$final_cutoff,
                cutoff_edge = one(new, "final_cutoff"),
                cutoff_edge_joint = one(joint, "final_cutoff"),
                constraint = one(joint, "cutoff_constraint"),
                chain_margin = round(one(joint, "cutoff_chain_margin"), 2)))
        }
    }
}
out$changed <- !(out$cutoff_stored == out$cutoff_edge_joint)
out <- out[order(out$sample, out$compartment, out$read_length), ]
write.table(out, sep = "\t", quote = FALSE, row.names = FALSE)
