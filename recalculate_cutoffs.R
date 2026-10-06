#!/usr/bin/env Rscript
# Compare the P-site cutoffs stored in existing RiboseQC results with the cutoffs the "onset" method gives.
#
#   Rscript recalculate_cutoffs.R <sample>_results_RiboseQC [more result files ...] > cutoffs_old_vs_new.tsv
#
# Nothing is rerun on the BAM files and no result file is changed. For every read length the stored 5'end profiles
# around start codons (25 nt, one row per transcript), the stored frame of each transcript and the stored aggregate
# 5'end profile are put back into the layout calc_cutoffs_from_profiles() expects, and the functions are taken from
# source/RiboseQC/R/riboseqc.R next to this script, so the installed package does not have to be rebuilt first.
#
# The per-transcript profiles stored in the results cover the positions -24..0 only. The aggregate profile
# (profiles_fivepr$five_prime_subcodon) covers the whole metagene and supplies position -25 and the first codons
# of the CDS, on which the onset tests look. In results without that aggregate profile (older or trimmed files),
# the persistence test of candidates -3 and 0 cannot be done and passes; the output says so in `note`.
#
# Columns: sample, compartment, read_length, reads_start (5'ends in the 25 nt), frame_preference, selected (used for
# P-sites in the stored run), cutoff_stored, cutoff_edge (0.99.2 method, read length on its own), cutoff_onset
# (onset method, read length on its own), onset_reads (5'ends at the onset, 0 if none found), cutoff_source
# (onset / edge_meta / vote), cutoff_joint (onset method, read lengths together), constraint (chain / anchor /
# imputed / conflict), chain_margin, changed (stored != joint), note.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("give one or more *_results_RiboseQC files")
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
fns <- new.env()
for (x in parse(file.path(here, "source", "RiboseQC", "R", "riboseqc.R"), keep.source = FALSE))
    if (is.call(x) && as.character(x[[1]]) %in% c("<-", "=") &&
        as.character(x[[2]]) %in% c("edge_from_profiles", "onset_from_profile", "calc_cutoffs_from_profiles",
                                    "constrain_cutoffs_readlengths")) eval(x, fns)

cols <- c(paste("5_UTR", 1:50, sep = "_"), paste("CDS", 1:99, sep = "_"), paste("3_UTR", 1:50, sep = "_"))

# stored results of one read length (+ the stored aggregate profile, or NULL) -> profile in the layout of
# RiboseQC_analysis (5_UTR_1..50, CDS_1..99, 3_UTR_1..50), one row per transcript plus one row with the rest of the aggregate
rebuild_profile <- function(res, agg) {
    win <- as.matrix(res$profiles_start)
    wincols <- 51 + as.numeric(colnames(win))                       # position 0 is CDS_1, column 51
    prof <- matrix(0, nrow = nrow(win), ncol = length(cols), dimnames = list(rownames(win), cols))
    prof[, wincols] <- win
    # frame of the 5'ends over the CDS, per transcript, weighted by the transcript's reads in the window
    # (on CDS_31..33: away from the first codons, on which the onset tests look)
    fr <- as.matrix(res$frames)
    frcols <- c("CDS_31", "CDS_32", "CDS_33")
    if (!is.null(rownames(fr)) && !is.null(rownames(win)) && any(rownames(win) %in% rownames(fr))) {
        tx <- intersect(rownames(win), rownames(fr))
        prof[tx, frcols] <- fr[tx, , drop = FALSE] * rowSums(win[tx, , drop = FALSE])
    } else {
        prof[, frcols] <- matrix(colMeans(fr), nrow = nrow(prof), ncol = 3, byrow = TRUE) * rowSums(win)
    }
    if (!is.null(agg)) {
        rest <- rep(0, length(cols)); names(rest) <- cols
        a <- as.numeric(agg[cols]); a[is.na(a)] <- 0; rest[] <- a
        rest[c(wincols, match(frcols, cols))] <- 0                   # already in the rows above
        prof <- rbind(prof, aggregate_rest = rest)
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
        aggs <- res_all$profiles_fivepr$five_prime_subcodon[[comp]]
        new <- list(); edge <- list(); note <- list()
        for (rl in rls) {
            agg <- NULL; note[[rl]] <- ""
            if (!is.null(aggs[[rl]]) && is.null(dim(aggs[[rl]]))) agg <- aggs[[rl]]
            if (is.null(agg)) note[[rl]] <- "no aggregate profile: persistence of -3 and 0 not tested"
            prof <- rebuild_profile(stored[[rl]], agg)
            res <- tryCatch(fns$calc_cutoffs_from_profiles(prof, length_max = as.numeric(rl), method = "onset"),
                            error = function(e) NULL)
            if (!is.null(res)) new[[rl]] <- res
            rese <- tryCatch(fns$calc_cutoffs_from_profiles(prof, length_max = as.numeric(rl), method = "edge"),
                             error = function(e) NULL)
            if (!is.null(rese)) edge[[rl]] <- rese
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
                cutoff_edge = one(edge, "final_cutoff"),
                cutoff_onset = one(new, "final_cutoff"),
                onset_reads = one(new, "onset_reads"),
                cutoff_source = one(new, "cutoff_source"),
                cutoff_joint = one(joint, "final_cutoff"),
                constraint = one(joint, "cutoff_constraint"),
                chain_margin = round(one(joint, "cutoff_chain_margin"), 2),
                note = note[[rl]]))
        }
    }
}
out$changed <- !(out$cutoff_stored == out$cutoff_joint)
out <- out[order(out$sample, out$compartment, out$read_length), ]
write.table(out, sep = "\t", quote = FALSE, row.names = FALSE)
