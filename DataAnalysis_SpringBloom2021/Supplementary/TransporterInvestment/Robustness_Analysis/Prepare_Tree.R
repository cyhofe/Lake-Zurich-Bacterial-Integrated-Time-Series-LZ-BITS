############################################################
# PREPARE PHYLOGENETIC TREE FOR PGLS / PIC
#
# ROOTED =
#   Archaeon NZE17MAR2102_MAG_1102 defines the bacterial
#   root, then the archaeal tip is removed.
#
# CLEANED =
#   1. Remove the two invalid Chloroflexota candidate roots.
#   2. Remove the archaeal taxon AFTER using it for rooting.
#   3. Restore "_MAG_" -> ".MAG_" for Species_ID matching.
#
# NOT DONE HERE =
#   - no filtering to retained TPF / Bray-Curtis species
#   - no ultrametric conversion
#   - no branch-length modification
#   - no topology manipulation
############################################################


############################################################
# 0. SETUP
############################################################

rm(list = ls())
graphics.off()

setwd(
  "/SynologyDrive/PhD/Projects/SpringBloom2021/DataAnalysis_SpringBloom2021/Supplementary/TransporterInvestment/Robustness_Analysis/"
)

suppressPackageStartupMessages(library(ape))

taxonomy_file <- "perc-32-37-marked_labels.txt"

old_candidates <- c(
  "GCA_005881435_1",
  "GCA_035709505_1"
)

archaeal_outgroup <- "NZE17MAR2102_MAG_1102"

out_dir <- "Prepared_Tree"
dir.create(out_dir, showWarnings = FALSE)

clean_id <- function(x) sub("_MAG_", ".MAG_", x, fixed = TRUE)

section <- function(x) {
  cat("\n============================================================\n",
      x,
      "\n============================================================\n",
      sep = "")
}


############################################################
# 1. LOAD AND CHECK IQ-TREE TREE
############################################################

tree_files <- list.files(pattern = "\\.treefile$", full.names = TRUE)

if (length(tree_files) != 1) {
  stop("Expected exactly one *.treefile; found ", length(tree_files))
}

tree_file <- tree_files[1]
tree_raw <- read.tree(tree_file)

if (!inherits(tree_raw, "phylo")) stop("Tree was not read as a phylo object.")
if (is.null(tree_raw$edge.length)) stop("Tree has no branch lengths.")
if (any(!is.finite(tree_raw$edge.length))) stop("Non-finite branch lengths found.")
if (any(tree_raw$edge.length < 0)) stop("Negative branch lengths found.")
if (anyDuplicated(tree_raw$tip.label)) stop("Duplicated tree tip labels found.")

required <- c(old_candidates, archaeal_outgroup)

if (!all(required %in% tree_raw$tip.label)) {
  stop(
    "Missing required taxa: ",
    paste(setdiff(required, tree_raw$tip.label), collapse = ", ")
  )
}

section("RAW TREE QC")

cat("Input:          ", tree_file, "\n")
cat("Tips:           ", Ntip(tree_raw), "\n")
cat("Internal nodes: ", Nnode(tree_raw), "\n")
cat("Rooted:         ", is.rooted(tree_raw), "\n")
cat("Binary:         ", is.binary(tree_raw), "\n")
cat("Zero branches:  ", sum(tree_raw$edge.length == 0), "\n\n")

print(summary(tree_raw$edge.length))

write.tree(
  tree_raw,
  file.path(out_dir, "SpringBloom_original_unrooted.nwk")
)


############################################################
# 2. LOAD TAXONOMY AND EXTRACT PHYLUM
#
# Actual taxonomy format:
#
# NZE...MAG_346   d_Bacteria_p_Bacteroidota_c_Bacteroidia...
#
# Also supports standard GTDB p__Phylum; format.
############################################################

if (!file.exists(taxonomy_file)) {
  stop("Taxonomy file not found: ", taxonomy_file)
}

tax_lines <- trimws(readLines(taxonomy_file, warn = FALSE))
tax_lines <- tax_lines[nzchar(tax_lines)]

tax <- data.frame(
  Tree_ID = sub("^([^[:space:]]+).*$", "\\1", tax_lines),
  Taxonomy = sub("^[^[:space:]]+[[:space:]]+", "", tax_lines),
  stringsAsFactors = FALSE
)

if (anyDuplicated(tax$Tree_ID)) stop("Duplicated IDs in taxonomy file.")

extract_phylum <- function(x) {
  
  out <- rep(NA_character_, length(x))
  
  # Standard GTDB format: p__Bacteroidota;
  i <- grepl("p__", x, fixed = TRUE)
  out[i] <- sub("^.*p__([^;]+).*$", "\\1", x[i])
  
  # Your Dryad format: _p_Bacteroidota_c_
  j <- !i & grepl("_p_", x, fixed = TRUE) &
    grepl("_c_", x, fixed = TRUE)
  
  out[j] <- sub(
    "^.*?_p_(.*?)_c_.*$",
    "\\1",
    x[j],
    perl = TRUE
  )
  
  out
}

tax$Phylum <- extract_phylum(tax$Taxonomy)
tax$Species_ID <- clean_id(tax$Tree_ID)

section("ROOTING TAXA")

print(
  tax[
    tax$Tree_ID %in% required,
    c("Tree_ID", "Phylum", "Taxonomy")
  ],
  row.names = FALSE
)


############################################################
# 3. CREATE UNROOTED WORKING TREE
############################################################

tree_work <- tree_raw

if (is.rooted(tree_work)) {
  tree_work <- unroot(tree_work, keep.root.edge = FALSE)
}


############################################################
# 4. REMOVE INVALID CHLOROFLEXOTA CANDIDATE ROOTS
############################################################

tree_with_archaea <- drop.tip(
  tree_work,
  old_candidates
)

if (!archaeal_outgroup %in% tree_with_archaea$tip.label) {
  stop("Archaeal outgroup unexpectedly missing.")
}

section("REMOVE OLD CANDIDATE ROOTS")

cat("Starting tips: ", Ntip(tree_work), "\n")
cat("Removed:       ", paste(old_candidates, collapse = ", "), "\n")
cat("Remaining:     ", Ntip(tree_with_archaea), "\n")


############################################################
# 5. INSPECT ARCHAEAL OUTGROUP
############################################################

D <- cophenetic.phylo(tree_with_archaea)

archaeal_distances <- sort(
  D[
    archaeal_outgroup,
    setdiff(tree_with_archaea$tip.label, archaeal_outgroup)
  ]
)

n_neighbors <- min(20, length(archaeal_distances))

nearest_archaea <- data.frame(
  Rank = seq_len(n_neighbors),
  Neighbor = names(archaeal_distances)[seq_len(n_neighbors)],
  Distance = as.numeric(archaeal_distances[seq_len(n_neighbors)])
)

write.csv(
  nearest_archaea,
  file.path(out_dir, "Archaeal_outgroup_nearest_neighbors.csv"),
  row.names = FALSE
)

# Terminal branch lengths
tip_edges <- match(
  seq_len(Ntip(tree_with_archaea)),
  tree_with_archaea$edge[, 2]
)

terminal_lengths <- tree_with_archaea$edge.length[tip_edges]
names(terminal_lengths) <- tree_with_archaea$tip.label

# IMPORTANT: consistent variable name everywhere
archaeal_branch <- unname(
  terminal_lengths[archaeal_outgroup]
)

archaeal_percentile <- mean(
  terminal_lengths <= archaeal_branch
)

section("ARCHAEAL OUTGROUP QC")

cat("Terminal branch:           ", archaeal_branch, "\n")
cat("Terminal branch percentile:",
    round(100 * archaeal_percentile, 2), "%\n")
cat("Nearest bacterial taxon:   ", names(archaeal_distances)[1], "\n")
cat("Distance to nearest taxon: ", archaeal_distances[1], "\n\n")

print(head(nearest_archaea, 10))

if (archaeal_percentile >= 0.99) {
  warning(
    "The archaeal outgroup is among the longest terminal branches. ",
    "Inspect its placement carefully; rooting is still performed."
  )
}


############################################################
# 6. CREATE UNROOTED BACTERIAL COMPARISON TREE
#
# Archaeon is removed WITHOUT first establishing a root.
# This gives the comparison tree for the unrooted PDF.
############################################################

tree_unrooted <- drop.tip(
  tree_with_archaea,
  archaeal_outgroup,
  rooted = FALSE
)

if (is.rooted(tree_unrooted)) {
  tree_unrooted <- unroot(tree_unrooted, keep.root.edge = FALSE)
}

cat("\nUnrooted bacterial tree:\n")
cat("Tips:   ", Ntip(tree_unrooted), "\n")
cat("Rooted: ", is.rooted(tree_unrooted), "\n")


############################################################
# 7. ROOT USING THE ARCHAEAL OUTGROUP
#
# resolve.root = FALSE avoids deliberately adding an
# artificial zero-length root branch.
############################################################

tree_positioned <- root(
  tree_with_archaea,
  outgroup = archaeal_outgroup,
  resolve.root = FALSE,
  edgelabel = TRUE
)

# Save provenance tree containing the archaeal outgroup.
write.tree(
  tree_positioned,
  file.path(out_dir, "SpringBloom_archaeal_outgroup_positioned.nwk")
)


############################################################
# 8. REMOVE ARCHAEON AFTER ROOTING
#
# The archaeal TIP is removed.
# The ROOT POSITION it established is retained.
############################################################

tree_rooted <- drop.tip(
  tree_positioned,
  archaeal_outgroup
)

if (!is.rooted(tree_rooted)) {
  stop("Final bacterial tree is not recognized as rooted.")
}

if (!is.binary(tree_rooted)) {
  stop("Final rooted bacterial tree is not bifurcating.")
}


############################################################
# 9. CHECK ROOTING DID NOT CHANGE BACTERIAL DISTANCES
############################################################

common <- sort(tree_rooted$tip.label)

D_unrooted <- cophenetic.phylo(tree_unrooted)[common, common]
D_rooted   <- cophenetic.phylo(tree_rooted)[common, common]

max_distance_change <- max(abs(D_unrooted - D_rooted))

cat(
  "\nMaximum pairwise-distance change after rooting: ",
  max_distance_change,
  "\n",
  sep = ""
)

if (max_distance_change > 1e-8) {
  stop("Rooting unexpectedly changed bacterial patristic distances.")
}


############################################################
# 10. CLEAN Species_ID LABELS
#
# "Cleaned" means:
#
#   NZE05MAY2102_MAG_534
#          becomes
#   NZE05MAY2102.MAG_534
#
# No additional species filtering occurs here.
############################################################

tree_clean <- tree_rooted
tree_clean$tip.label <- clean_id(tree_clean$tip.label)

tree_unrooted_clean <- tree_unrooted
tree_unrooted_clean$tip.label <- clean_id(tree_unrooted_clean$tip.label)

if (anyDuplicated(tree_clean$tip.label)) {
  stop("Cleaning created duplicated Species_ID values.")
}

if (!setequal(tree_clean$tip.label, tree_unrooted_clean$tip.label)) {
  stop("Rooted and unrooted trees do not contain identical bacterial taxa.")
}


############################################################
# 11. MATCH PHYLUM TO TREE TIPS
#
# Phylum is plotting metadata ONLY.
# It is NOT added to analytical Newick tip names.
############################################################

phylum <- tax$Phylum[
  match(tree_clean$tip.label, tax$Species_ID)
]

n_missing_phylum <- sum(is.na(phylum))
phylum[is.na(phylum)] <- "Unknown_Phylum"

tree_phylum <- data.frame(
  Species_ID = tree_clean$tip.label,
  Phylum = phylum,
  stringsAsFactors = FALSE
)

write.csv(
  tree_phylum,
  file.path(out_dir, "SpringBloom_tree_phylum_mapping.csv"),
  row.names = FALSE
)

section("PHYLUM MATCHING")

cat("Tree tips:            ", Ntip(tree_clean), "\n")
cat("Matched phylum:       ", Ntip(tree_clean) - n_missing_phylum, "\n")
cat("Missing phylum:       ", n_missing_phylum, "\n\n")

print(sort(table(tree_phylum$Phylum), decreasing = TRUE))

if (n_missing_phylum > 0) {
  warning(n_missing_phylum, " tree tip(s) lack a matched phylum.")
}


############################################################
# 12. FINAL QC
############################################################

n_zero     <- sum(tree_clean$edge.length == 0)
n_negative <- sum(tree_clean$edge.length < 0)

expected_tips <- Ntip(tree_raw) - length(old_candidates) - 1

section("FINAL ROOTED + CLEANED TREE")

cat("Expected tips:     ", expected_tips, "\n")
cat("Observed tips:     ", Ntip(tree_clean), "\n")
cat("Internal nodes:    ", Nnode(tree_clean), "\n")
cat("Rooted:            ", is.rooted(tree_clean), "\n")
cat("Binary:            ", is.binary(tree_clean), "\n")
cat("Ultrametric:       ", is.ultrametric(tree_clean), "\n")
cat("Zero branches:     ", n_zero, "\n")
cat("Negative branches: ", n_negative, "\n\n")

print(summary(tree_clean$edge.length))

if (Ntip(tree_clean) != expected_tips) stop("Unexpected final number of tips.")
if (!is.rooted(tree_clean)) stop("Final analytical tree is not rooted.")
if (!is.binary(tree_clean)) stop("Final analytical tree is not binary.")
if (n_negative > 0) stop("Negative branches in final tree.")


############################################################
# 13. SAVE ANALYTICAL TREES
############################################################

final_tree_file <- file.path(
  out_dir,
  "SpringBloom_rooted_cleaned_for_analysis.nwk"
)

unrooted_tree_file <- file.path(
  out_dir,
  "SpringBloom_unrooted_cleaned_comparison.nwk"
)

write.tree(tree_clean, final_tree_file)
write.tree(tree_unrooted_clean, unrooted_tree_file)


############################################################
# 14. PREPARE PDF LABELS
#
# Analytical Newick:
#   NZE05MAY2102.MAG_534
#
# PDF only:
#   NZE05MAY2102.MAG_534 [Chloroflexota]
############################################################

plot_labels <- paste0(
  tree_phylum$Species_ID,
  " [",
  tree_phylum$Phylum,
  "]"
)

names(plot_labels) <- tree_phylum$Species_ID


############################################################
# 15. UNROOTED PDF
############################################################

plot_unrooted <- tree_unrooted_clean

plot_unrooted$tip.label <- unname(
  plot_labels[plot_unrooted$tip.label]
)

if (anyNA(plot_unrooted$tip.label)) {
  stop("Missing labels in unrooted plot.")
}

unrooted_pdf <- file.path(
  out_dir,
  "SpringBloom_UNROOTED_with_phylum.pdf"
)

pdf(
  unrooted_pdf,
  width = 32,
  height = 32,
  useDingbats = FALSE
)

par(mar = c(1, 1, 3, 1))

plot(
  plot_unrooted,
  type = "unrooted",
  use.edge.length = TRUE,
  cex = 0.26,
  lab4ut = "axial",
  no.margin = FALSE
)

title(
  "SpringBloom bacterial phylogeny - UNROOTED",
  cex.main = 1
)

add.scale.bar(cex = 0.7)

dev.off()


############################################################
# 16. ROOTED PDF
#
# ROOT is explicitly marked.
# Ladderization is only for visualization.
############################################################

plot_rooted <- ladderize(tree_clean, right = TRUE)

plot_rooted$tip.label <- unname(
  plot_labels[plot_rooted$tip.label]
)

if (anyNA(plot_rooted$tip.label)) {
  stop("Missing labels in rooted plot.")
}

root_node <- unique(
  setdiff(
    plot_rooted$edge[, 1],
    plot_rooted$edge[, 2]
  )
)

if (length(root_node) != 1) {
  stop("Could not identify exactly one root node.")
}

rooted_pdf <- file.path(
  out_dir,
  "SpringBloom_ROOTED_with_phylum.pdf"
)

pdf(
  rooted_pdf,
  width = 20,
  height = max(30, Ntip(plot_rooted) * 0.11),
  useDingbats = FALSE
)

par(mar = c(1, 2, 3, 1))

plot(
  plot_rooted,
  type = "phylogram",
  direction = "rightwards",
  use.edge.length = TRUE,
  show.tip.label = TRUE,
  cex = 0.36,
  label.offset = 0.003,
  no.margin = FALSE
)

# Graphical marker only; does NOT alter the tree.
nodelabels(
  text = "ROOT",
  node = root_node,
  frame = "rect",
  cex = 0.8
)

add.scale.bar(cex = 0.7)

title(
  paste0(
    "SpringBloom bacterial phylogeny - ROOTED\n",
    "Root established using archaeal outgroup ",
    clean_id(archaeal_outgroup),
    "; archaeal tip removed afterward"
  ),
  cex.main = 0.9
)

dev.off()


############################################################
# 17. SAVE QC TABLE
############################################################

qc <- data.frame(
  Metric = c(
    "Input_tree",
    "Original_tips",
    "Removed_Chloroflexota_candidates",
    "Archaeal_outgroup",
    "Archaeal_terminal_branch",
    "Archaeal_terminal_branch_percentile",
    "Archaeal_nearest_taxon",
    "Archaeal_nearest_taxon_distance",
    "Final_bacterial_tips",
    "Final_rooted",
    "Final_binary",
    "Final_ultrametric",
    "Final_zero_branches",
    "Final_negative_branches",
    "Maximum_distance_change_after_rooting",
    "Tips_with_known_phylum",
    "Tips_with_unknown_phylum"
  ),
  Value = c(
    basename(tree_file),
    Ntip(tree_raw),
    paste(old_candidates, collapse = ";"),
    archaeal_outgroup,
    archaeal_branch,
    archaeal_percentile,
    names(archaeal_distances)[1],
    archaeal_distances[1],
    Ntip(tree_clean),
    is.rooted(tree_clean),
    is.binary(tree_clean),
    is.ultrametric(tree_clean),
    n_zero,
    n_negative,
    max_distance_change,
    Ntip(tree_clean) - n_missing_phylum,
    n_missing_phylum
  ),
  stringsAsFactors = FALSE
)

write.csv(
  qc,
  file.path(out_dir, "SpringBloom_tree_QC.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "R_sessionInfo.txt")
)


############################################################
# 18. FINAL REPORT
############################################################

section("DONE")

cat("MASTER ROOTED ANALYSIS TREE:\n", final_tree_file, "\n\n", sep = "")
cat("UNROOTED COMPARISON TREE:\n", unrooted_tree_file, "\n\n", sep = "")
cat("ROOTED PDF:\n", rooted_pdf, "\n\n", sep = "")
cat("UNROOTED PDF:\n", unrooted_pdf, "\n\n", sep = "")

cat(
  "Final tree has NOT been filtered to retained TPF or ",
  "Bray-Curtis species.\n"
)