# CRC-ImmuneMap | Step 06: Tumor patient-level association and robustness
# INPUT: Step 05 saved per-patient scores; Step 04 optional sample denominators.
# Research question: SPP1_TAM vs FAP_CAF vs CD8_Dysfunction across tumor patients.
# PRIMARY: tumor patients with >=20 cells in EACH of the 3 compartments.
# Unit of replication is PATIENT, never individual cells.
# This is observational / exploratory, not causal inference.
# Use identical frozen candidate signatures for independent validation.

required <- c("dplyr", "tidyr", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

project_dir <- "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap"
if (!dir.exists(project_dir)) stop("Project directory not found: edit project_dir")
setwd(project_dir)
input_dir <- "results/tables/programs"
out_dir <- "results/tables/associations"
fig_dir <- "results/figures/associations"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

infile <- file.path(input_dir, "candidate_program_scores_by_patient_condition.csv")
manifest_file <- file.path(input_dir, "candidate_signature_manifest_v1.csv")
coverage_file <- file.path(input_dir, "signature_gene_coverage.csv")
for (f in c(infile, manifest_file, coverage_file)) {
  if (!file.exists(f)) stop("Required Step 05 file missing: ", f)
}
sc <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)
manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)
coverage <- read.csv(coverage_file, stringsAsFactors = FALSE)
needed <- c("patient_id", "condition", "program", "lineage", "n_cells",
            "mean_score", "median_score")
if (!all(needed %in% names(sc))) stop("Missing columns: ", paste(setdiff(needed,names(sc)),collapse=", "))
programs <- c("SPP1_TAM", "FAP_CAF", "CD8_Dysfunction", "CD8_Cytotoxicity")
if (!setequal(unique(manifest$program), programs)) stop("Unexpected candidate signatures in manifest")
if (nrow(manifest) != 29 || anyDuplicated(manifest[c("program","gene")])) {
  stop("Unexpected candidate manifest structure; inspect before proceeding")
}
if (!all(coverage$present) || nrow(coverage) != nrow(manifest)) {
  stop("Signature gene coverage not complete; stop rather than silently altering signatures")
}
if (anyDuplicated(sc[c("patient_id","condition","program")])) {
  stop("Duplicate patient-condition-program rows; investigate sample mapping")
}
if (!all(programs %in% unique(sc$program))) stop("One or more program scores missing")
if (any(!sc$condition %in% c("Tumor", "Normal"))) stop("Unexpected condition names")
if (anyNA(sc$patient_id) || anyNA(sc$n_cells) || any(sc$n_cells < 0)) {
  stop("Invalid identifiers / cell counts")
}
if (any(sc$n_cells == 0 & !is.na(sc$mean_score))) {
  stop("Absent compartments must have NA scores, not zero")
}
if (any(sc$n_cells > 0 & (!is.finite(sc$mean_score) | !is.finite(sc$median_score)))) {
  stop("Non-finite scores in nonempty compartments")
}

primary_pairs <- data.frame(
  comparison = c("SPP1_TAM__FAP_CAF", "SPP1_TAM__CD8_Dysfunction",
                 "FAP_CAF__CD8_Dysfunction"),
  x = c("SPP1_TAM", "SPP1_TAM", "FAP_CAF"),
  y = c("FAP_CAF", "CD8_Dysfunction", "CD8_Dysfunction"),
  stringsAsFactors = FALSE
)
# Rank correlation can be estimated on a handful of patients, but claims should
# remain exploratory. Monte-Carlo two-sided permutation p-values avoid assuming
# conventional asymptotic p-values are precise at n ~ 11.
set.seed(20261014)
B <- 19999L
cor_stats <- function(df, x, y, B = 0L) {
  keep <- is.finite(df[[x]]) & is.finite(df[[y]])
  u <- df[keep, , drop=FALSE]
  n <- nrow(u)
  if (n < 5L || length(unique(u[[x]])) < 3L || length(unique(u[[y]])) < 3L) {
    return(data.frame(n_patients=n, rho=NA_real_, p_asymptotic=NA_real_,
                      p_permutation=NA_real_))
  }
  rx <- rank(u[[x]], ties.method="average")
  ry <- rank(u[[y]], ties.method="average")
  rho <- unname(cor(rx, ry, method="pearson"))
  approx_p <- suppressWarnings(cor.test(u[[x]], u[[y]],
                            method="spearman", exact=FALSE)$p.value)
  perm_p <- NA_real_
  if (B > 0) {
    permuted <- replicate(B, cor(rx, sample(ry), method="pearson"))
    perm_p <- (1 + sum(abs(permuted) >= abs(rho) - 1e-12)) / (B + 1)
  }
  data.frame(n_patients=n, rho=rho, p_asymptotic=approx_p,
             p_permutation=perm_p)
}

# Convert input scores to wide patient-row layout; retain counts for thresholds.
wide <- sc %>%
  select(patient_id, condition, program, n_cells, mean_score, median_score) %>%
  pivot_wider(names_from=program, values_from=c(n_cells, mean_score, median_score),
              names_sep="__") %>%
  arrange(condition, patient_id)
tumor <- wide %>% filter(condition == "Tumor")
if (nrow(tumor) != 23L) {
  warning("Expected 23 tumor patients in previous audit; observed ", nrow(tumor),
          ". Results continue but verify input and cohort identity.")
}

eligible <- function(data, min_cells = 20L, metric = "mean_score", complete_three=TRUE) {
  if (complete_three) {
    for (p in c("SPP1_TAM", "FAP_CAF", "CD8_Dysfunction")) {
      data <- data %>% filter(.data[[paste0("n_cells__",p)]] >= min_cells)
    }
  }
  out <- data %>% select(patient_id, condition)
  for (p in programs) {
    count <- data[[paste0("n_cells__",p)]]
    v <- data[[paste0(metric,"__",p)]]
    out[[p]] <- ifelse(count >= min_cells, v, NA_real_)
  }
  out
}

# Main analysis: same 11 patients for all three primary comparisons.
primary <- eligible(tumor, 20L, "mean_score", complete_three=TRUE)
if (nrow(primary) < 5) stop("Fewer than 5 complete tumor patients; primary correlations not interpretable")
write.csv(primary, file.path(out_dir,"primary_complete_case_patients_ge20.csv"),row.names=FALSE)

primary_results <- bind_rows(lapply(seq_len(nrow(primary_pairs)), function(i) {
  p <- primary_pairs[i,]
  stats <- cor_stats(primary, p$x, p$y, B=B)
  cbind(p, stats)
})) %>% mutate(p_BH_permutation=p.adjust(p_permutation, method="BH"),
              p_BH_asymptotic=p.adjust(p_asymptotic, method="BH"),
              analysis="Primary | complete 3-program tumor cohort | mean scores | >=20 cells")
write.csv(primary_results, file.path(out_dir,"primary_three_correlations_ge20.csv"),row.names=FALSE)

# Extra coverage: pairwise available patients, not the primary complete-case cohort.
# Pairwise sample sets are different for each comparison; report separately.
pairwise_available <- eligible(tumor,20L,"mean_score",complete_three=FALSE)
pairwise_exploratory <- bind_rows(lapply(seq_len(nrow(primary_pairs)), function(i) {
  p <- primary_pairs[i,]
  cbind(p,cor_stats(pairwise_available,p$x,p$y,B=0L))
})) %>% mutate(analysis="Exploratory | pairwise available tumor patients | >=20 cells")
write.csv(pairwise_exploratory,
          file.path(out_dir,"exploratory_pairwise_available_ge20.csv"),row.names=FALSE)

# Thresholds / median aggregation sensitivity, still with the SAME 3-program
# complete-case requirement. Each row explicitly reports its changing sample size.
sensitivity <- bind_rows(lapply(c(10L,20L,30L,50L),function(threshold) {
  bind_rows(lapply(c("mean_score", "median_score"),function(metric) {
    dd <- eligible(tumor, threshold, metric, complete_three=TRUE)
    bind_rows(lapply(seq_len(nrow(primary_pairs)),function(i) {
      p <- primary_pairs[i,]
      cbind(data.frame(min_cells=threshold,score_summary=metric),p,
            cor_stats(dd,p$x,p$y,B=0L))
    }))
  }))
}))
write.csv(sensitivity,file.path(out_dir,"threshold_and_median_sensitivity.csv"),row.names=FALSE)

# Leave-one-patient-out robustness: only 11-patient primary cohort.
leave_one_out <- bind_rows(lapply(seq_len(nrow(primary)),function(i) {
  dd <- primary[-i,,drop=FALSE]
  bind_rows(lapply(seq_len(nrow(primary_pairs)),function(j) {
    p <- primary_pairs[j,]
    cbind(data.frame(excluded_patient=primary$patient_id[i]),p,
          cor_stats(dd,p$x,p$y,B=0L))
  }))
}))
write.csv(leave_one_out,file.path(out_dir,"leave_one_patient_out_ge20.csv"),row.names=FALSE)
loo_summary <- leave_one_out %>% group_by(comparison) %>% summarise(
  n_omissions=n(), min_rho=if(all(is.na(rho))) NA_real_ else min(rho,na.rm=TRUE),
  max_rho=if(all(is.na(rho))) NA_real_ else max(rho,na.rm=TRUE),
  pct_positive=round(100*mean(rho>0,na.rm=TRUE),1),.groups="drop")
write.csv(loo_summary,file.path(out_dir,"leave_one_out_summary.csv"),row.names=FALSE)

# Orthogonal descriptive comparison: CD8 cytotoxic and dysfunction share a
# compartment and may co-occur. This does NOT imply exhausted cell function.
cd8_exploratory <- cor_stats(pairwise_available,"CD8_Dysfunction","CD8_Cytotoxicity",B=0L)
write.csv(cd8_exploratory,file.path(out_dir,"cd8_dysfunction_cytotoxicity_exploratory.csv"),row.names=FALSE)

# Tumor-vs-Normal within-patient descriptive comparisons only. Patient identity
# matching in GEO must be confirmed as genuinely paired clinically before inferential tests.
paired_descriptive <- bind_rows(lapply(programs,function(p) {
  tmp <- sc %>% filter(program == p, n_cells >= 20L) %>%
    select(patient_id,condition,n_cells,mean_score) %>%
    pivot_wider(names_from=condition,values_from=c(n_cells,mean_score),names_sep="__") %>%
    filter(is.finite(mean_score__Tumor) & is.finite(mean_score__Normal)) %>%
    mutate(program=p,delta_tumor_minus_normal=mean_score__Tumor-mean_score__Normal)
  tmp
})) %>% select(program,everything())
write.csv(paired_descriptive,file.path(out_dir,"potential_paired_tumor_normal_descriptive.csv"),row.names=FALSE)
paired_summary <- paired_descriptive %>% group_by(program) %>% summarise(
  n_matched_patients=n(), median_delta=median(delta_tumor_minus_normal),
  .groups="drop")
write.csv(paired_summary,file.path(out_dir,"potential_paired_program_coverage.csv"),row.names=FALSE)

# Chart 1: scatterplots - strictly within primary three-program complete cohort.
scatter_long <- bind_rows(lapply(seq_len(nrow(primary_pairs)),function(i) {
  p <- primary_pairs[i,]
  data.frame(patient_id=primary$patient_id,comparison=p$comparison,
             x_program=p$x, y_program=p$y,
             x=primary[[p$x]],y=primary[[p$y]])
}))
label_df <- primary_results %>% transmute(comparison,
 label=paste0("n=",n_patients, " | Spearman rho=",sprintf("%.2f",rho),
              " | permutation p=",formatC(p_permutation,format="f",digits=3),
              " | BH q=",formatC(p_BH_permutation,format="f",digits=3)))
p_scatter <- ggplot(scatter_long,aes(x,y)) +
  geom_point(size=2.8,alpha=.85) +
  geom_smooth(method="lm",se=FALSE,linewidth=.45,linetype="dashed",color="grey50") +
  facet_wrap(~comparison,scales="free",nrow=1) +
  geom_text(aes(label=patient_id),check_overlap=TRUE,vjust=-.75,size=2.5) +
  geom_text(data=label_df,aes(x=-Inf,y=Inf,label=label),
            inherit.aes=FALSE,hjust=-.02,vjust=1.1,size=2.8) +
  labs(title="Discovery: exploratory tumor patient-level program associations",
       subtitle="Same complete-case patients for all 3 comparisons; >=20 cells/compartment; dashed lines are visual guides",
       x="Program mean UCell score (x)",y="Program mean UCell score (y)") +
  theme_bw(base_size=10) + theme(strip.text=element_text(size=9),plot.title=element_text(face="bold"))
ggsave(file.path(fig_dir,"primary_three_associations_ge20.png"),p_scatter,
       width=17,height=5.3,dpi=300)

# Chart 2: how the estimated correlations vary with coverage cutoffs.
p_sensitivity <- sensitivity %>% filter(score_summary=="mean_score",!is.na(rho)) %>%
  ggplot(aes(x=min_cells,y=rho,group=comparison)) +
  geom_hline(yintercept=0,linetype="dotted",color="grey60") +
  geom_line() + geom_point(aes(size=n_patients)) + facet_wrap(~comparison,nrow=1) +
  scale_x_continuous(breaks=c(10,20,30,50)) + ylim(-1,1) +
  labs(title="Sensitivity to minimum cells per compartment",
       subtitle="Complete 3-program tumor patients; n changes with the cutoff",
       x="Minimum cells in each lineage",y="Spearman rho",size="Patients") +
  theme_bw(base_size=10)
ggsave(file.path(fig_dir,"threshold_sensitivity.png"),p_sensitivity,
       width=13,height=4.7,dpi=300)

cat("\n========== STEP 06 SUMMARY ==========\n")
cat("Tumor patients total:",nrow(tumor),"\n")
cat("Primary complete-case patients >=20:",nrow(primary),"\n")
cat("Primary cohort:",paste(primary$patient_id,collapse=", "),"\n")
cat("\nPRIMARY CORRELATIONS (permutation p, BH across 3 hypotheses):\n")
print(primary_results[,c("comparison","n_patients","rho","p_permutation","p_BH_permutation")],row.names=FALSE)
cat("\nPAIRWISE AVAILABLE (different patients per pair; exploratory only):\n")
print(pairwise_exploratory[,c("comparison","n_patients","rho","p_asymptotic")],row.names=FALSE)
cat("\nLEAVE-ONE-OUT RANGE:\n")
print(loo_summary,row.names=FALSE)
cat("\nPOTENTIALLY MATCHED TUMOR/NORMAL (>=20 cells in both):\n")
print(paired_summary,row.names=FALSE)
cat("\nCD8 dysfunction vs cytotoxicity (exploratory):\n")
print(cd8_exploratory,row.names=FALSE)
cat("\nSaved tables:",normalizePath(out_dir,winslash="/"),"\n")
cat("Saved figures:",normalizePath(fig_dir,winslash="/"),"\n")
cat("NOTE: These are exploratory patient-level associations, not proof of causal signaling.\n")
