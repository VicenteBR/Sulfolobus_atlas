###############################
## Reproducibility + Bias QC ##
## Full script with stats log##
###############################

suppressPackageStartupMessages({
  library(ggplot2)
  library(DESeq2)
  library(pheatmap)
  library(dplyr)
  library(Gviz)
  library(RColorBrewer)
  library(vsn)
})

## ---------- Config ----------
setwd("")

stat_file <- "analysis_stats.txt"

# init the log file
sink(stat_file); cat("Statistical Testing Log\n=======================\n\n"); sink()

log_stat <- function(section, description, result) {
  sink(stat_file, append = TRUE)
  cat("Section:", section, "\n")
  cat("Description:", description, "\n")
  if (is.matrix(result) || is.data.frame(result)) {
    print(result)
  } else {
    capture.output(print(result), file = stat_file, append = TRUE)
  }
  cat("\n---\n\n")
  sink()
}

safe_require <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Package '", pkg, "' not installed. Skipping features that need it.")
    return(FALSE)
  }
  TRUE
}

## ---------- TPM normalization ----------
tpm_norm <- function(df) {
  stopifnot("Length" %in% colnames(df))
  count_mat <- df[ , 6:ncol(df), drop = FALSE]
  gene_length <- df$Length
  tpm <- function(counts, gene_length) {
    rpk <- counts / (gene_length / 1000)
    scaling_factor <- sum(rpk) / 1e6
    rpk / scaling_factor
  }
  tpm_mat <- as.data.frame(sapply(count_mat, function(x) tpm(x, gene_length)))
  colnames(tpm_mat) <- colnames(count_mat)
  tpm_mat <- round(tpm_mat, 2) + 1
  tpm_mat <- cbind(GeneID = rownames(df), tpm_mat)
  tpm_mat
}

## ---------- Load data + build length map ----------
length_map <- NULL  # will hold GeneID -> Length across any table that has it

grab_lengths <- function(df) {
  if (!("Length" %in% colnames(df))) return(NULL)
  data.frame(GeneID = rownames(df), Length = df$Length, row.names = NULL, check.names = FALSE)
}

# -- NUDIX
table_in <- read.table("NUDIX_ko", sep = "\t", skip = 1, header = TRUE, row.names = 1)
if (is.null(length_map)) length_map <- grab_lengths(table_in)
tpm_df <- tpm_norm(table_in)
nudix_df <- tpm_df[ , grepl("WT1L|WT2L|GeneID", colnames(tpm_df)), drop = FALSE]
colnames(nudix_df) <- c("GeneID","Nudix_1","Nudix_2")

# -- Heat shock
table_in <- read.table("heat-shock-2024", sep = "\t", skip = 1, header = TRUE, row.names = 1)
if (is.null(length_map)) length_map <- grab_lengths(table_in) else length_map <- distinct(rbind(length_map, grab_lengths(table_in)))
tpm_df <- tpm_norm(table_in)
hs_df <- tpm_df[ , grepl("75|GeneID", colnames(tpm_df)), drop = FALSE]
colnames(hs_df) <- c("GeneID","HS_1","HS_2","HS_3","HS_4")

# -- Growth
table_in <- read.table("Growth-pH", sep = "\t", skip = 1, header = TRUE, row.names = 1)
if (is.null(length_map)) length_map <- grab_lengths(table_in) else length_map <- distinct(rbind(length_map, grab_lengths(table_in)))
tpm_df <- tpm_norm(table_in)
growth_df <- tpm_df[ , grepl("75_2.4_log|GeneID", colnames(tpm_df)), drop = FALSE]
colnames(growth_df) <- c("GeneID","Growth_1","Growth_2")

# -- Starvation
table_in <- read.table("Starvation", sep = "\t", skip = 1, header = TRUE, row.names = 1)
if (is.null(length_map)) length_map <- grab_lengths(table_in) else length_map <- distinct(rbind(length_map, grab_lengths(table_in)))
tpm_df <- tpm_norm(table_in)
starv_df <- tpm_df[ , grepl("t0_|GeneID", colnames(tpm_df)), drop = FALSE]
colnames(starv_df) <- c("GeneID","Starvation_1","Starvation__2","Starvation_3","Starvation_4","Starvation_5","Starvation_6")

## ---------- Merge all reference TPMs ----------
all_tpm <- Reduce(function(x,y) merge(x,y, by = "GeneID", all = TRUE),
                  list(nudix_df, hs_df, growth_df, starv_df))
all_tpm[is.na(all_tpm)] <- 0

## ---------- Detection sets (TPM > 1 in ANY replicate) ----------
detected <- lapply(list(
  Nudix      = grep("^Nudix_", colnames(all_tpm), value=TRUE),
  HS         = grep("^HS_",    colnames(all_tpm), value=TRUE),
  Growth     = grep("^Growth_",colnames(all_tpm), value=TRUE),
  Starvation = grep("^Starvation_",colnames(all_tpm), value=TRUE)
), function(cols) {
  all_tpm$GeneID[rowSums(all_tpm[, cols, drop=FALSE] > 1) > 0]
})
names(detected) <- c("Nudix","HS","Growth","Starvation")

## ---------- Replicate-level correlations (matrix + p-values) ----------
log_mat <- log2(as.matrix(all_tpm[ , -1, drop = FALSE]) + 1)
cors <- cor(log_mat, method = "pearson", use = "pairwise.complete.obs")
log_stat("Replicate-level Correlations",
         "Pearson correlation matrix across all replicates (log2 TPM+1).", cors)

# p-value matrix
cor_pvals <- matrix(NA_real_, ncol(cors), nrow(cors),
                    dimnames = list(colnames(cors), colnames(cors)))
for (i in seq_len(ncol(cors))) {
  for (j in seq_len(ncol(cors))) {
    ct <- suppressWarnings(cor.test(log_mat[,i], log_mat[,j], method="pearson"))
    cor_pvals[i,j] <- ct$p.value
  }
}
log_stat("Replicate-level Correlations",
         "Pairwise Pearson correlation p-values across all replicates.", cor_pvals)

# Heatmap (visual)
pheatmap::pheatmap(cors, main = "Replicate-level correlations (Pearson R)")

## ---------- Average replicates per experiment ----------
avg_profiles <- data.frame(
  GeneID     = all_tpm$GeneID,
  Nudix      = rowMeans(all_tpm[, grep("^Nudix_", colnames(all_tpm)), drop = FALSE]),
  HS         = rowMeans(all_tpm[, grep("^HS_",    colnames(all_tpm)), drop = FALSE]),
  Growth     = rowMeans(all_tpm[, grep("^Growth_",colnames(all_tpm)), drop = FALSE]),
  Starvation = rowMeans(all_tpm[, grep("^Starvation_",colnames(all_tpm)), drop = FALSE])
)
log_avg <- cbind(GeneID = avg_profiles$GeneID, log2(avg_profiles[,-1] + 1))

## ---------- Scatter with R + p (helper) ----------
scatter_with_R <- function(x, y, xlab, ylab, title, section) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  ct <- cor.test(x, y, method = "pearson")
  log_stat(section, paste0("Scatter (Pearson): ", title), ct)
  
  plot(x, y, pch=16, cex=0.6, xlab=xlab, ylab=ylab,
       main=title, cex.lab=1.2, cex.main=1.3)
  abline(0,1,col="red",lwd=2)
  r <- unname(ct$estimate); p <- ct$p.value
  p_txt <- ifelse(p < 2.2e-16, "< 2.2e-16", signif(p, 3))
  usr <- par("usr")
  text(usr[1] + 0.05*(usr[2]-usr[1]),
       usr[4] - 0.05*(usr[4]-usr[3]),
       labels=paste0("R = ", round(r, 3), "\nP ", p_txt),
       adj=c(0,1), col="blue", cex=1.2, font=2)
}
## Helper for consistent axis labels with subscripted 2
axis_label <- function(name) bquote(.(name)~log[2]*"(TPM+1)")

## ---------- Pairwise scatter plots (ALL GENES) + stats table ----------
exps <- colnames(log_avg)[-1]
pairs_all <- combn(exps, 2, simplify = FALSE)

# Save & set bigger fonts (axis titles, ticks, main title)
.par_old <- par(no.readonly = TRUE)
par(mfrow=c(2,3), mar=c(5,5,3.5,1),
    cex.lab=1.6,   # axis titles
    cex.axis=1.3,  # tick labels
    cex.main=1.5,  # plot title
    mgp=c(2.6, 0.8, 0))  # spacing: title/labels/ticks baseline

all_pair_tests <- do.call(rbind, lapply(pairs_all, function(p) {
  x <- log_avg[[p[1]]]; y <- log_avg[[p[2]]]
  scatter_with_R(x, y,
                 xlab = axis_label(p[1]),
                 ylab = axis_label(p[2]),
                 title = paste(p[1], "vs", p[2]),
                 section = "All genes: averaged profiles")
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  ct <- cor.test(x, y, method="pearson")
  data.frame(X=p[1], Y=p[2], PearsonR=unname(ct$estimate), Pvalue=ct$p.value)
}))
par(mfrow=c(1,1))
all_pair_tests$AdjP_BH <- p.adjust(all_pair_tests$Pvalue, method="BH")
log_stat("All genes: averaged profiles",
         "Pairwise Pearson tests (R, p, BH-adjusted) for log\u2082(TPM+1).", all_pair_tests)

## ---------- asRNA/ncRNA only ----------
idx_as_nc <- grepl("asRNA_|ncRNA_", log_avg$GeneID)
log_avg_as <- log_avg[idx_as_nc, , drop = FALSE]
if (nrow(log_avg_as) == 0) stop("No rows match 'asRNA_' or 'ncRNA_' in GeneID.")

pairs_as <- combn(exps, 2, simplify = FALSE)
par(mfrow=c(2,3), mar=c(5,5,3.5,1),
    cex.lab=1.6, cex.axis=1.3, cex.main=1.5, mgp=c(2.6,0.8,0))

as_pair_tests <- do.call(rbind, lapply(pairs_as, function(p) {
  x <- log_avg_as[[p[1]]]; y <- log_avg_as[[p[2]]]
  scatter_with_R(x, y,
                 xlab = axis_label(p[1]),
                 ylab = axis_label(p[2]),
                 title = paste("asRNA/ncRNA:", p[1], "vs", p[2]),
                 section = "asRNA/ncRNA: averaged profiles")
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  ct <- cor.test(x, y, method="pearson")
  data.frame(X=p[1], Y=p[2], PearsonR=unname(ct$estimate), Pvalue=ct$p.value)
}))
par(mfrow=c(1,1))

# Restore original graphics settings
par(.par_old)

as_pair_tests$AdjP_BH <- p.adjust(as_pair_tests$Pvalue, method="BH")
log_stat("asRNA/ncRNA: averaged profiles",
         "Pairwise Pearson tests (R, p, BH-adjusted) for log\u2082(TPM+1) in asRNA/ncRNA subset.",
         as_pair_tests)



## ---------- Detection overlap (Jaccard) + overlap significance ----------
# Universe = all asRNA/ncRNA GeneIDs present in at least one experiment
idx_as_nc_all <- grepl("asRNA_|ncRNA_", all_tpm$GeneID)
universe <- all_tpm$GeneID[idx_as_nc_all]

detected_as <- lapply(exps, function(exp){
  cols <- grep(paste0("^", exp, "_"), colnames(all_tpm), value=TRUE)
  sel  <- grepl("asRNA_|ncRNA_", all_tpm$GeneID)
  all_tpm$GeneID[sel & rowSums(all_tpm[sel, cols, drop=FALSE] > 1) > 0]
})
names(detected_as) <- exps

jaccard <- matrix(NA_real_, length(exps), length(exps),
                  dimnames=list(exps, exps))
overlap_p <- matrix(NA_real_, length(exps), length(exps),
                    dimnames=list(exps, exps))

for (i in seq_along(exps)) for (j in seq_along(exps)) {
  A <- detected_as[[i]]; B <- detected_as[[j]]
  uN <- length(universe); aN <- length(A); bN <- length(B)
  k  <- length(intersect(A,B))
  jaccard[i,j] <- if (length(unique(c(A,B))) == 0) NA_real_ else k / length(unique(c(A,B)))
  # Hypergeometric: P(overlap >= k)
  overlap_p[i,j] <- if (all(c(uN,aN,bN) > 0)) phyper(q = k-1, m = aN, n = uN - aN, k = bN, lower.tail = FALSE) else NA_real_
}
log_stat("Detection overlap (asRNA/ncRNA)",
         "Jaccard index matrix across experiments (TPM>1 in any replicate).", jaccard)
log_stat("Detection overlap (asRNA/ncRNA)",
         "Hypergeometric p-values for overlap enrichment across experiments.", overlap_p)

# Heatmap of Jaccard
image(1:length(exps), 1:length(exps), t(jaccard[nrow(jaccard):1,]),
      xaxt="n", yaxt="n", xlab="Experiment", ylab="Experiment",
      main="Jaccard index (asRNA/ncRNA detection)")
axis(1, at=1:length(exps), labels=exps, las=2)
axis(2, at=1:length(exps), labels=rev(exps), las=2)
for (i in seq_along(exps)) for (j in seq_along(exps)) {
  val <- jaccard[i,j]; if (!is.na(val)) text(j, length(exps)-i+1, sprintf("%.2f", val))
}

## ---------- MA + correlation-vs-abundance with bin-wise stats (asRNA/ncRNA) ----------
ma_and_binnedR <- function(x, y, xname, yname) {
  A <- 0.5 * (x + y)
  M <- y - x
  
  # MA scatter
  plot(A, M, pch=16, cex=0.5,
       xlab="A = average log2(TPM+1)", ylab="M = diff (Y - X)",
       main=paste("MA:", xname, "vs", yname))
  abline(h=0, col="red", lwd=2)
  
  # Bin-wise correlations
  nb <- 10
  br <- quantile(A, probs=seq(0,1,length.out=nb+1), na.rm=TRUE)
  r_by_bin <- mid <- p_by_bin <- rep(NA_real_, nb)
  for (i in 1:nb) {
    idx <- A >= br[i] & A < br[i+1]
    mid[i] <- mean(br[i:(i+1)])
    if (sum(idx) >= 5) {
      ct <- suppressWarnings(cor.test(x[idx], y[idx], method="pearson"))
      r_by_bin[i] <- unname(ct$estimate)
      p_by_bin[i] <- ct$p.value
    }
  }
  # correlation vs abundance plot
  plot(mid, r_by_bin, type="b", pch=16, ylim=c(0,1),
       xlab="A (bin centers, log2 TPM+1)", ylab="Pearson R",
       main=paste("Correlation vs abundance:", xname,"vs",yname))
  abline(h=0:1, col="grey85", lty=3)
  
  # Log stats
  out <- data.frame(A_center=round(mid,2), R=r_by_bin, Pvalue=p_by_bin)
  out$AdjP_BH <- p.adjust(out$Pvalue, method="BH")
  log_stat("asRNA/ncRNA: correlation vs abundance",
           paste0("Bin-wise Pearson tests for ", xname, " vs ", yname), out)
}

pairs_as <- combn(exps, 2, simplify=FALSE)
op <- par(mfrow=c(2,3), mar=c(4.5,4.5,3,1))
for (p in pairs_as) {
  x <- log_avg_as[[p[1]]]; y <- log_avg_as[[p[2]]]
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  ma_and_binnedR(x, y, p[1], p[2])
}
par(op)

## ---------- PCA (asRNA/ncRNA replicates) ----------
cols_repl <- grep("^(Nudix|HS|Growth|Starvation)_", colnames(all_tpm), value=TRUE)
log_mat_as <- log2(as.matrix(all_tpm[idx_as_nc_all, cols_repl, drop=FALSE]) + 1)
pca <- prcomp(t(log_mat_as), scale.=FALSE)
pc <- pca$x[,1:2]
expl <- round(100 * (pca$sdev^2)/sum(pca$sdev^2), 1)
labs <- colnames(all_tpm)[cols_repl]
groups <- sub("_.*","", labs)
cols <- as.integer(as.factor(groups))

plot(pc[,1], pc[,2], pch=16, cex=1,
     xlab=paste0("PC1 (", expl[1], "%)"),
     ylab=paste0("PC2 (", expl[2], "%)"),
     main="PCA: Reference replicates (asRNA/ncRNA)")
abline(h=0,v=0,col="grey85")
text(pc[,1], pc[,2], labels=labs, pos=3, cex=0.8)
legend("topright", legend=levels(as.factor(groups)),
       col=seq_along(levels(as.factor(groups))), pch=16, bty="n")

# Log PCA variance explained
pca_var <- data.frame(PC = paste0("PC", seq_along(pca$sdev)),
                      VarianceExplained = round(100*(pca$sdev^2)/sum(pca$sdev^2), 2))
log_stat("PCA (asRNA/ncRNA)",
         "Variance explained by principal components (replicates).", pca_var)

## ---------- LENGTH bias tests across ALL pairs ----------
length_map = data.frame(GeneID = all_tpm$GeneID, Length = table_in$Length)


# Ensure we have a usable length_map
if (!is.null(length_map) && "GeneID" %in% colnames(length_map) && "Length" %in% colnames(length_map)) {
  # Deduplicate length_map by first occurrence
  length_map <- length_map[!duplicated(length_map$GeneID), , drop = FALSE]
  
  # Build per-gene asRNA avg table + length
  avg_as <- avg_profiles[idx_as_nc, , drop=FALSE]
  as_len <- merge(avg_as[, c("GeneID", exps)], length_map, by="GeneID", all.x=TRUE)
  as_len$Length_kb <- as_len$Length / 1000
  
  # For |M| vs length stats we need log2 avg TPM+1 subset with GeneID
  log_as <- cbind(GeneID = log_avg_as$GeneID, log_avg_as[, exps, drop=FALSE])
  
  len_kw_tbl <- list()
  len_pw_list <- list()
  m_vs_len_tbl <- list()
  k <- 1
  
  op <- par(mfrow=c(2,3), mar=c(4.5,4.5,3,1))
  for (p in pairs_as) {
    Aname <- p[1]; Bname <- p[2]
    
    # Detection grouping (Both / A only / B only / None)
    Aset <- detected_as[[Aname]]
    Bset <- detected_as[[Bname]]
    as_len$group <- ifelse(as_len$GeneID %in% intersect(Aset,Bset), "Both",
                           ifelse(as_len$GeneID %in% Aset, paste0(Aname, " only"),
                                  ifelse(as_len$GeneID %in% Bset, paste0(Bname, " only"), "None")))
    
    # Boxplot
    boxplot(Length_kb ~ group, data=as_len,
            main=paste("Length by detection:", Aname, "vs", Bname),
            ylab="Length (kb)", las=2)
    
    # Stats: Kruskal–Wallis & pairwise Wilcoxon (BH)
    kw <- suppressWarnings(kruskal.test(Length_kb ~ group, data=as_len))
    pw <- suppressWarnings(pairwise.wilcox.test(as_len$Length_kb, as_len$group, p.adjust.method="BH"))
    len_kw_tbl[[k]] <- data.frame(X=Aname, Y=Bname, KruskalP=kw$p.value)
    len_pw_df <- as.data.frame(as.table(pw$p.value))
    names(len_pw_df) <- c("Group1","Group2","AdjP_BH")
    len_pw_df$X <- Aname; len_pw_df$Y <- Bname
    len_pw_list[[k]] <- len_pw_df
    
    # |M| vs Length (Spearman)
    x <- log_as[[Aname]]; y <- log_as[[Bname]]
    df_m <- merge(data.frame(GeneID = log_as$GeneID, M = abs(y - x)),
                  length_map, by="GeneID", all.x=TRUE)
    plot(df_m$Length/1000, df_m$M, pch=16, cex=0.5,
         xlab="Length (kb)", ylab="|M| = |Δ log2(TPM+1)|",
         main=paste("|M| vs Length:", Aname, "vs", Bname))
    lines(stats::lowess(df_m$Length/1000, df_m$M, f=0.4), lwd=2, col="red")
    
    sp <- suppressWarnings(cor.test(df_m$Length/1000, df_m$M, method="spearman", use="complete.obs"))
    m_vs_len_tbl[[k]] <- data.frame(X=Aname, Y=Bname,
                                    SpearmanRho=unname(sp$estimate),
                                    Pvalue=sp$p.value)
    k <- k + 1
  }
  par(op)
  
  len_kw_tbl <- do.call(rbind, len_kw_tbl)
  len_kw_tbl$AdjP_BH <- p.adjust(len_kw_tbl$KruskalP, "BH")
  log_stat("Length bias (asRNA/ncRNA)",
           "Kruskal–Wallis tests for length differences among detection groups (per pair).",
           len_kw_tbl)
  
  len_pw_tbl <- do.call(rbind, len_pw_list)
  log_stat("Length bias (asRNA/ncRNA)",
           "Pairwise Wilcoxon tests (BH-adjusted) for length among detection groups (per pair).",
           len_pw_tbl)
  
  m_vs_len_tbl <- do.call(rbind, m_vs_len_tbl)
  m_vs_len_tbl$AdjP_BH <- p.adjust(m_vs_len_tbl$Pvalue, "BH")
  log_stat("|M| vs Length (asRNA/ncRNA)",
           "Spearman tests for association between |Δ log2(TPM+1)| and length (per pair).",
           m_vs_len_tbl)
  
} else {
  log_stat("Length bias", "Length map not available; length-based tests skipped.", "N/A")
}

#########################
## End of analysis run ##
#########################

message("Done. All statistical outputs saved to: ", normalizePath(stat_file))
