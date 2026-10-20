# 08_anes_pairs_bounds.R
#
# Year-by-item-pair Pearson correlations for the ANES issue items used by
# Baldassarri & Gelman (2008, AJS) and Kozlowski & Murphy (2021, SSR), plus
# the attainable bounds r_min / r_max for every pair-year.
#
# Item selection, missing-value handling and recodes follow Kozlowski &
# Murphy's makeBaldGelData.R exactly (git submodule at
# replication/reproduce_partisans_ssr). Their weighted correlations use
# WGCNA::cor, which weights means by w but cross-products by w^2; we compute
# that version (r_wgcna) only to check the replication, and use the standard
# design-weighted Pearson r (r_w) for the analysis, because the bounds must be
# computed from the same weighted marginals as the correlation they bound.
#
# Not yet wired into 00_run_pipeline.R: needs the ANES Cumulative Data File,
# which ANES terms do not allow us to redistribute (data/ANES/ is gitignored;
# download from electionstudies.org).
#
# Reads:  data/ANES/csv/anes_timeseries_cdf_csv_20260205.csv
# Writes: output/data/anes_pairs_bounds.rds

suppressMessages(library(data.table))

cdf_path <- "data/ANES/csv/anes_timeseries_cdf_csv_20260205.csv"

# ---------------------------------------------------------------------------
# Variables (from makeBaldGelData.R)
# ---------------------------------------------------------------------------
ovars_core   <- c("vcf0004", "vcf0803", "vcf0301")
names_core   <- c("year", "ideoid", "partyid")
ovars_design <- c("vcf0017", "vcf0009x", "vcf0009y", "vcf0009z")
names_design <- c("mode", "faceWT", "webWT", "pooledWT")

ovars_econ <- c("vcf0839", "vcf0886", "vcf0809", "vcf0806", "vcf0887", "vcf9050", "vcf0893",
                "vcf0894", "vcf9046", "vcf0890", "vcf9049", "vcf0889", "vcf0891", "vcf9047")
names_econ <- c("fedspend_govserv", "fedspend_poor", "gov_jobguar", "healthinsur",
                "fedspend_childcare", "fedspend_asstblacks", "fedspend_homeless",
                "fedspend_welfare", "fedspend_foodstamp", "fedspend_schools",
                "fedspend_socsec", "fedspend_aids", "fedspend_college", "fedspend_envi")
ovars_civr <- c("vcf9016", "vcf0867a", "vcf9018", "vcf9042", "vcf9037", "vcf0830", "vcf9014",
                "vcf9040", "vcf9013", "vcf9017", "vcf0817", "vcf0816", "vcf9039", "vcf9041",
                "vcf0814", "vcf9015", "vcf0813")
names_civr <- c("more_equal", "affirmact", "equal_treat", "blacks_deservemore", "blacks_fair",
                "aid2blacks", "rights2much", "blacks_specfav", "equal_opp", "worryLessEqual",
                "school_bus", "school_integ", "hard4blacks", "blacks_shldTry", "civrights2fast",
                "chancesNEQ", "blacks_poschange")
ovars_moral <- c("vcf0877", "vcf0876a", "vcf0852", "vcf0854", "vcf0878", "vcf0853", "vcf0851",
                 "vcf0834", "vcf0837", "vcf0838", "vcf9043")
names_moral <- c("gayer_military", "protect_homo", "adjust_morals", "tolerate_difvals",
                 "gays_adoptkids", "tradvals", "newLstylesBad", "women_equals", "abortion1",
                 "abortion2", "school_prayer")
ovars_frgn <- c("vcf0843", "vcf0811", "vcf0841", "vcf0892", "vcf0888", "vcf9048")
names_frgn <- c("defense_spend", "urban_unrest", "ruskies", "fedspend_frgnaid",
                "fedspend_crime", "fedspend_space")

all_ovars <- c(ovars_core, ovars_design, ovars_econ, ovars_civr, ovars_moral, ovars_frgn)
all_names <- c(names_core, names_design, names_econ, names_civr, names_moral, names_frgn)

econ_issues <- c("fedspend_govservRC", "fedspend_poor", "gov_jobguar", "healthinsur",
                 "fedspend_childcare", "fedspend_asstblacks", "fedspend_homeless",
                 "fedspend_welfare", "fedspend_foodstamp", "fedspend_schools",
                 "fedspend_socsec", "fedspend_aids", "fedspend_college", "fedspend_envi")
civr_issues <- c("more_equalRC", "affirmact", "equal_treat", "blacks_deservemore", "blacks_fair",
                 "aid2blacks", "rights2muchRC", "blacks_specfavRC", "equal_opp",
                 "worryLessEqualRC", "school_bus", "school_integ", "hard4blacks",
                 "blacks_shldTryRC", "civrights2fast", "chancesNEQ", "blacks_poschange")
moral_issues <- c("gayer_military", "protect_homo", "adjust_morals", "tolerate_difvals",
                  "gays_adoptkids", "tradvalsRC", "newLstylesBadRC", "women_equals",
                  "abortion4RC", "school_prayer")
frgn_issues <- c("defense_spend", "urban_unrest", "ruskies", "fedspend_frgnaid",
                 "fedspend_crimeRC", "fedspend_space")
issue_vars <- c(econ_issues, civr_issues, moral_issues, frgn_issues)
domain_of <- c(setNames(rep("econ", length(econ_issues)), econ_issues),
               setNames(rep("civr", length(civr_issues)), civr_issues),
               setNames(rep("moral", length(moral_issues)), moral_issues),
               setNames(rep("frgn", length(frgn_issues)), frgn_issues),
               partyid = "partyid", ideoid = "ideoid")

# ---------------------------------------------------------------------------
# Load and recode (as in makeBaldGelData.R)
# ---------------------------------------------------------------------------
cat("Reading ANES CDF...\n")
header <- names(fread(cdf_path, nrows = 0))
d <- fread(cdf_path, select = header[tolower(header) %in% all_ovars], showProgress = FALSE)
setnames(d, tolower(names(d)))
d <- d[vcf0004 >= 1972, ..all_ovars]
setnames(d, all_ovars, all_names)

for (v in setdiff(all_names, c("year", names_design))) {
  set(d, j = v, value = ifelse(d[[v]] %in% 1:7, d[[v]], NA_integer_))
}
d[, affirmact := ifelse(affirmact %in% 1:5, affirmact, NA_integer_)]
for (v in c("fedspend_asstblacks", "fedspend_foodstamp", "fedspend_socsec", "fedspend_space")) {
  set(d, j = v, value = ifelse(d[[v]] %in% 1:3, d[[v]], NA_integer_))
}
d[, `:=`(fedspend_govservRC = -fedspend_govserv + 8,
         more_equalRC       = -more_equal + 6,
         rights2muchRC      = -rights2much + 6,
         blacks_specfavRC   = -blacks_specfav + 6,
         worryLessEqualRC   = -worryLessEqual + 6,
         blacks_shldTryRC   = -blacks_shldTry + 6,
         tradvalsRC         = -tradvals + 6,
         newLstylesBadRC    = -newLstylesBad + 6,
         abortion1RC        = -abortion1 + 5,
         abortion2RC        = -abortion2 + 5,
         fedspend_crimeRC   = -fedspend_crime + 4)]
d[, abortion4RC := ifelse(is.na(abortion1RC), abortion2RC, abortion1RC)]

# Face-to-face respondents only, with face-to-face weights (as K&M)
d <- d[mode == 0 & !is.na(faceWT) & faceWT > 0]

# ---------------------------------------------------------------------------
# Correlations and bounds
# ---------------------------------------------------------------------------
# Exact bounds from two marginal distributions (scores ascending, proportions).
# r_max: comonotone (quantile) coupling, merging both items' cumulative
# cut-points. r_min(X, Y) = -r_max(X, -Y): the countermonotone coupling is the
# comonotone coupling of X with Y reversed.
rmax_from_props <- function(sx, px, sy, py) {
  cx <- cumsum(px); cy <- cumsum(py)
  cx[length(cx)] <- 1; cy[length(cy)] <- 1
  br <- sort(unique(c(0, cx, cy)))
  mid <- (head(br, -1) + tail(br, -1)) / 2
  w <- diff(br)
  qx <- sx[pmin(findInterval(mid, cx) + 1, length(sx))]
  qy <- sy[pmin(findInterval(mid, cy) + 1, length(sy))]
  mx <- sum(px * sx); my <- sum(py * sy)
  sdx <- sqrt(sum(px * sx^2) - mx^2); sdy <- sqrt(sum(py * sy^2) - my^2)
  (sum(w * qx * qy) - mx * my) / (sdx * sdy)
}
bounds_from_props <- function(sx, px, sy, py) {
  c(r_min = -rmax_from_props(sx, px, -rev(sy), rev(py)),
    r_max =  rmax_from_props(sx, px, sy, py))
}

pair_stats <- function(x, y, w) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]; w <- w[ok]
  if (length(x) < 30 || length(unique(x)) < 2 || length(unique(y)) < 2) return(NULL)
  wn <- w / sum(w)
  mx <- sum(wn * x); my <- sum(wn * y)
  r_w <- sum(wn * (x - mx) * (y - my)) / sqrt(sum(wn * (x - mx)^2) * sum(wn * (y - my)^2))
  # WGCNA::cor weighted version (means weighted by w, products by w^2)
  r_wgcna <- sum(w^2 * (x - mx) * (y - my)) / sqrt(sum(w^2 * (x - mx)^2) * sum(w^2 * (y - my)^2))
  px <- tapply(wn, x, sum); py <- tapply(wn, y, sum)
  b <- bounds_from_props(as.numeric(names(px)), as.numeric(px), as.numeric(names(py)), as.numeric(py))
  list(n = length(x), K1 = length(px), K2 = length(py), r_w = r_w, r_wgcna = r_wgcna,
       r_min = unname(b["r_min"]), r_max = unname(b["r_max"]))
}

vars <- c("ideoid", "partyid", issue_vars)
pairs <- t(combn(vars, 2))
years <- sort(unique(d$year))
cat("Computing", nrow(pairs), "pairs x", length(years), "years...\n")

res <- rbindlist(lapply(years, function(yr) {
  dy <- d[year == yr]
  rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    s <- pair_stats(dy[[pairs[i, 1]]], dy[[pairs[i, 2]]], dy$faceWT)
    if (is.null(s)) return(NULL)
    c(list(year = yr, var1 = pairs[i, 1], var2 = pairs[i, 2]), s)
  }))
}))

res[, type := fifelse(var1 %in% c("partyid", "ideoid") & var2 %in% c("partyid", "ideoid"), "party_ideo",
              fifelse(var1 == "partyid" | var2 == "partyid", "issue_party",
              fifelse(var1 == "ideoid" | var2 == "ideoid", "issue_ideo", "issue_issue")))]
res[, `:=`(dom1 = domain_of[var1], dom2 = domain_of[var2])]

viol <- res[r_w < r_min - 1e-9 | r_w > r_max + 1e-9]
cat(sprintf("Pair-years: %d | years: %s\n", nrow(res), paste(range(res$year), collapse = "-")))
cat(sprintf("Bound violations (should be 0): %d\n", nrow(viol)))
cat(sprintf("max |r_w - r_wgcna|: %.4f; mean: %.4f\n",
            max(abs(res$r_w - res$r_wgcna)), mean(abs(res$r_w - res$r_wgcna))))

saveRDS(res, "output/data/anes_pairs_bounds.rds")
cat("Saved to output/data/anes_pairs_bounds.rds\n")

# Item marginals by year (weighted proportions over valid responses), used to
# measure how far response distributions move between waves.
marg <- rbindlist(lapply(years, function(yr) {
  dy <- d[year == yr]
  rbindlist(lapply(vars, function(v) {
    ok <- !is.na(dy[[v]])
    if (sum(ok) < 30) return(NULL)
    p <- tapply(dy$faceWT[ok], dy[[v]][ok], sum)
    data.table(year = yr, item = v, score = as.numeric(names(p)), prop = as.numeric(p) / sum(p),
               n = sum(ok))
  }))
}))
saveRDS(marg, "output/data/anes_item_marginals.rds")
cat("Saved to output/data/anes_item_marginals.rds\n")
