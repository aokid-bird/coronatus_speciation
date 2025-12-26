# =================#
#### libraries ####
# =================#
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(readr)
library(vcfR)
library(data.table)
library(future.apply)

# ==================#
#### parameters ####
# ==================#
bamlist.path <- snakemake@input[["bamlist"]]
vcf.path <- snakemake@input[["vcf"]]
out <- snakemake@output[["raxmlcatg"]]
.na_val <- snakemake@config[["raxml"]][["na_fill"]] %>% as.numeric()
threads <- snakemake@threads
var_option <- snakemake@config[["raxml"]][["var_option"]]
maxSize <- snakemake@params[["maxSize"]] %>% as.numeric()
out.path <- dirname(out)

options(future.globals.maxSize = maxSize)

# ============#
#### data ####
# ============#
pos.lib <- tibble(
  alleles = c(
    "A/C", "A/G", "A/T", "C/G", "C/T", "G/T",
    "C/A", "G/A", "T/A", "G/C", "T/C", "T/G"
  ),
  position_ref = c(
    "1,5,2", "1,6,3", "1,7,4", "2,8,3", "2,9,4", "3,10,4",
    "2,5,1", "3,6,1", "4,7,1", "3,8,2", "4,9,2", "4,10,3"
  )
)

ambig.lib <- tibble(
  alleles = c(
    "A/A", "C/C", "G/G", "T/T",
    "A/C", "A/G", "A/T", "C/G", "C/T", "G/T",
    "C/A", "G/A", "T/A", "G/C", "T/C", "T/G",
    "."
  ),
  genotype = c(
    "A", "C", "G", "T",
    "M", "R", "W", "S", "Y", "K",
    "M", "R", "W", "S", "Y", "K",
    "N"
  )
)

# =================#
#### functions ####
# =================#
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall = k))

gp2ambig.dt <- function(gt_GP, gt_GT, position_ref, .decimal = 15, na_val = 0.00, ...) {
  if (is.na(gt_GT)) {
    p <- str_c(rep(as.character(na_val), 10), collapse = ",")
  } else {
    v <- data.table(pos = 1:10, val = "0.0")
    val <-
      c(position_ref, gt_GP) %>%
      str_split(",", simplify = TRUE) %>%
      t() %>%
      as.data.table() %>%
      .[, .(pos = as.integer(V1), GP = as.numeric(V2))] %>%
      .[, GP := as.numeric(specify_decimal(GP, .decimal))]

    val[GP == max(GP), GP := 1 - sum(GP[GP != max(GP)])]

    val[, GP := format(GP, scientific = FALSE)]

    v <- merge(v, val, by = "pos", all.x = TRUE)
    v[, val := fifelse(is.na(GP), val, GP)] # valの置換
    v[, val.num := as.numeric(val)]
    v[, val := fifelse(val.num == 0, "0.0", val)]

    p <- paste(v$val, collapse = ",")
  }

  p
}

# ===========#
#### Run ####
# ===========#
# read bamlist
tmp.bamlist <- read_delim(str_interp("${bamlist.path}"), delim = " ", col_names = FALSE)
tmp.bamlist.line <- tmp.bamlist %>% pull(X1)

cat("bamlist was successfully read\n\n")

ID <-
  tmp.bamlist.line %>%
  str_remove_all(., ".*/|.bam")

# read vcf.gz by vcfr
vcf <- read.vcfR(str_interp("${vcf.path}"), verbose = FALSE)
# convert to tidy format
Z <- vcfR2tidy(vcf, format_fields = c("GT", "GP"))

fix.lib <-
  Z$fix %>%
  select(ChromKey, POS, REF, ALT) %>%
  mutate(
    REF_REF_state = str_c(REF, "/", REF),
    REF_ALT_state = str_c(REF, "/", ALT),
    ALT_ALT_state = str_c(ALT, "/", ALT)
  ) %>%
  select(ChromKey, POS, REF_ALT_state) %>%
  left_join(., pos.lib, by = c("REF_ALT_state" = "alleles"))

comb.dt <- comb.df <-
  Z$gt %>%
  left_join(., fix.lib, by = c("ChromKey", "POS")) %>%
  left_join(., ambig.lib, by = c("gt_GT_alleles" = "alleles"))

# rm(list = c("Z", "vcf"));gc(reset=TRUE)

setDT(comb.dt)
# a <- Sys.time()

plan(multisession, workers = threads)
comb.dt[, gtambig := future_mapply(gp2ambig.dt, gt_GP, gt_GT, position_ref, na_val = .na_val, MoreArgs = list(.decimal = 15))]

cat("gp2ambig.dt successfully worked\n\n")

setattr(comb.dt, "class", c("tbl_df", "tbl", "data.frame"))

concat_GP.df <-
  comb.dt %>%
  select(ChromKey, POS, Indiv, gtambig) %>%
  pivot_wider(names_from = Indiv, values_from = gtambig) %>%
  select(ChromKey, POS, tmp.bamlist.line) %>%
  unite("concat_GP", tmp.bamlist.line, sep = "\t", remove = TRUE)
# b <- Sys.time()
# b-a
write_csv(concat_GP.df, str_interp("${out.path}/concat_GP.csv"))

concat_geno.df <-
  comb.df %>%
  select(ChromKey, POS, Indiv, genotype) %>%
  pivot_wider(names_from = Indiv, values_from = genotype) %>%
  select(ChromKey, POS, tmp.bamlist.line) %>%
  unite("concat_GP", tmp.bamlist.line, sep = "", remove = TRUE)

# check if the site is variable at the genotype inference
check.var <-
  concat_geno.df$concat_GP %>%
  lapply(str_split, pattern = "") %>%
  lapply(unlist) %>%
  lapply(unique) %>%
  lapply(length) %>%
  unlist()

to_be_removed <- which(check.var == 1)
# remove those estimated as monomorphic at the genotype inference if the option is "onlyvar"
if (length(to_be_removed) > 0 && var_option == "onlyvar") {
  concat_geno.df <- concat_geno.df[-c(to_be_removed), ]
}

# combine
catg <-
  left_join(concat_geno.df, concat_GP.df, by = c("ChromKey", "POS")) %>%
  select(!c(ChromKey, POS)) %>%
  unite("all", everything(), sep = "\t", remove = TRUE)

catg.full <-
  bind_rows(
    tibble(all = str_c(length(tmp.bamlist.line), nrow(catg), sep = "\t")),
    tibble(all = str_c(ID, collapse = "\t")),
    catg
  )

write_delim(catg.full, str_interp("${out}"), col_names = FALSE, quote = "none")
cat("raxml CATG input was successfully generated!\n\n")
