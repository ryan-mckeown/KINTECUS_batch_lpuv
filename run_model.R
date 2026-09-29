# =====================================================================
# run_model.R -- batch driver: substep loop around Kintecus
# =====================================================================
# For each experiment row: seed concentrations, then repeatedly
#   (1) recompute rate constants from the CURRENT concentrations
#   (2) run Kintecus for dt seconds
#   (3) carry EVERY species forward into the next substep
# Concentrations are written straight into SPECIES.DAT -- there are no
# per-species text placeholders, so a species can never be silently
# dropped or reset to zero.

library(readxl)
source("photolysis.R")

KINTECUS <- "C:/KINTECUS/kintecus.exe"
DT       <- 1.0          # substep length (s)
PINNED   <- c("H+", "OH-")   # held constant by Kintecus for the whole run

if (.Platform$OS.type == "windows")
  try(system("taskkill /F /IM kintecus.exe", show.output.on.console = FALSE), silent = TRUE)
if (file.exists("C:/KINTECUS/KINKEY.TXT") && !file.exists("KINKEY.TXT"))
  file.copy("C:/KINTECUS/KINKEY.TXT", "KINKEY.TXT")

num <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || is.na(x)) default else x
}

# --- SPECIES.DAT: the template is just the species list -----------------
sp_lines <- readLines("SPECIES_template.DAT", warn = FALSE)
sp_lines <- trimws(sp_lines[!grepl("^\\s*#", sp_lines) & trimws(sp_lines) != ""])
SPECIES  <- trimws(sub("\t.*$", "", sp_lines))
SPECIES  <- SPECIES[SPECIES != "END"]
stopifnot(!any(duplicated(SPECIES)))

write_species <- function(conc) {
  rows <- sprintf("%s\t0\t%.6E\tY\t0\tNo\t%s",
                  SPECIES,
                  conc[SPECIES],
                  ifelse(SPECIES %in% PINNED, sprintf("%.6E", conc[SPECIES]), "No"))
  writeLines(c("# Species\tResidence\tInitial\tDisplay\tExternal\tSpecial\tConstant File?",
               rows, "END"), "SPECIES.DAT")
}

fill_template <- function(file_in, file_out, subs) {
  txt <- readLines(file_in, warn = FALSE)
  for (k in names(subs)) txt <- gsub(k, subs[[k]], txt, fixed = TRUE)
  writeLines(txt, file_out)
}

# --- main loop --------------------------------------------------------
experiments <- read_excel("kintecus_experimental_batch.xlsx")
col <- function(nm, i, default = 0)
  if (nm %in% names(experiments)) num(experiments[[nm]][i], default) else default

for (i in seq_len(nrow(experiments))) {
  id   <- experiments$Exp_ID[i]
  hrt  <- col("HRT", i, 1); E_avg <- col("E_avg", i)
  pH   <- col("pH", i);     path  <- col("path_length", i)
  T_K  <- if ("Temp_K" %in% names(experiments)) {
    col("Temp_K", i, 298.15)
  } else if ("Temp_C" %in% names(experiments)) {
    col("Temp_C", i, 25) + 273.15
  } else 298.15
  if (path <= 0) stop(sprintf("Exp %s: path_length must be > 0", id))

  # every species starts at zero, then the dosed ones are filled in
  conc <- setNames(numeric(length(SPECIES)), SPECIES)
  seed <- c(free_chlorine(pH, col("FAC_mgL", i)),
            carbonate(pH, col("Alk_mgL", i)),
            ammonia(pH, col("NH3_mgL_N", i)),
            chloramines(col("NH2Cl_mgL", i), col("NHCl2_mgL", i)),
            c(H2O2    = col("H2O2_0", i),
              DIOXANE = col("DIOXANE_0", i),
              H2O     = 55.56,
              `H+`    = 10^-pH,
              `OH-`   = 10^(pH - CONST$pKw),
              `Cl-`   = col("Cl_mgL",  i) / 1000 / CONST$MW_Cl,
              TOC     = col("TOC_mgL", i) / 1000 / CONST$MW_C,
              O2      = col("O2_mgL",  i) / 1000 / CONST$MW_O2,
              `NO3-`  = col("NO3_mgL", i) / 1000 / CONST$MW_N,
              `NO2-`  = col("NO2_mgL", i) / 1000 / CONST$MW_N))
  missing <- setdiff(names(seed), SPECIES)
  if (length(missing))
    stop("Seeded species absent from SPECIES_template.DAT: ", paste(missing, collapse = ", "))
  conc[names(seed)] <- seed

  # Fixed background absorbance = measured UVT minus what the model absorbs at t0.
  uvt   <- if ("UVT" %in% names(experiments)) col("UVT", i, 100) else col("UVT_pct", i, 100)
  A_sp0 <- sum(PHOTO$eps * conc[PHOTO$species])
  A_bg  <- max(-log10(max(uvt, 1e-4) / 100) - A_sp0, 0)

  cat(sprintf("\n=== %s | %.1f s | pH %.2f | %.2f cm | %.1f K | A_bg %.5f ===\n",
              id, hrt, pH, path, T_K, A_bg))

  history <- list(); t_now <- 0
  while (t_now < hrt - 1e-9) {
    dt <- min(DT, hrt - t_now)

    ph <- k_photo(conc, E_avg, path, A_bg, A_sp0)
    bp <- k_breakpoint(T_K, conc)

    fill_template("PARM_template.DAT", "PARM.DAT", list(
      "_HRT_"      = sprintf("%.6E", dt),
      "_SAMP_INT_" = sprintf("%.6E", dt / 10),
      "_TEMP_K_"   = sprintf("%d", round(T_K))))

    fill_template("MODEL_template.DAT", "MODEL.DAT", c(
      setNames(as.list(sprintf("%.4E", ph$k)),
               c("_K_H2O2_", "_K_HOCL_", "_K_OCL_", "_K_NH2CL_", "_K_NHCL2_")),
      setNames(as.list(sprintf("%.4E", unlist(bp))),
               paste0("_K_", sub("^k", "", names(bp)), "_"))))

    write_species(conc)

    if (file.exists("CONC.TXT")) file.remove("CONC.TXT")
    system2(KINTECUS, "-ig:mass", stdout = "", stderr = "", wait = TRUE)
    if (!file.exists("CONC.TXT"))
      stop(sprintf("Kintecus produced no output at t=%.1f s (Exp %s)", t_now, id))

    out <- read.table("CONC.TXT", header = TRUE, check.names = FALSE)
    names(out) <- trimws(names(out))
    tcol <- grep("^time", names(out), ignore.case = TRUE)[1]
    out[[tcol]] <- out[[tcol]] + t_now

    # carry EVERY species forward -- nothing is dropped or zeroed
    got <- intersect(SPECIES, names(out))
    conc[got] <- as.numeric(unlist(out[nrow(out), got]))

    out$I_ratio    <- ph$I_ratio
    out$E_now      <- ph$E_now
    out$UVT_pct    <- ph$UVT_pct
    out$A_path     <- ph$A_path
    for (s in PHOTO$species) out[[paste0("k_photo_", s)]] <- ph$k[[s]]

    history[[length(history) + 1]] <- if (length(history)) out[-1, ] else out
    t_now <- t_now + dt
  }

  run <- do.call(rbind, history)
  tcol <- grep("^time", names(run), ignore.case = TRUE)[1]
  # Fluence, as the target defines it: the CURRENT irradiance times elapsed
  # time, not a cumulative integral of a rising irradiance.  Verified against
  # the target workbook -- its dose column divided by (E_avg * t) equals its
  # I ratio on every row.  Because I_ratio rises, this product sits ~1% above
  # the integral of the same trajectory.
  run$Fluence_mJ_cm2 <- E_avg * run$I_ratio * run[[tcol]]
  c0 <- num(experiments[["DIOXANE_0"]][i])
  run$LRV <- ifelse(run$DIOXANE > 0, log10(c0 / run$DIOXANE), NA)

  f <- sprintf("CONC_Exp_%s.csv", id)
  write.csv(run, f, row.names = FALSE)
  cat(sprintf("   %s | LRV %.3f | fluence %.0f mJ/cm2 | I ratio %.4f\n",
              f, tail(run$LRV, 1), tail(run$Fluence_mJ_cm2, 1), tail(run$I_ratio, 1)))
}

cat("\nDone.\n")
