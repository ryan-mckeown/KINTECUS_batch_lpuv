# =====================================================================
# photolysis.R -- chemistry helpers for the UV/H2O2/chloramine model
# =====================================================================
# Everything tunable lives in CONST / PHOTO below. Nothing else in the
# pipeline hard-codes a constant.

CONST <- list(
  pKw       = 14,      # textbook 25C = 14.00; target model uses 13.967 (Kw=1.079e-14)
  pK1       =  6.32,   # H2CO3* <-> HCO3- + H+   (target 6.32; textbook 6.35)
  pK2       = 10.20,   # HCO3-  <-> CO3-- + H+   (target 10.20; textbook 10.33)
  pKa_HOCl  =  7.53,   # HOCl   <-> OCl- + H+
  pKa_NH4   =  9.30,   # NH4+   <-> NH3  + H+
  MW_Cl2    = 70.90,   # free chlorine & NH2Cl, mg/L as Cl2 -> M
  MW_NHCl2  = 141.80,  # NHCl2 carries 2 Cl → 2 × 70.90 per mole as Cl2
  MW_CaCO3  = 100.09,
  MW_N      = 14.01,
  MW_Cl     = 35.453,
  MW_C      = 12.011,
  MW_O2     = 31.998,
  MW_H2O2   = 34.0147,  # mg/L as H2O2 -> M
  MW_DIOX   = 88.106    # 1,4-dioxane C4H8O2, ug/L -> M
)

# 254 nm molar absorptivity (L/mol/cm) and quantum yield (mol/Einstein)
PHOTO <- data.frame(
  species = c("H2O2", "HOCl", "OCl-",  "NH2Cl", "NHCl2"),
  eps     = c(  18.6,   58.0,   60.0,    371.0,   126.0),
  qy      = c(  0.50,   0.55,   0.55,    0.294,    0.82),
  stringsAsFactors = FALSE
)

# --- speciation -------------------------------------------------------
# Each returns molar concentrations for the species it owns.

free_chlorine <- function(pH, FAC_mgL) {
  if (!isTRUE(FAC_mgL > 0)) return(c(HOCl = 0, `OCl-` = 0))
  f  <- 10^-pH / (10^-pH + 10^-CONST$pKa_HOCl)   # HOCl fraction
  Ct <- FAC_mgL / 1000 / CONST$MW_Cl2
  c(HOCl = Ct * f, `OCl-` = Ct * (1 - f))
}

# Alkalinity (mg/L as CaCO3) -> full carbonate speciation at this pH.
#   Alk_eq = Ct*(a1 + 2*a2) + [OH-] - [H+]
# so Ct = (Alk_eq + [H+] - [OH-]) / (a1 + 2*a2).
# The +[H+] term is NOT optional here: at pH 5.6 with low-alkalinity RO
# permeate it is ~4% of the alkalinity, and dropping it put total carbonate
# 4.3% below the target model's. With it, t0 carbonate matches to <0.1%.
carbonate <- function(pH, Alk_mgL) {
  if (!isTRUE(Alk_mgL > 0)) return(c(H2CO3 = 0, `HCO3-` = 0, `CO3--` = 0))
  H  <- 10^-pH; OH <- 10^(pH - CONST$pKw)
  K1 <- 10^-CONST$pK1; K2 <- 10^-CONST$pK2
  den <- H^2 + H*K1 + K1*K2
  a0 <- H^2/den; a1 <- H*K1/den; a2 <- K1*K2/den
  Alk_eq <- Alk_mgL / 1000 / CONST$MW_CaCO3 * 2
  Ct <- max(Alk_eq + H - OH, 0) / (a1 + 2*a2)
  c(H2CO3 = Ct*a0, `HCO3-` = Ct*a1, `CO3--` = Ct*a2)
}

ammonia <- function(pH, NH3_mgL_N) {
  if (!isTRUE(NH3_mgL_N > 0)) return(c(NH3 = 0, `NH4+` = 0))
  f  <- 10^-pH / (10^-pH + 10^-CONST$pKa_NH4)    # NH4+ fraction
  Ct <- NH3_mgL_N / 1000 / CONST$MW_N
  c(NH3 = Ct * (1 - f), `NH4+` = Ct * f)
}

chloramines <- function(NH2Cl_mgL, NHCl2_mgL) {
  c(NH2Cl = if (isTRUE(NH2Cl_mgL > 0)) NH2Cl_mgL/1000/CONST$MW_Cl2   else 0,
    NHCl2 = if (isTRUE(NHCl2_mgL > 0)) NHCl2_mgL/1000/CONST$MW_NHCl2 else 0)
}

# --- photolysis -------------------------------------------------------
# Optically-dilute (Taylor-expanded) Beer-Lambert, matching the target
# model's 2.303*QY*eps*E*1000 form:  k = ln(10)*QY*I0*eps*L / K
#
# Irradiance ratio.  As absorbers are consumed the cell transmits more, so
# the rate constants rise.  The ratio is built from the PATH-AVERAGED
# fluence rate -- what a well-mixed cell actually sees -- not from the
# intensity at the exit plane:
#
#     Ebar(A, L) = (1 - 10^-(A*L)) / (A*L)
#
# For small A*L this goes as 1 + dA*L*ln10/2, whereas the exit-plane form
# 10^(dA*L) goes as 1 + dA*L*ln10 -- exactly twice the exponent.  Using the
# exit-plane form put I_ratio at 1.0793 where the target has 1.0385.
#
# Only the photolysing species enter this ratio; background absorbance is
# deliberately excluded (the target does the same).  Under the old
# exit-plane form A_bg cancelled out; here it would NOT, so it is kept out
# on purpose and used for reporting UVT only.
Ebar <- function(A, L) if (A * L < 1e-12) log(10) else (1 - 10^(-A * L)) / (A * L)

k_photo <- function(conc, E_avg, path_cm, A_bg, A_sp0 = NA) {
  A_sp <- sum(PHOTO$eps * vapply(PHOTO$species, function(s) conc[[s]], 0))
  if (is.na(A_sp0)) A_sp0 <- A_sp
  I_ratio <- Ebar(A_sp, path_cm) / Ebar(A_sp0, path_cm)
  I0      <- E_avg * I_ratio / (470697.6 * path_cm)   # Einstein/L/s
  k <- setNames(log(10) * PHOTO$qy * I0 * PHOTO$eps * path_cm, PHOTO$species)
  A_path <- (A_bg + A_sp) * path_cm                   # reporting only
  list(k = k, A_path = A_path, I_ratio = I_ratio, A_sp = A_sp,
       UVT_pct = 100 * 10^-A_path, E_now = E_avg * I_ratio)
}

# --- breakpoint / chloramine rate constants ---------------------------
# Formula numbers follow the source spreadsheet's "FORMULA" placeholders.
# 214 excludes [NHCl2] (Kintecus supplies it as the reaction's reactant).
# 207: k207 = 167*[OH-], now resolved -- OH- is used directly from conc.
k_breakpoint <- function(T_K, conc) {
  g  <- function(s) { v <- conc[[s]]; if (is.null(v) || is.na(v)) 0 else v }
  H  <- max(g("H+"), 1e-14); OH <- max(g("OH-"), 1e-14)
  HCO3 <- g("HCO3-"); H2CO3 <- g("H2CO3"); CO3 <- g("CO3--")
  HOCl <- g("HOCl");  OCl   <- g("OCl-")

  list(
    k201 = 2.04e9 * exp(-1887 / T_K),                      # HOCl + NH3 -> NH2Cl
    k202 = 1.38e8 * exp(-8800 / T_K),                      # NH2Cl -> HOCl + NH3
    k203 = 3.0e5  * exp(-2010 / T_K),                      # HOCl + NH2Cl -> NHCl2
    k205 = (3.78e10 * exp(-2169  / T_K) / 3600) * H +      # 2NH2Cl -> NHCl2 + NH3
           (0.87    * exp(-503   / T_K) / 3600) * HCO3 +
           (2.52e25 * exp(-16860 / T_K) / 3600) * H2CO3,
    k207 = 167 * OH,                                        # NHCl2 -> NOH + 2H+ + 2Cl-
    k211 = 3.28e9 * OH + 9.00e4 * (10^-7.5 * HOCl / H) + 6.00e6 * CO3,
    k212 = 5.56e10 * OH,                                   # NHCl2 + NCl3
    k213 = 1.39e10 * OH,                                   # NH2Cl + NCl3
    k214 = 66.0 * OCl,                                     # NHCl2 -> NO3-
    k215 = 1.60e-6 + 8.0*OH + 890.0*OH^2 + 65.0*HCO3*OH    # NCl3 -> NHCl2 + HOCl
  )
}
