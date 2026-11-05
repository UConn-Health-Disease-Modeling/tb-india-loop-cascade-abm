"""Calibration targets.

Change from the earlier single-year calibration: the three cross-sectional 2020 burden targets are
replaced by a four-point WHO incidence TRAJECTORY plus the 2020 prevalence and
mortality anchors, so the calibration constrains the direction of the epidemic,
not only its level.

Incidence series: WHO Global TB Programme estimates for India, retrieved from
the World Bank mirror of the WHO series (indicator SH.TBS.INCD), per 100,000/yr:

    2000 322 | 2005 311 | 2010 276 | 2015 237 | 2020 197 | 2023 192

Four points (2000, 2010, 2015, 2020) are used as targets; 2005 and 2023 are held
out so the fit can be checked against unused years.  sigma is set at 17% of the
point estimate, matching the ratio used for the earlier 2020 target
(32 / 188) and consistent with the width of WHO's own uncertainty intervals for
India.

NOTE ON VINTAGE: the earlier calibration used incidence_2020 = 188 (WHO Global
TB Report 2021 vintage); the current WHO series puts 2020 at 197.  The whole
trajectory here comes from one vintage so the trend is internally consistent;
the difference is well inside one sigma.

Split:
  CALIBRATION (emulator-gated) = the 6 burden targets below
  VALIDATION  (best-fit 3-sigma check only) = the 7 cascade targets
"""
from __future__ import annotations

# WHO / World Bank SH.TBS.INCD, India, per 100,000/yr
WHO_INCIDENCE = {2000: 322, 2005: 311, 2010: 276, 2015: 237, 2020: 197, 2023: 192}
SIGMA_FRACTION = 0.17

TREND_YEARS = [2000, 2010, 2015, 2020]
HOLDOUT_YEARS = [2005, 2023]

BURDEN_TARGETS = {
    f"incidence_{y}": (WHO_INCIDENCE[y], round(WHO_INCIDENCE[y] * SIGMA_FRACTION, 1))
    for y in TREND_YEARS
}
BURDEN_TARGETS.update({
    "prevalence_2020": (304, 36.1),     # National TB Prevalence Survey 2019-21
    "mortality_2020":  (33,  6.0),      # WHO 2020
})

CASCADE_TARGETS = {
    "noloop_frac":               (0.20, 0.06),
    "testloop_frac":             (0.50, 0.10),
    "testtxfail_frac":           (0.30, 0.07),
    "tests_per_case_overall":    (3.5,  0.78),
    "tests_per_case_noloop":     (1.0,  0.42),
    "tests_per_case_testloop":   (3.0,  0.92),
    "tests_per_case_testtxfail": (4.5,  1.28),
}

TARGETS = {**BURDEN_TARGETS, **CASCADE_TARGETS}
BURDEN = list(BURDEN_TARGETS)
CASCADE = list(CASCADE_TARGETS)

# Derived diagnostics reported alongside the fit, not gated on
DERIVED = {
    # annualised % change in incidence implied by the targets
    "incidence_annual_change_2000_2020":
        ((WHO_INCIDENCE[2020] / WHO_INCIDENCE[2000]) ** (1 / 20) - 1) * 100,
}


if __name__ == "__main__":
    print("Calibration (emulator-gated) targets:")
    for k, (z, s) in BURDEN_TARGETS.items():
        print(f"  {k:22s} {z:>7} +/- {s}")
    print("\nValidation (3-sigma check) targets:")
    for k, (z, s) in CASCADE_TARGETS.items():
        print(f"  {k:26s} {z:>6} +/- {s}")
    print("\nWHO implied trend 2000-2020: %+.2f%%/yr"
          % DERIVED["incidence_annual_change_2000_2020"])
    print("earlier model without secular improvement, 2000-2020: +0.05%/yr")
