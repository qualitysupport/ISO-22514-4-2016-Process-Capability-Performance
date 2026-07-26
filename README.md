# ISO 22514-4 Capability & Performance Explorer

An R Shiny tool for calculating and teaching process capability and performance measures following the procedures described in **ISO 22514-4:2016 — Statistical methods in process management — Capability and performance — Part 4: Process capability estimates and performance measures**.

Built as course material for teaching capability/performance concepts alongside the AIAG-VDA / ISO 22514 family of standards. Load your own data, or use the built-in example dataset to explore how capability (Cp/Cpk) and performance (Pp/Ppk) diverge as a process departs from statistical control.

## What it does

- **Concepts explorer (Clause 3):** interactive visualizations of dispersion, reference limits (the 0.135th/99.865th percentiles), the reference interval, and — critically — why short-term (inherent) and total (long-term) dispersion are kept as separate, clearly labeled quantities rather than one number.
- **Capability indices (Clause 4):** Cp, CPU, CPL, and Cpk, calculated from a short-term standard deviation estimated via the subgroup-range method or the moving-range method (Annex A), plus proportion out-of-specification (4.8).
- **Performance indices (Clause 5):** Pp, PPU, PPL, and Ppk, calculated from the total (long-term) standard deviation of all the data — no requirement that the process be shown to be in statistical control.
- **Other target-based indices (Clause 4.7):** Process Capability Fraction (PCF), Mean Square Error (MSE), the Qk index, and Cpm.
- **Non-normal data (Annex C):** skewness/kurtosis diagnostics, a Shapiro-Wilk normality test, the probability-paper percentile method, and maximum-likelihood distribution fitting (log-normal and Weibull) as alternatives to the normal-distribution formulas.
- **Confidence intervals (Annex D):** normal-approximation (formula-method) intervals for Cp, Cpk, Pp, and Ppk, at a confidence level you choose.
- **Downloadable report (Clause 6):** a summary report and the raw data, ready to export as CSV.

Data can be loaded from a built-in example, an uploaded CSV/Excel file, or pasted in manually.

## Running it

Requires R with the following packages:

```r
install.packages(c("shiny", "bslib", "ggplot2", "DT", "readxl", "gridExtra"))
```

`bslib` must be version 0.5.0 or later (needed for the persistent sidebar layout) — the app checks this on startup and will stop with a clear message if your installed version is too old. `MASS` is also used but ships with base R, so it doesn't need a separate install.

Then:

```r
shiny::runApp("app.R")
```

## A note on capability vs. performance

Cp/Cpk and Pp/Ppk are often reported as if they're interchangeable — they aren't. Capability indices assume the process has been demonstrated to be in statistical control and use only the inherent, short-term variation; performance indices make no such assumption and use the total variation actually observed, including any drift or instability. Performance will typically be lower than capability for the same data, and a big gap between the two is itself diagnostic — it usually means the process isn't as stable as a control-chart-only view would suggest. The app calculates both side by side for exactly this reason.

## A note on non-normal data

The normal-distribution formulas in Clause 4/5 assume your data follows a bell curve. If it doesn't, ISO 22514-4 (Clause 4.5/5.3, Annex C) provides alternatives: reading percentiles off a probability plot, or fitting a specific distribution family and using its theoretical percentiles instead. This app implements both — the probability-paper method (via empirical percentiles) and maximum-likelihood fits for the log-normal and Weibull families, the two most common in capability work. Percentile estimates from either method get unstable on small samples (well under ~100 points), so treat them as indicative rather than precise until your sample size grows.

## Disclaimer

This tool implements calculation methods described in ISO 22514-4:2016 but is not affiliated with, endorsed by, or reviewed by ISO. It does not reproduce the standard itself. For authoritative guidance, consult the official ISO 22514-4:2016 document, available for purchase from [ISO](https://www.iso.org/) or your national standards body.

Provided as-is, for educational and reference use. The "Capable / Marginal / Not capable" labels shown next to Cp/Cpk/Pp/Ppk (using the common 1.33 / 1.00 thresholds) are industry convention, not values mandated by the standard itself — ISO 22514-4 does not itself set minimum acceptable index values; those are a matter of agreement between the parties involved.

## Author

Dan Lay Jr.
Metrologist | ASQ Certified Calibration Technician | Calibration Support LLC
[www.calibrationsupport.com](https://www.calibrationsupport.com) · [LinkedIn](https://linkedin.com/in/dlayjr)

## License

MIT — see [LICENSE](https://github.com/qualitysupport/iso-22514-4-capability-performance/blob/main/LICENSE).
