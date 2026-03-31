
# Organic Matter Standing Stock Analysis

Data analysis for organic matter (coarse and fine organic matter + periphyton scrubbed from rocks) ash-free dry mass, chlorophyll-a, natural abundance of carbon and nitrogen stable isotope values, and C:N content ratio.

## 📊 What This Analysis Does

This repository runs linear mixed-effects models with fixed season for the following response variables:

- Ash-free dry mass of CBOM, FBOM, and periphyton (g/m²)
- Chlorophyll-a concentration of periphyton scrubbed from rocks (g/m²)
- Natural abundance of carbon stable isotopes — δ¹³C (‰)
- Natural abundance of nitrogen stable isotopes — δ¹⁵N (‰)
- Carbon to nitrogen content ratio — C:N
- Summary statistics for all above variables plus % carbon and % nitrogen

## 🚀 Quick Start

**Prerequisites**

Install required R packages:

`install.packages(c("lme4", "lmerTest", "dplyr", "ggplot2"))`

**Running the Analysis**

1. Clone or download this repository
2. Open the R project file in RStudio
3. Update file paths to point to your data location
4. Run all chunks or knit the document

## 🔬 Research Background

**Question**

How do organic matter standing stocks, isotopic compositions, and C:N change across a mining-induced salinity gradient?

**Hypotheses**

1. CBOM and FBOM quantity would not differ because all watersheds were predominantly forested with similar stream slopes, but periphyton would increase along the SC gradient from salt subsidies (SO₄²⁻, Ca²⁺, HCO₃⁻).
2. δ¹³C and δ¹⁵N would increase due to more metabolic and salt-tolerant microbial activity and/or different nutrient sources associated with mining activities.
3. C:N would decrease across the SC gradient if auto- and heterotrophic microbial biofilms were subsidized by mining-associated salts.

**Variables**

| Variable | Units | CBOM | FBOM | Periphyton |
|----------|-------|------|------|------------|
| Ash-free dry mass | g/m² | ✓ | ✓ | ✓ |
| Chlorophyll-a | g/m² | | | ✓ |
| δ¹³C | ‰ | ✓ | ✓ | ✓ |
| δ¹⁵N | ‰ | ✓ | ✓ | ✓ |
| C:N ratio | — | ✓ | ✓ | ✓ |
| Specific conductivity | μS/cm | | | |

**Abbreviations**

| Code | Meaning |
|------|---------|
| OM | Organic matter |
| CBOM | Coarse benthic organic matter |
| FBOM | Fine benthic organic matter |
| ALGAE | Periphyton (used in code) |

## 📚 Citation

If you use this code in your research, please cite:

Tabor L.M., Adey A.K., Brown T., Hotchkiss E.R., McLaughlin D., Meehan C., Reid R., Schoelholtz S., Sinning K., Underwood M.E., Zipper C., Entrekin S.A.

## 📧 Contact

For questions about this repository, please open an issue or contact the corresponding authors.

## 📜 License

Please use responsibly and cite appropriately.

---

**Happy analyzing!** 🌊📊🪨
