# Revealed Preference under Uniform Choice

R code for the fixed numerical comparisons in **“Revealed Preference under Uniform Choice: Observation Rates and Budget Design”**, by **Hideyasu Sato**.

The code evaluates a two-budget design under independent uniform expenditure shares. It computes exact threshold-crossing probabilities and certifies the minimum integer number of repetitions needed to reach the specified rejection probability. It uses deterministic arithmetic, not simulated choices.

## Files

| File | Purpose |
| --- | --- |
| [`reproduce.R`](reproduce.R) | Calculation, certification, and validation script |
| [`inputs/fixed_spec.json`](inputs/fixed_spec.json) | The 12 fixed parameter cases and design rules |
| [`expected/exact_reference.json`](expected/exact_reference.json) | Exact reference results for validation |
| [`DEPENDENCIES.txt`](DEPENDENCIES.txt) | Tested software and package versions |
| [`README.txt`](README.txt) | Detailed instructions and numerical methods |
| [`.gitignore`](.gitignore) | Excludes generated output and local R files |
| [`LICENSE`](LICENSE) | MIT license for the repository materials |

## Requirements

Tested with **R 4.6.1**, **gmp 0.7.5.1**, and **jsonlite 2.0.0** on Windows. No empirical data or commercial software are required. Other operating systems and R versions have not been tested.

Install missing packages in R:

```r
install.packages(c("gmp", "jsonlite"))
```

The reproduction script does not install packages automatically.

## Run

Download or clone the repository, preserving its folder structure. From its root, run in a terminal:

```sh
Rscript --vanilla reproduce.R --representative
```

This reproduces the six main-text rows (C02 and C08) and the tolerance-boundary case C11: nine design-rule rows and eight distinct probability evaluations.

For the complete fixed comparison:

```sh
Rscript --vanilla reproduce.R --all
```

The complete run covers 12 cases, 36 design-rule rows, and the four specified integer-design enumerations: **102 distinct evaluations**. Running without an argument also selects `--all`. Paths resolve relative to the script.

## Outputs and manuscript mapping

Results are written to `output/representative/` or `output/all/`:

| Output | Contents |
| --- | --- |
| `comparison.csv` | Probabilities displayed to 80 significant digits and exact integer repetition counts |
| `exact_results.json` | Exact rational probabilities and certificates for the minimum repetition counts |
| `validation.json` | Reference comparisons and boundary checks |
| `session_info.txt` | R and package versions used for the run |

The main-text comparison uses **C02 and C08**. The complete supplementary comparison uses all 12 cases. The integer-design enumeration uses **C01, C02, C04, and C07**.

The target rejection probability is **α = 4/5**. The fixed-repeat comparison uses **M = 10 observations at each budget**, giving **T = 20** in total. C12 lies outside the theorem's global exponent-matching range and is a comparison within the two-budget family only.

Generated output is excluded by `.gitignore`; it can be recreated with the commands above. To preserve long integers and displayed digits when reading the CSV in R, use `colClasses = "character"`. Exact values are stored as strings in JSON.

## Verification and arithmetic

All 102 exact probabilities, fixed-repeat rejection probabilities, and minimum integer repetition counts were checked against the original fixed calculations and matched. The script also checks analytical endpoints and 16 small cases using an alternative summation implementation.

Probabilities use `gmp` exact rational arithmetic. Minimum repetition counts are certified using outward rational bounds for logarithms; rounded displays do not determine the integer result. The script stops if its precision is insufficient for certification. Reference results are used only for validation, not as inputs to the probability calculations.

See [`README.txt`](README.txt) for the formulas, precision rules, and verification limits.

## Interpretation

The probabilities describe the specified independent uniform-share benchmark. They are not estimates of human choice behavior or guarantees of power under an unknown behavioral distribution. Reproducing the numerical comparisons does not prove the paper's theorems or establish globally optimal designs at a finite number of goods.

## License

The author-provided code, documentation, fixed specifications, and reference values are released under the [MIT License](LICENSE). Third-party R packages retain their respective licenses.
