Revealed Preference under Uniform Choice: Observation Rates and Budget Design
R reproduction package

Author: Hideyasu Sato

Scope
-----
This package reproduces the fixed finite comparisons in the manuscript. It
does not estimate human behavior, simulate choices, optimize arbitrary price
arrays, or introduce new parameter cases.

The input contains 12 cases, three design rules per case, and exhaustive
integer-n comparisons only in C01, C02, C04, and C07. Accounting for duplicate
rules, there are 102 distinct (case, n) probability evaluations. The target is
alpha=4/5; the displayed fixed-repeat probability uses M=10 at each of the two
budgets, for T=20 observations in total.

Requirements
------------
R, gmp, and jsonlite. Tested with R 4.6.1 (Windows UCRT), gmp 0.7.5.1, and
jsonlite 2.0.0. No commercial software and no Rmpfr dependency are required.
The code does not install packages automatically. If needed, install the two
packages in an R session using a CRAN mirror selected by the user:

  install.packages(c("gmp", "jsonlite"))

On systems that cannot install a binary gmp package, its usual build
requirements apply. See the package's official installation documentation.
No claim of testing on Linux, macOS, or older R versions is made.

Execution
---------
Extract or clone this folder, preserving its relative file layout. From this
folder, run either command in a terminal:

  Rscript --vanilla reproduce.R --representative
  Rscript --vanilla reproduce.R --all

No argument defaults to --all. The script can also be invoked by its path
from another working directory; inputs and outputs resolve relative to the
script, not the current working directory.

--representative reproduces C02 and C08 (the six main-text rows), plus C11
(the exact tolerance endpoint). It computes eight distinct probabilities
and nine rule rows, including duplicated n choices at the endpoint.

--all reproduces all 12 fixed cases, all 36 rule rows, and the four specified
integer-n enumerations: 102 distinct evaluations in total. This is the
complete frozen numerical scope, not an expanded comparison exercise.

Outputs
-------
output/representative/ or output/all/ contains:
  comparison.csv     80-significant-digit displays and exact integer M and T
  exact_results.json exact q, exact fixed-M probability, and M certificates
  validation.json    comparisons with fixed exact references and endpoints
  session_info.txt   actual R and package versions used in that run

Outputs are generated, are excluded by .gitignore, and need not be committed.
To read a CSV in R, use colClasses="character" if preserving all displayed
digits matters. Spreadsheet software may silently round long integers. The
JSON encodes exact integers and rational values as strings.

Manuscript mapping
------------------
Main finite-comparison table: C02 and C08, all three rules (six rows).
Supplementary complete comparison table: all 12 cases, three rules each.
Integer-n enumeration table: C01, C02, C04, and C07.
Two-changed-goods and no-neutral-goods formulas: endpoint checks in every
selected case. Sixteen small cases compare the prefix-sum implementation
with a literal two-dimensional positive sum. The strict epsilon=b endpoint
is separately checked to give q=0.

The roles "main", "boundary", and "outside" describe the relation to the
theorem. C12 (R=6/5, epsilon=1/20) is outside the global exponent-matching
range and remains a comparison within the two-budget family only.

Arithmetic and certification
----------------------------
All input fractions are strings converted to gmp bigq. The positive
Gamma-tail sum, q, and [1-(1-q)^10]^2 use exact rational arithmetic. Neutral
mass m=0 is treated separately, without a Gamma or Dirichlet density of
shape zero.

Let C=log(5+2*sqrt(5)) and D=-log(1-q). An integer M reaches alpha=4/5 iff
M*D >= C. Both logarithms are enclosed by exact rational intervals using

  log(x)=2*sum(k>=0, z^(2*k+1)/(2*k+1)), z=(x-1)/(x+1), x>=1,

with the omitted positive tail bounded by its first denominator and a
geometric series. The square root is enclosed using an exact integer square
root. C is range-reduced as 3*log(2)+log((5+2*sqrt(5))/8). Interval endpoints
are rounded outwards to 80 decimal places using integer arithmetic.

The candidate is ceil(C_upper / D_lower), evaluated exactly. It is accepted
only when BOTH strict rational inequalities hold:

  M*D_lower > C_upper
  (M-1)*D_upper < C_lower.

Thus no rounded floating logarithm or rounded probability decides M. If the
fixed precision cannot certify an input, the script stops rather than
silently reporting a minimum. The 80-significant-digit scientific displays
use exact rounding to nearest, ties to even; they are not used in decisions.

Differences from the original implementation
--------------------------------------------
The original finite comparison was written in Python using Fraction and
Decimal. The R version keeps the same positive-sum formula and rational
logarithm enclosure, but replaces Decimal-based candidate generation with
ceil(C_upper/D_lower) in rational arithmetic. Display formatting is also
implemented with exact integer operations. The final two certificate
inequalities are unchanged. All 102 exact q values and all minimum integer
M values are checked against the frozen reference output.

Reference files and verification limits
--------------------------------------
inputs/fixed_spec.json holds only the fixed numerical specifications.
expected/exact_reference.json holds minimal exact reference q and M values,
plus enumeration argmax/argmin sets. The references are checks, never inputs
to the probability formula or M certificate. They were extracted from the
fixed calculations used by the manuscript; internal working notes are not
part of this package.

Matching the references verifies the implementation against those results;
it does not prove the paper's theorems or establish global finite-s design
optimality. The independent endpoint identities provide additional checks
on the probability formula.

License and citation
--------------------
The author-provided code, documentation, fixed input specifications, and
reference values in this repository are distributed under the MIT License.
Copyright (c) 2026 Hideyasu Sato. See LICENSE. Third-party R packages retain
their own licenses and are installed separately.

Related manuscript: Hideyasu Sato, "Revealed Preference under Uniform Choice:
Observation Rates and Budget Design" (2026, manuscript).
No journal publication, DOI, or release identifier is asserted here. When
citing a run, record the repository URL and the exact commit or release used.

Generated output is excluded by .gitignore and may be recreated locally.
The package contains no empirical dataset, original-paper PDF, private review
report, credential, absolute local path, or automatic upload command.
