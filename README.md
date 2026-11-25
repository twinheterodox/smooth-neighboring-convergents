# Smooth convergents and the size of the next denominator

twinheterodox<br>
Prepared with Codex

Suppose an irrational has infinitely many convergents $P/q$ for which the prime factors of $|P|$, $q$, and the next denominator $Q$ belong to a fixed set $S$. What is the largest exponent in Marques’s lower bound for $\log Q$? For two primes the sharp exponent is 1. For three primes the supremum is at most 1/2. The square-root scale occurs infinitely often, and smaller minima at other powers could lower the exponent further. [Read the paper (PDF)](output/pdf/smooth-neighboring-convergents.pdf)

At the fixed convergent $65{,}536/59{,}049$, the least possible $Q$ has 5,921 digits if only 2 is allowed. Allowing 5 as well reduces the minimum to 19 digits. Figure 1 shows the exponent pairs for this example. [Figure source](paper/figures/signed-product-exponents.tex) · [Supplementary continued-fraction diagram (PDF)](output/pdf/continued-fraction-branches.pdf)

## Reproduce the finite checks

Run the following commands from the repository root with Python 3. They use only the standard library.

```sh
python3 checks/check_examples.py
python3 checks/check_parents.py
python3 checks/independent-checks.py
```

The programs enumerate smooth products, reconstruct consecutive convergents, and check the arithmetic in the examples. [Verification and review scope](review/verification.md)

## Build the paper

Run pdfLaTeX twice from the repository root with a TeX Live installation that includes TikZ. The second pass resolves the cross-references. [Manuscript source](paper/smooth-neighboring-convergents.tex)

```sh
mkdir -p build
pdflatex -no-shell-escape -halt-on-error -output-directory=build paper/smooth-neighboring-convergents.tex
pdflatex -no-shell-escape -halt-on-error -output-directory=build paper/smooth-neighboring-convergents.tex
```
