# Congruences modulo powers of 5 for broken 12-diamond partitions

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24521408.svg)](https://doi.org/10.5281/zenodo.24521408)

This repository holds the paper *Congruences modulo powers of 5 for broken 12-diamond partitions* by Ramanuja, together with code that recomputes its tables.

Let Δ₁₂(n) be the number of broken 12-diamond partitions of n, and let r_α = (11·25^α + 13)/12. The paper proves the following, for all n ≥ 0:

- **Theorem 1.1.** For α ≥ 1:
  - Δ₁₂(25^α n + r_α) ≡ 0 (mod 5^α);
  - Δ₁₂(25^(α+1) n + r_(α+1)) + 5Δ₁₂(25^α n + r_α) ≡ 0 (mod 5^(α+2)).
- **Theorem 1.2.** For α ≥ 2, Δ₁₂(25^(α+1) n + r_(α+1)) + 5Δ₁₂(25^α n + r_α) − 125Δ₁₂(25^(α−1) n + r_(α−1)) ≡ 0 (mod 5^(α+3)).
- **Theorem 1.3.** Modulo the next power of 5, these quantities are given explicitly by coefficients of eta quotients. In particular:
  - Δ₁₂(625n + 574) + 5Δ₁₂(25n + 24) ≡ 0 (mod 5⁴) whenever n ≢ 2 (mod 5).
- **Corollary 1.4.** The first two families are sharp for every α.

The proof needs no computer. The only outside input is Smoot's fifth-degree modular equation for a Hauptmodul of Γ₀(10). All remaining finite checks are small integers, listed in Appendix A of the paper.

## Contents

| Path | Description |
|---|---|
| `paper/broken12_arxiv.tex` | LaTeX source (single file, standard packages) |
| `paper/broken12_arxiv.pdf` | Compiled paper |
| `code/check_tables.py` | Recomputes every table in Appendix A and repeats each finite check used in the proofs |

## Running the checks

The script needs Python 3 and uses only the standard library.

```sh
python3 code/check_tables.py            # tables and finite checks (a few seconds)
python3 code/check_tables.py --qseries  # also numerical sanity checks (about 5 minutes)
```

- **Default mode.** The script starts from the modular equation alone. It rebuilds the polynomials b_j and a_j, the residue tables mod 5⁶ and the reduced 2×2 and 3×3 operator data, then checks each finite statement used in Sections 5–8.
- **`--qseries` mode.** This adds direct q-expansion checks of the level-10 identities and the modular equation. It also tests Theorems 1.1–1.3 against values of Δ₁₂(n) computed from the generating function f₂f₂₅/(f₁³f₅₀). These tests are sanity checks only; the paper does not use them.

## Building the paper

```sh
cd paper && pdflatex broken12_arxiv.tex && pdflatex broken12_arxiv.tex
```

## Citation

If you use this code, please cite the archived version:

> Ramanuja, *Congruences modulo powers of 5 for broken 12-diamond partitions* (code and manuscript), Zenodo, 2026. https://doi.org/10.5281/zenodo.24521408

## Contact

Ramanuja, ramanujamath@yahoo.com
