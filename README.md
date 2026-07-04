# Exact-arithmetic certificates for three autoconvolution inequalities

**Kevin Russell** — *ProjectForty2 / CHRONOS agent*

This repository contains a short computer-assisted note and its full
verification package: machine-verified, exact-rational bounds for the three
autoconvolution-inequality constants popularized by the recent AI-search
literature (AlphaEvolve problems B.1–B.3; Green's "100 open problems" §6),
plus exact re-evaluations of four published constructions. Writing
`g = f∗f` for the autoconvolution:

> **Theorem 1** (first autocorrelation inequality). Let `C₁` be the largest
> constant with `max_t (f∗f)(t) ≥ C₁ (∫f)²` for all `f ≥ 0` supported on
> `[-1/4, 1/4]`. Then
>
> ```
> C₁ ≤ Q₁ < 1.5028503020710078765191888 ,
> Q₁ = 1661053907432921000898809384050165075655877436050568807950600634368000
>      ──────────────────────────────────────────────────────────────────────
>      1105269037870172510060979851604850227379926078898766312058428453057449
> ```
>
> **Theorem 2** (second autocorrelation inequality). Let `C₂` be the smallest
> constant with `‖g‖₂² ≤ C₂ ‖g‖₁ ‖g‖∞` for all `f ≥ 0`. Then
>
> ```
> C₂ ≥ Q₂ > 0.9629011010961756912885013 ,
> Q₂ = 140651861665566489683881393353250795846281833
>      ─────────────────────────────────────────────
>      146070932420211259869783468438333325818535926
> ```
>
> **Theorem 3** (third autocorrelation inequality, signed `f`). Let `C₃` be
> the largest constant with `|max_t (f∗f)(t)| ≥ C₃ (∫f)²` for all
> `f: [-1/4,1/4] → ℝ` (sign changes allowed) with `∫f ≠ 0`. Then
>
> ```
> C₃ ≤ Q₃ < 1.452304333183157960701885 ,
> Q₃ = 11753128449293701953238517385067272445617294540800000
>      ─────────────────────────────────────────────────────
>      8092744874989952471246071559466128309374865340943729
> ```

All three rationals are exact evaluations — Python integer / `Fraction`
arithmetic, Kronecker-substitution convolution up to `n = 524288` cells,
**no floating-point operation on any certified path** — of explicit step
constructions, via a short lemma (piecewise linearity of `g` for step `f`)
that makes the continuum functional a finite exact computation with zero
discretization loss.

Context, per constant:

| | our certificate | previous status |
|---|---|---|
| `C₁` | `≤ 1.50285031` | best *proven*: `π/2 = 1.5707963…` (Schinzel–Schmidt); every sub-1.51 value in the literature is float-only |
| `C₂` | `≥ 0.96290110` | best prior record `0.94136` (Jaech–Joseph) — float-only in-paper, **upgraded to proven here** by exact re-certification of their published vector |
| `C₃` | `≤ 1.45230434` | no rigorously certified nontrivial bound existed in either direction, to our knowledge |

## Attribution

Credit belongs where the mathematics originated.

- **The three certified constructions** were produced by anonymous AI search
  agents competing on the EinsteinArena platform — **"Hyra"** (Theorems 1
  and 2) and **"OrganonAgent"** (Theorem 3) — within an open ecosystem of
  agents iteratively optimizing these functionals. We make no priority claim
  over any construction. The certificates evaluate the *existing* leader
  vectors: they tie, and do not beat, the corresponding floating-point
  leaderboard scores (agreement to roughly one unit in the last place); the
  contribution is that float records become theorems.
- **Re-certified published constructions:** Matolcsi–Vinuesa (2010, `n=150`
  signed example → exactly `1.4581056857680625`), AlphaEvolve (2025, `n=400`
  vector → exactly `1.4556427953745406`), Boyer–Li (2025, 575-integer vector
  → their printed exact fraction, reproduced to the last digit), and
  Jaech–Joseph (2025, 2540-interval vector → exactly
  `0.9413628949737643… ≥ 0.94136`).
- **Method lineage:** the certification lemma + exact-rational evaluation
  recipe follows our note on the Erdős minimum-overlap constant
  ([erdos-minimum-overlap-bound](https://github.com/techno-optimist/erdos-minimum-overlap-bound)).

This note was prepared computer-assisted with **CHRONOS**, ProjectForty2's
autonomous research agent, under the author's direction; the author reviewed
the mathematics and takes responsibility for all claims.

## What this package does *not* prove

Honest scoping is the point of the note, so it is worth repeating here:

- **Every bound is one-sided.** Theorems 1 and 3 are *upper* bounds only
  (the lower sides stand as published: `C₁ ≥ 1.28` Cloninger–Steinerberger,
  and the trivial `C₃ ≥ 1`); Theorem 2 is a *lower* bound only (on the upper
  side nothing beyond the trivial Hölder `C₂ ≤ 1` is known).
- **No leaderboard is beaten.** The certificates evaluate the existing arena
  leader vectors exactly; they tie the float scores to ~1 ulp.
- **Priority claims are hedged.** "First exact-certified upper bound in the
  sub-1.51 regime" (Theorem 1) and "first rigorously certified bound for the
  signed variant" (Theorem 3) hold *to our knowledge*, per the targeted
  literature survey in the note; tabulated 2026 values of a Together AI team
  have unverified certification status (and our values are on the favorable
  side of them regardless).
- **The `C₃` cousin constant is different.** Values `1.4993 / 1.4688` in the
  literature refer to `max_t |f∗f(t)|` (absolute value *inside* the max) —
  a different problem; do not compare.
- See the pending-verification ledger (Section 7 of the note) for the full
  list of items not machine-verified at the time of writing.

## Reproduce it

Requires **Python ≥ 3.9**, standard library only, for every certified path.
(`numpy` is needed only by `exact_board4.py`'s *non-certified* float
cross-check section and the figure; `pip install -r requirements.txt`.)
All commands run from the repository root; input vectors ship in `vectors/`.

```bash
# Theorem 1 — exact certification of the n=90000 vector (~6 s on Apple silicon)
python3 scripts/exact_board2.py | grep -E "^C_exact_(num|den) " \
    | diff - <(grep -E "^C_exact_(num|den) " certs/exact_board2_out.txt) && echo OK1

# Theorem 2 — n=524288, big-integer squaring dominates (~20 s)
python3 scripts/exact_board3.py | grep -E "^(S1 |Linf |S2 |exact_num|exact_den)" \
    | diff - <(grep -E "^(S1 |Linf |S2 |exact_num|exact_den)" certs/exact_board3_out.txt) && echo OK2

# Theorem 3 — signed Kronecker split, n=100000 (~15 s; writes exact_board4_out.txt to cwd)
python3 scripts/exact_board4.py >/dev/null && diff \
    <(grep -E "^(num|den)=" exact_board4_out.txt) \
    <(grep -E "^(num|den)=" certs/exact_board4_out.txt) && echo OK3

# Re-certifications: Boyer–Li + Jaech–Joseph (seconds)
python3 scripts/certify.py        # expect: "EXACT MATCH to published fraction? True"
                                  #         ">= 0.94136 (their claim) ? True"

# Re-certifications: Matolcsi–Vinuesa n=150 + AlphaEvolve n=400 (seconds)
python3 scripts/cert_published_c3.py   # expect floats 1.4581056857680625 / 1.4556427953745406
```

The cert files in `certs/` include wall-clock lines from the original runs,
so the checks above compare the invariant lines (the exact numerators and
denominators) rather than byte-diffing whole files. The exact fractions must
match to the last digit; they were produced independently on two machines
(the reference certs on a workstation-class GPU host, the checks above
re-verified on an Apple-silicon Mac before publication).

## Build the PDF

The prebuilt `autocorrelation_certificates.pdf` is included. To rebuild
(e.g. with [tectonic](https://tectonic-typesetting.github.io/)):

```bash
python3 make_figure.py                    # regenerates fig_ac_constructions.pdf
tectonic autocorrelation_certificates.tex
```

## Layout

```
.
├── autocorrelation_certificates.tex   the note (source of truth)
├── autocorrelation_certificates.pdf   prebuilt PDF (13 pp.)
├── fig_ac_constructions.pdf           Figure 1 (three constructions; display-only)
├── make_figure.py + fig_data_ac.json  regenerate Figure 1
├── requirements.txt                   optional deps (figure + one float cross-check)
├── scripts/
│   ├── exact_board2.py                Theorem 1: exact Q₁ (Kronecker, n=90000)
│   ├── exact_board3.py                Theorem 2: exact Q₂ (Kronecker, n=524288)
│   ├── exact_board4.py                Theorem 3: exact Q₃ (signed Kronecker, n=100000)
│   ├── certify.py                     Boyer–Li + Jaech–Joseph re-certifications
│   └── cert_published_c3.py           Matolcsi–Vinuesa + AlphaEvolve re-certifications
├── certs/                             reference outputs (exact fractions, enclosures)
└── vectors/                           all input vectors, frozen
    ├── ac_board2_leader.json          Hyra n=90000 (EinsteinArena board copy)
    ├── ac_board3_leader.json          Hyra n=524288 (EinsteinArena board copy)
    ├── ac_board4_leader.json          OrganonAgent n=100000 (EinsteinArena board copy)
    ├── mv150_printed.txt              Matolcsi–Vinuesa n=150 (arXiv:0907.1379 appendix)
    ├── AE-2025-third-400.json         AlphaEvolve n=400 (public results notebook)
    ├── coeffBL.txt                    Boyer–Li 575 integers (their repository)
    └── jj_function.py                 Jaech–Joseph 2540 floats (their repository)
```

## License

Code and certificate data are released under the MIT License (see `LICENSE`).
The note text and figures (`*.tex`, `*.pdf`) are © 2026 Kevin Russell,
released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Third-party construction vectors in `vectors/` are redistributed with
attribution: the three EinsteinArena leader vectors as public leaderboard
submissions (Hyra, OrganonAgent), and the published vectors from their
respective public sources (Matolcsi–Vinuesa arXiv:0907.1379; AlphaEvolve
results notebook; Boyer–Li and Jaech–Joseph repositories, as cited in the
note).
