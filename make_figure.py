#!/usr/bin/env python3
"""Figure for the autoconvolution-certificates note.

Three certified step constructions (display normalized to integral 1) and
their autoconvolutions. Data: fig_data_ac.json (min/max envelope downsample
of the exact leader vectors; display is floating point -- every certified
quantity in the note is exact)."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = json.load(open("fig_data_ac.json"))

TITLES = {
    "board2": r"$C_1$ construction ($n=90000$, Hyra)",
    "board3": r"$C_2$ construction ($n=524288$, Hyra)",
    "board4": r"$C_3$ construction ($n=100000$, OrganonAgent)",
}
QLBL = {
    "board2": r"$\sup_t g = Q_1 \approx 1.5028503$",
    "board3": r"$\|g\|_2^2/(\|g\|_1\|g\|_\infty) = Q_2 \approx 0.9629011$",
    "board4": r"$\max_t g = Q_3 \approx 1.4523043$",
}

fig, axes = plt.subplots(2, 3, figsize=(11.6, 5.6))
for j, key in enumerate(["board2", "board3", "board4"]):
    b = d[key]
    ax = axes[0][j]
    ax.fill_between(b["x"], b["flo"], b["fhi"], color="#33507a", lw=0.4)
    ax.plot(b["x"], b["fhi"], color="#33507a", lw=0.3)
    ax.set_xlim(-0.25, 0.25)
    ax.set_title(TITLES[key], fontsize=9.5)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_ylabel(r"$f$" if j == 0 else "")
    ax.tick_params(labelsize=8)

    ax = axes[1][j]
    ax.fill_between(b["xg"], b["glo"], b["ghi"], color="#7a3333", lw=0.4)
    ax.plot(b["xg"], b["ghi"], color="#7a3333", lw=0.3)
    ax.set_xlim(-0.5, 0.5)
    ax.axhline(0, color="k", lw=0.5)
    if key != "board3":
        ax.plot([b["xstar"]], [b["gmax"]], marker="o", ms=4, mfc="none",
                mec="k", mew=0.9)
    ax.set_ylabel(r"$g=f*f$" if j == 0 else "")
    ax.set_xlabel(QLBL[key], fontsize=9)
    ax.tick_params(labelsize=8)

fig.tight_layout(pad=0.8)
fig.savefig("fig_ac_constructions.pdf")
print("wrote fig_ac_constructions.pdf")
