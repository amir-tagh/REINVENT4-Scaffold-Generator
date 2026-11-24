#!/usr/bin/env python

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


COLORS = {
    "Full library": "#9CA3AF",
    "CMOS V1": "#4477AA",
    "CMOS V2": "#228833"
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "axes.linewidth": 1.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.facecolor": "white"
})


def save_figure(fig, output, name):

    for extension in ["pdf", "svg", "png"]:

        path = output / f"{name}.{extension}"

        fig.savefig(path, dpi=600)

    plt.close(fig)

    print("Saved:", name)


def boxplot(ax, groups, column, ylabel):

    values = [
        pd.to_numeric(
            df[column], errors="coerce"
        ).dropna().values
        for df in groups.values()
    ]

    bp = ax.boxplot(
        values,
        patch_artist=True,
        showfliers=False,
        widths=0.6,
        medianprops={
            "color": "black",
            "linewidth": 1.4
        }
    )

    for patch, name in zip(
        bp["boxes"], groups
    ):
        patch.set_facecolor(COLORS[name])
        patch.set_alpha(0.75)

    ax.set_xticks(
        range(1, len(groups) + 1)
    )

    ax.set_xticklabels(
        list(groups),
        rotation=20,
        ha="right"
    )

    ax.set_ylabel(ylabel)


def panel_label(ax, letter):

    ax.text(
        -0.12, 1.05,
        letter,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold"
    )


def molecular_properties(groups, output):

    properties = [
        ("QED", "QED"),
        ("SA", "Synthetic accessibility"),
        ("SCScore", "Synthetic complexity"),
        ("MW", "Molecular weight (Da)"),
        ("cLogP", "cLogP"),
        ("TPSA", "TPSA ($\\AA^2$)"),
        ("HBD", "Hydrogen bond donors"),
        ("HBA", "Hydrogen bond acceptors"),
        ("ParentSimilarity", "Parent similarity")
    ]

    fig, axes = plt.subplots(
        3, 3,
        figsize=(11, 10),
        constrained_layout=True
    )

    for i, (column, label) in enumerate(
        properties
    ):

        ax = axes.flat[i]

        boxplot(ax, groups, column, label)

        panel_label(ax, chr(65 + i))

    save_figure(
        fig, output,
        "Figure_01_Molecular_Properties"
    )


def solubility(groups, output):

    properties = [
        ("LogS_ESOL", "ESOL"),
        ("LogS_ADMETAI", "ADMET-AI"),
        ("LogS_Consensus", "Consensus")
    ]

    fig, axes = plt.subplots(
        1, 3,
        figsize=(12, 4.5),
        constrained_layout=True
    )

    for i, (column, title) in enumerate(
        properties
    ):

        ax = axes[i]

        boxplot(
            ax, groups,
            column,
            "Predicted logS (mol/L)"
        )

        ax.set_title(title)

        ax.axhline(
            np.log10(100e-6),
            color="#555555",
            linestyle="--",
            linewidth=1,
            label="100 µM"
        )

        ax.axhline(
            np.log10(300e-6),
            color="#AA3377",
            linestyle=":",
            linewidth=1.4,
            label="300 µM"
        )

        panel_label(ax, chr(65 + i))

    axes[-1].legend(loc="lower right")

    save_figure(
        fig, output,
        "Figure_02_Solubility"
    )


def multiparameter_performance(groups, output):

    properties = [
        ("QED", "QED"),
        ("Score_CNS_Partial", "CNS"),
        ("Score_BBB", "BBB"),
        ("Score_ADME", "ADME"),
        ("Score_Toxicity", "Toxicity")
    ]

    fig, ax = plt.subplots(
        figsize=(8, 4.5)
    )

    x = np.arange(len(properties))
    width = 0.24

    for i, (name, df) in enumerate(
        groups.items()
    ):

        means = [
            df[column].mean()
            for column, _ in properties
        ]

        ax.bar(
            x + (i - 1) * width,
            means,
            width,
            label=name,
            color=COLORS[name]
        )

    ax.set_xticks(x)

    ax.set_xticklabels([
        label for _, label in properties
    ])

    ax.set_ylabel("Mean desirability score")
    ax.set_ylim(0, 1.05)
    ax.legend(ncol=3, loc="upper center")

    save_figure(
        fig, output,
        "Figure_03_Multiparameter_Performance"
    )


def toxicity(groups, output):

    fig, axes = plt.subplots(
        1, 2,
        figsize=(10, 4.5),
        constrained_layout=True
    )

    boxplot(
        axes[0],
        groups,
        "Score_Toxicity",
        "Toxicity desirability"
    )

    panel_label(axes[0], "A")

    x = np.arange(6)
    width = 0.24

    for i, (name, df) in enumerate(
        groups.items()
    ):

        counts = (
            df["Toxicity_Alert_Count"]
            .value_counts(normalize=True)
            .reindex(x, fill_value=0)
            * 100
        )

        axes[1].bar(
            x + (i - 1) * width,
            counts.values,
            width,
            label=name,
            color=COLORS[name]
        )

    axes[1].set_xticks(x)

    axes[1].set_xlabel(
        "Number of predicted toxicity alerts"
    )

    axes[1].set_ylabel("Molecules (%)")
    axes[1].legend()

    panel_label(axes[1], "B")

    save_figure(
        fig, output,
        "Figure_04_Toxicity"
    )


def structural_diversity(benchmark, output):

    df = pd.read_csv(
        benchmark / "diversity_comparison.csv"
    ).set_index("Group")

    properties = [
        (
            "Unique_Murcko_Scaffolds",
            "Unique Murcko scaffolds"
        ),
        (
            "Mean_Pairwise_Tanimoto",
            "Mean pairwise Tanimoto"
        ),
        (
            "Mean_Nearest_Neighbor_Tanimoto",
            "Mean nearest-neighbor Tanimoto"
        )
    ]

    fig, axes = plt.subplots(
        1, 3,
        figsize=(11, 4),
        constrained_layout=True
    )

    names = ["CMOS V1", "CMOS V2"]

    for i, (column, ylabel) in enumerate(
        properties
    ):

        ax = axes[i]

        values = df.loc[
            ["CMOS_V1", "CMOS_V2"],
            column
        ].values

        bars = ax.bar(
            names,
            values,
            color=[COLORS[n] for n in names],
            width=0.6
        )

        labels = [
            f"{v:.0f}" if i == 0
            else f"{v:.4f}"
            for v in values
        ]

        ax.bar_label(
            bars,
            labels=labels,
            padding=3
        )

        ax.set_ylim(
            0, max(values) * 1.25
        )

        ax.set_ylabel(ylabel)

        panel_label(ax, chr(65 + i))

    save_figure(
        fig, output,
        "Figure_05_Structural_Diversity"
    )


def selection_overlap(benchmark, output):

    df = pd.read_csv(
        benchmark / "selection_overlap.csv"
    )

    df = df[df["Parent"] != "ALL"]

    fig, ax = plt.subplots(
        figsize=(9, 4.5)
    )

    x = np.arange(len(df))

    shared = df["Shared"].values
    exclusive = df["V2_Exclusive"].values

    ax.bar(
        x,
        shared,
        color="#888888",
        label="Shared with V1"
    )

    ax.bar(
        x,
        exclusive,
        bottom=shared,
        color=COLORS["CMOS V2"],
        label="V2 exclusive"
    )

    for i, (s, e) in enumerate(
        zip(shared, exclusive)
    ):

        ax.text(
            i, s / 2, str(s),
            ha="center",
            va="center",
            color="white"
        )

        ax.text(
            i, s + e / 2, str(e),
            ha="center",
            va="center",
            color="white"
        )

    ax.set_xticks(x)

    ax.set_xticklabels(
        df["Parent"].str.replace(
            "Hit_", "Hit ", regex=False
        ),
        rotation=30,
        ha="right"
    )

    ax.set_ylabel("Selected molecules")
    ax.set_ylim(0, 35)
    ax.legend(ncol=2)

    save_figure(
        fig, output,
        "Figure_06_Selection_Overlap"
    )


def sensitivity_analysis(path, output):

    df = pd.read_csv(path)

    fig, axes = plt.subplots(
        1, 2,
        figsize=(11, 4.5),
        constrained_layout=True
    )

    for _, row in df.iterrows():

        name = row["Profile"].replace(
            "cns_", ""
        ).capitalize()

        color = (
            "#228833"
            if name == "Safety"
            else "#4477AA"
        )

        axes[0].scatter(
            row["Solubility_uM"],
            row["Toxicity"],
            s=75,
            color=color
        )

        axes[0].annotate(
            name,
            (
                row["Solubility_uM"],
                row["Toxicity"]
            ),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=9
        )

        axes[1].scatter(
            row["Pairwise_Diversity"],
            row["Scaffolds"],
            s=75,
            color=color
        )

        axes[1].annotate(
            name,
            (
                row["Pairwise_Diversity"],
                row["Scaffolds"]
            ),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=9
        )

    axes[0].set_xlabel(
        "Consensus solubility (µM)"
    )

    axes[0].set_ylabel(
        "Toxicity desirability"
    )

    axes[1].set_xlabel(
        "Mean pairwise diversity"
    )

    axes[1].set_ylabel(
        "Unique Murcko scaffolds"
    )

    axes[1].yaxis.set_major_locator(
        MaxNLocator(integer=True)
    )

    for i, ax in enumerate(axes):

        ax.margins(x=0.20, y=0.15)

        panel_label(ax, chr(65 + i))

    save_figure(
        fig, output,
        "Figure_07_Weight_Sensitivity"
    )


def parent_enrichment(benchmark, output):

    df = pd.read_csv(
        benchmark / "per_parent_comparison.csv"
    )

    properties = [
        ("QED", "QED"),
        ("LogS_Consensus", "Consensus logS"),
        ("Score_CNS_Partial", "CNS"),
        ("Score_ADME", "ADME"),
        ("Score_Toxicity", "Toxicity"),
        ("Toxicity_Alert_Count", "Toxicity alerts")
    ]

    fig, axes = plt.subplots(
        2, 3,
        figsize=(13, 7),
        constrained_layout=True
    )

    for i, (column, label) in enumerate(
        properties
    ):

        ax = axes.flat[i]

        sub = df[
            df["Property"].eq(column)
        ].sort_values("ParentID")

        values = sub["V2_vs_V1"].values

        favorable = (
            values < 0
            if column == "Toxicity_Alert_Count"
            else values > 0
        )

        colors = [
            "#228833" if good else "#CC6677"
            for good in favorable
        ]

        ax.bar(
            np.arange(len(sub)),
            values,
            color=colors
        )

        ax.axhline(
            0,
            color="black",
            linewidth=0.8
        )

        ax.set_xticks(
            np.arange(len(sub))
        )

        ax.set_xticklabels(
            sub["ParentID"].str.replace(
                "Hit_", "", regex=False
            ),
            rotation=45
        )

        ax.set_ylabel("V2 − V1")
        ax.set_title(label)

        panel_label(ax, chr(65 + i))

    save_figure(
        fig, output,
        "Figure_08_Parent_Enrichment"
    )


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--benchmark",
        required=True
    )

    parser.add_argument(
        "--sensitivity"
    )

    parser.add_argument(
        "--output",
        required=True
    )

    args = parser.parse_args()

    benchmark = Path(args.benchmark)
    output = Path(args.output)

    output.mkdir(
        parents=True,
        exist_ok=True
    )

    groups = {
        "Full library": pd.read_csv(
            benchmark /
            "library_selection_membership.csv"
        ),
        "CMOS V1": pd.read_csv(
            benchmark /
            "CMOS_V1_Selected_Annotated.csv"
        ),
        "CMOS V2": pd.read_csv(
            benchmark /
            "CMOS_V2_Selected_Annotated.csv"
        )
    }

    print("\nInput libraries:")

    for name, df in groups.items():

        print(f"{name}: {len(df)}")

    molecular_properties(groups, output)
    solubility(groups, output)
    multiparameter_performance(groups, output)
    toxicity(groups, output)
    structural_diversity(benchmark, output)
    selection_overlap(benchmark, output)

    if args.sensitivity:

        sensitivity_analysis(
            Path(args.sensitivity),
            output
        )

    parent_enrichment(benchmark, output)

    print("\nFIGURE GENERATION COMPLETE")
    print("Output:", output)


if __name__ == "__main__":
    main()
