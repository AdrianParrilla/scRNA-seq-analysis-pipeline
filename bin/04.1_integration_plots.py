#!/opt/env/bin/python

import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text


def integration_dot_plot(metrics):

    metrics_sub = metrics.loc[:,['Embedding','Batch correction','Bio conservation']]
    metrics_sub = metrics_sub[metrics_sub['Embedding'] != 'Metric Type']

    metrics_sub["Batch correction"] = metrics_sub["Batch correction"].astype(float)
    metrics_sub["Bio conservation"] = metrics_sub["Bio conservation"].astype(float)
    metrics_sub["Embedding"] = metrics_sub["Embedding"].astype(str)

    plt.figure(figsize=(6, 5)) 

    ax = sns.scatterplot(
        data=metrics_sub,
        x="Batch correction",
        y="Bio conservation",
        hue="Embedding",
        palette="husl",
        s=100,
        legend=False
    )

    texts = []

    for _, row in metrics_sub.iterrows():
        texts.append(plt.text(
            row["Batch correction"], 
            row["Bio conservation"], 
            row["Embedding"], 
            fontsize=14
        ))

    adjust_text(
        texts, 
        expand=(2, 2),
        arrowprops=dict(arrowstyle='->', color='grey', lw=0.5)
    )

    ax.grid(False)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_xlabel("Batch correction", fontsize=16)
    ax.set_ylabel("Bio conservation", fontsize=16)

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

    plt.savefig("integration_dot_plot.png", dpi=300, bbox_inches="tight")
    
def radar_plot(metrics):

    metrics_sub = metrics.loc[metrics['Embedding'] != 'Metric Type',
        ~metrics.columns.isin(['Batch correction', 'Bio conservation', 'Total'])]

    scores = metrics_sub.columns[1:]
    metrics_sub[scores] = metrics_sub[scores].apply(pd.to_numeric)

    print(scores)
    N = len(scores)

    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]

    fig = plt.figure(figsize=(8, 8))
    ax = plt.subplot(111, polar=True)

    bio_metrics = {"Isolated labels", "KMeans NMI", "KMeans ARI", "Silhouette label", "cLISI"}
    batch_metrics = {"BRAS", "iLISI", "Graph connectivity", "PCR comparison"}

    label_colors = []
    for m in scores:
        if m in bio_metrics:
            label_colors.append("#2e7d32")
        elif m in batch_metrics:
            label_colors.append("#1565c0")
        else:
            label_colors.append("black")

    methods = metrics_sub["Embedding"].unique()
    husl_palette = sns.color_palette("husl", len(methods))
    method_colors = dict(zip(methods, husl_palette))

    for _, row in metrics_sub.iterrows():
        method = row["Embedding"]
        values = row[scores].tolist()
        values += values[:1]
        color = method_colors.get(method, "gray")
        ax.plot(angles, values, linewidth=2, alpha=0.8, label=method, color=color)
        ax.fill(angles, values, alpha=0.02, color=color)

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(scores, fontsize=12, fontweight='bold')
    for tick_label, color in zip(ax.get_xticklabels(), label_colors):
        tick_label.set_color(color)

    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"], fontsize=10)
    ax.set_ylim(0, 1.05)

    leg = ax.legend(loc="upper left", bbox_to_anchor=(1.05, 1.05), frameon=False)
    for text in leg.get_texts():
        text.set_fontsize(12)

    plt.savefig("radar_plot.png", dpi=300, bbox_inches="tight")


def main(metrics_csv):

    metrics = pd.read_csv(metrics_csv)

    integration_dot_plot(metrics)
    radar_plot(metrics)
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Integrate scRNA seq data from different batches using scVI')
    parser.add_argument('-i', '--metrics_csv', type=str, required=True, help='CSV output file from scib-metrics')

    args = parser.parse_args()

    main(
     metrics_csv = args.metrics_csv
     )

