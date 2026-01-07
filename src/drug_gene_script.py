import os
import json
import logging
from datetime import datetime

import pandas as pd
import requests
import matplotlib.pyplot as plt

# ──────────────────────────────────────────────
# Configuration (default paths – will be overridden in main())
# ──────────────────────────────────────────────
INPUT_FILE = "output/hub_genes.csv"        # CSV must contain a 'Gene' column
OUTPUT_DIR = "output"
OUTPUT_FULL_CSV = os.path.join(OUTPUT_DIR, "drug_gene_interactions.csv")
OUTPUT_BAR = os.path.join(OUTPUT_DIR, "drug_gene_interaction_barplot.png")

TOP_GENES = 15  # number of top genes to display in stacked bar chart

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ──────────────────────────────────────────────
# Logging
# ──────────────────────────────────────────────
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)-7s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger("DGIdbPipeline")

logger = setup_logging()

# ──────────────────────────────────────────────
# Load Genes
# ──────────────────────────────────────────────
def load_gene_names(file_path):
    try:
        df = pd.read_csv(file_path)
        if "Gene" not in df.columns:
            raise ValueError("CSV is missing a 'Gene' column.")
        genes = df["Gene"].dropna().unique().tolist()
        logger.info("Loaded %d unique genes.", len(genes))
        return genes
    except Exception as e:
        logger.error("Failed to read input file: %s", e)
        return []

# ──────────────────────────────────────────────
# Query DGIdb GraphQL
# ──────────────────────────────────────────────
def query_dgidb(genes):
    url = "https://dgidb.org/api/graphql"
    query = f"""
    {{
      genes(names: {json.dumps(genes)}) {{
        nodes {{
          name
          interactions {{
            drug {{ name conceptId }}
            interactionScore
            interactionTypes {{ type directionality }}
            interactionAttributes {{ name value }}
            publications {{ pmid }}
            sources {{ sourceDbName }}
          }}
        }}
      }}
    }}
    """
    try:
        r = requests.post(url, json={"query": query}, timeout=60)
        if r.status_code != 200:
            logger.error("DGIdb API error (%d): %s", r.status_code, r.text)
            return []
        return r.json().get("data", {}).get("genes", {}).get("nodes", [])
    except Exception as e:
        logger.error("DGIdb request failed: %s", e)
        return []

# ──────────────────────────────────────────────
# Extract Data (flatten)
# ──────────────────────────────────────────────
def extract_interactions(data):
    records = []
    for item in data:
        gene = item.get("name", "")
        for inter in item.get("interactions", []):
            drug_info = inter.get("drug", {})
            records.append({
                "Gene": gene,
                "Drug": drug_info.get("name", ""),
                "Concept_ID": drug_info.get("conceptId", ""),
                "Score": inter.get("interactionScore", 0),
                "Interaction_Types": "; ".join([
                    f"{t.get('type','')} ({t.get('directionality','')})"
                    for t in inter.get("interactionTypes", [])
                ]) or None,
                "Interaction_Attributes": "; ".join([
                    f"{a.get('name','')}={a.get('value','')}"
                    for a in inter.get("interactionAttributes", [])
                ]),
                "Publications_PMIDs": "; ".join([
                    str(p.get("pmid", "")) for p in inter.get("publications", [])
                ]),
                "Sources": "; ".join([
                    s.get("sourceDbName", "") for s in inter.get("sources", [])
                ])
            })
    df = pd.DataFrame(records)
    logger.info("Extracted %d drug–gene interactions.", len(df))
    return df

# ──────────────────────────────────────────────
# Save Full CSV
# ──────────────────────────────────────────────
def export_full_csv(df):
    try:
        df.to_csv(OUTPUT_FULL_CSV, index=False)
        logger.info("Full CSV saved → %s", OUTPUT_FULL_CSV)
    except Exception as e:
        logger.error("Failed to save full CSV: %s", e)

# ──────────────────────────────────────────────
# Stacked Bar Plot
# ──────────────────────────────────────────────
def draw_stacked_bar(df):
    # Ensure copy and fill missing interaction types
    stacked_df = df.copy()
    stacked_df["Interaction_Types"] = stacked_df["Interaction_Types"].fillna("Unknown")

    # Count unique drugs per gene & interaction type
    gene_interaction_counts = (
        stacked_df.groupby(["Gene", "Interaction_Types"])["Drug"]
        .nunique()
        .reset_index()
    )

    # Pivot for stacked bar chart
    pivot_df = (
        gene_interaction_counts.pivot(
            index="Gene", columns="Interaction_Types", values="Drug"
        )
        .fillna(0)
    )

    # Keep only top genes by number of unique drugs
    top_genes = (
        stacked_df.groupby("Gene")["Drug"]
        .nunique()
        .sort_values(ascending=False)
        .head(TOP_GENES)
        .index
    )
    pivot_df = pivot_df.loc[pivot_df.index.intersection(top_genes)]

    # Plot
    ax = pivot_df.plot(
        kind="bar",
        stacked=True,
        figsize=(12, 7),
        colormap="tab20"
    )
    ax.grid(True, which="both", axis="both", linestyle="--", alpha=0.6)
    plt.title("Drugs per Gene by Interaction Type", fontsize=14, weight="bold")
    plt.xlabel("Genes", fontsize=12)
    plt.ylabel("Number of Drugs", fontsize=12)
    plt.xticks(rotation=75, ha="right")
    plt.legend(title="Interaction Type", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(OUTPUT_BAR, dpi=300)
    plt.close()
    logger.info("Stacked bar plot saved → %s", OUTPUT_BAR)

# ──────────────────────────────────────────────
# Pipeline Runner
# ──────────────────────────────────────────────
def run_pipeline():
    start = datetime.now()
    logger.info("=== DGIdb Drug–Gene Interaction Pipeline Started ===")

    genes = load_gene_names(INPUT_FILE)
    if not genes:
        logger.warning("No valid gene names. Exiting.")
        return

    logger.info("Querying DGIdb for %d genes...", len(genes))
    data = query_dgidb(genes)
    if not data:
        logger.warning("No interaction data retrieved.")
        return

    df = extract_interactions(data)
    if df.empty:
        logger.warning("No interactions to process.")
        return

    export_full_csv(df)
    draw_stacked_bar(df)

    logger.info("Pipeline completed in %s.", datetime.now() - start)

# ──────────────────────────────────────────────
# Main callable for pipeline integration
# ──────────────────────────────────────────────
def main(hub_genes_path="hub_genes.csv", output_dir="output"):
    """
    Run DGIdb drug–gene interaction pipeline.
    Parameters
    ----------
    hub_genes_path : str
        Path to CSV file containing a 'Gene' column (hub genes).
    output_dir : str
        Directory where results (CSV + plot) will be saved.
    """
    global INPUT_FILE, OUTPUT_DIR, OUTPUT_FULL_CSV, OUTPUT_BAR
    INPUT_FILE = hub_genes_path
    OUTPUT_DIR = output_dir
    OUTPUT_FULL_CSV = os.path.join(OUTPUT_DIR, "drug_gene_interactions.csv")
    OUTPUT_BAR = os.path.join(OUTPUT_DIR, "drug_gene_interaction_barplot.png")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    run_pipeline()