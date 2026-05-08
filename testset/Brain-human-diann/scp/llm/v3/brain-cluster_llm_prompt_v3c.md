# Cluster marker summary for cell type annotation

## Experiment context

- **Experiment**: BrainSCPpx
- **Species**: human
- **Description**: fresh prenatal human brain samples at gestational weeks (GWs) 13, 15 and 19. For GW15 and GW19, we microdissected the germinal zone (GZ) and cortical plate (CP) from the prefrontal cortex. For the GW13 sample, we used the entire nascent frontal telencephalon (including ganglionic eminences

## Run metadata (auto-populated at runtime)

- **Prompt version**: v3c
- **Model**: {model_name}
- **Temperature**: 0
- **Run timestamp**: {iso_timestamp}
- **Prompt hash (SHA256)**: {sha256_of_prompt_text}

*These fields are populated by the calling script and stored alongside the
raw LLM output for audit trail purposes.*

## Developmental stage constraints

This experiment uses **prenatal human brain tissue at GW13–19**. The following stage-specific rules apply and take priority over generic cell-type vocabulary:

- **Mature astrocytes (GFAP+/AQP4+/ALDH1L1+/GJA1+) do NOT exist at GW13–19.** Astrogliogenesis begins ~GW20–24 and peaks postnatally. For any S100B+ or SIRT2+ population, use **"astroglial progenitor"** or **"glial progenitor"** — NOT "astrocyte".
- **S100B at this stage** marks astroglial progenitors or early OPC-lineage cells, not mature astrocytes. Mature astrocyte identity requires co-detection of AQP4, ALDH1L1, and GJA1; if any of these are absent, do NOT call "astrocyte".
- **DCX co-detection** in any S100B+ cluster is inconsistent with mature astrocyte identity at any developmental stage — DCX marks neuroblasts/immature neurons.
- **GFAP depletion relative to background** in a nominally astrocytic cluster is a contra-indicator for astrocyte identity at this stage.
- **Cycling progenitors** (PCNA+, TOP2A+, MCM-family+) should be labelled "cycling progenitor" plus the lineage qualifier (e.g., "cycling neural progenitor") — not a differentiated cell type.
- **Outer radial glia (oRG)**: HOPX+, TNC+, PTPRZ1+ — these cells are abundant at this stage and should be labelled "outer radial glia" or "oRG".
- **Do not use mature differentiated cell type labels** (neuron, astrocyte, oligodendrocyte) without ≥2 mature lineage markers AND absence of progenitor markers (PCNA, TOP2A, SOX2, NES).

You are a single-cell proteomics expert. Below is a per-cluster summary from a label-free single-cell proteomics (SCP) experiment. For each cluster you are given statistical evidence across three complementary modalities: detection specificity (Fisher exact), intensity fold-change (Wilcoxon, detected-only), and scplainer (linear mixed-effects, correcting for technical variation). Consensus markers are ranked across all three modalities.

**Your task for each cluster:**

1. Propose the most likely **cell type label**, drawing on the marker evidence and the experiment context above.

2. List the **key supporting markers** -- prioritise proteins with concordant evidence across >=2 modalities and/or strong detection specificity (det_rate_in >> det_rate_out, low q-value).

3. List any **contradictory or confounding signals** that weaken the call (e.g. immune markers in a nominally stromal cluster, low n_detected, absence of expected canonical markers).

4. Assign a **structured confidence score** using the rubric below.

5. If confidence is Low or if the cluster is ambiguous, explicitly state what additional marker evidence would resolve the annotation, citing specific proteins by name.

**Confidence rubric (use exactly these labels):**

| Label | Criteria |
|---|---|
| High | >=3 canonical lineage markers concordant across >=2 modalities; no major contradictions |
| Medium-high | 2-3 canonical markers concordant; minor contradictions explainable by biology |
| Medium | 1-2 canonical markers; contradictions present but cell type still most parsimonious |
| Low | Canonical markers absent or contradicted; cluster is ambiguous or likely a doublet/artefact |

**Output format -- return a Markdown table with exactly these columns:**

| Cluster | Suggested cell type | Key supporting markers | Contradictions / caveats | Confidence | Resolving markers (if Low/ambiguous) |

Then, after the table, provide a **per-cluster narrative** (one paragraph each) with quantitative citations (e.g. "Cd68 det_rate_in=1.00 vs 0.25, q=0.0000"). Do not assert a cell type without citing at least one specific statistic.

**Important caveats for this dataset:**

- Proteins detected across >=90% of all clusters should be treated as constitutively expressed or background and are not lineage-informative on their own.
- Depleted detection markers (proteins specifically *absent* in a cluster) are biologically informative and should be used to exclude lineages, not only positive markers to assign them.

**Contamination vs ambient carryover -- be conservative:**

- In single-cell proteomics, keratins (KRT5, KRT14, KRT10, TGM1), serum proteins (ALB, TF, HP, A2M), and hemoglobin (HBB, HBA1) are commonly detected as **ambient carryover** from sample preparation rather than true contaminant cell populations.
- **Do NOT label a cluster as 'Contaminant' solely because it expresses keratins, ALB, or HBB.** These proteins are frequently detected in genuine cell populations alongside their true lineage markers.
- Only call contamination when ALL of: (a) the cluster lacks any coherent lineage-specific markers across all modalities, AND (b) contaminant markers are overwhelmingly dominant (>80% of top markers), AND (c) no alternative biological explanation is plausible given the experiment context.
- When contaminant-associated proteins are present alongside genuine lineage markers, annotate the cluster by its **biological identity** and note ambient contamination as a caveat in the Contradictions column.
- For ambiguous clusters with mixed contaminant and lineage signals, prefer 'Low confidence' biological annotation over a 'Contaminant' label. The researcher can decide post-hoc whether to filter these clusters.

## Mechanistic annotation rules

- **Do NOT assign a mechanistic qualifier** (e.g., "IFN-stressed", "IFN-leaning", "hypoxic", "senescent", "apoptotic") based on a **single marker at <10% detection rate**.
- Mechanistic labels require **≥2 concordant markers from the same pathway** detected at ≥10% in the cluster.
  - "IFN-stressed" or "IFN-leaning" requires ≥2 of {STAT1, ISG15, IFIT1, MX1, OAS1, DDX58} detected at ≥10%. **A single ISG marker (e.g., DDX58 alone at 6%) does NOT qualify.**
  - "hypoxic" requires ≥2 of {HIF1A, LDHA, SLC2A1, VEGFA, ENO1↑} at ≥10%.
- **CRITICAL — label prohibition**: If the ≥2 concordant marker threshold is NOT met, the pathway name (e.g., "IFN", "interferon", "NF-κB", "hypoxic", "senescent") **MUST NOT appear anywhere in the cell type label string**. Place the single-marker observation in the Contradictions/caveats column ONLY.
- If only **1 pathway marker** is detected, describe the phenotype **factually** in the label: "metabolically depleted", "biosynthetically inactive", or "low-complexity proteome" — NOT a mechanistic pathway label.
- When a cluster shows coordinated **depletion** of biosynthetic machinery (ribosomal proteins, mitochondrial complex subunits, replication factors such as MCM proteins, PCNA) combined with low median_n_detected, the primary explanation is **dissociation/handling artefact**. The cell type label MUST include one of: "metabolically depleted", "biosynthetically inactive", "probable dissociation artefact", or "low-complexity". Example: "Cycling neural progenitor (metabolically depleted / probable dissociation artefact)". Do NOT describe such a cluster as "IFN-leaning", "stressed", or any pathway-active state.

## Marker coverage score

For each cluster, the calling script will pre-compute a **marker coverage score** defined as:

    coverage = (number of recommended panel markers detected at det_rate >= 0.20
                in this cluster) / (total recommended panel markers for the
                assigned cell type)

This score (0-1) will be appended to the annotation table post-hoc by the pipeline. You do not need to compute it, but you **must** assign a cell type label that maps unambiguously to one of the panels defined in the recommended markers output so the script can look it up. If the best label is ambiguous between two panels, name both separated by " / " (e.g. "Macrophage (Mrc1+) / Inflammatory myeloid").

## Pass structure

This is **Pass 1**. Your output will be reviewed and a supplemental marker query will be generated for any cluster where:
- Confidence is Low or Medium, or
- A contradictory signal was flagged, or
- The suggested cell type requires confirmation of a specific marker not yet in the data.

Pass 2 will re-present the same cluster data with additional researcher-requested marker statistics appended. In Pass 2, update only clusters where new evidence changes the call or confidence; carry forward unchanged calls explicitly stating "No change from Pass 1."

## Data columns provided per cluster

- **n_cells**: number of single cells in the cluster
- **median_n_detected**: median number of proteins detected per cell
- **detection markers**: proteins most specifically *detected* in this cluster (Fisher exact test on presence/absence). Columns: det_rate_in, det_rate_out, qvalue
- **depleted detection markers**: proteins specifically *absent/under-detected* in this cluster vs others
- **intensity markers (up-regulated)**: proteins with highest *intensity* in detected cells (Wilcoxon, detected-only). Columns: log2FC (positive), qvalue
- **depleted intensity markers (down-regulated)**: proteins with *lowest* intensity in this cluster vs others (negative log2FC)
- **scplainer markers (up/down-regulated)**: linear mixed-effects model, accounting for technical variation. Columns: log2FC, qvalue
- **consensus markers**: proteins ranked consistently across all three modalities. Columns: score, n_modalities (how many of 3 agreed)
- **AUCell pathways**: most active MSigDB Hallmark pathways. Columns: AUC_in, delta_vs_rest
- **researcher markers of interest**: biologically relevant markers provided by the researcher with per-cluster stats (may not reach significance everywhere)

---

## Cluster C0

- **n_cells**: 61
- **median_n_detected**: 636.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | UGDH | 0.64 | 3/3 |
  | S100B | 0.60 | 3/3 |
  | CKB | 0.58 | 3/3 |
  | EZR | 0.37 | 3/3 |
  | TAGLN2 | 0.35 | 3/3 |
  | CLIC1 | 0.33 | 3/3 |
  | VIM | 0.31 | 3/3 |
  | MSI1 | 0.24 | 1/3 |
  | FLNA | 0.23 | 3/3 |
  | SIRT2 | 0.22 | 2/3 |
  | SELENBP1 | 0.22 | 3/3 |
  | CNN3 | 0.22 | 3/3 |
  | PGM3 | 0.22 | 3/3 |
  | CTSD | 0.21 | 3/3 |
  | HSPB1 | 0.20 | 3/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | SRM | 0.97 | 0.73 | 6.3e-05 |
  | PPP1R14B | 0.97 | 0.70 | 2.0e-05 |
  | SNRPD1 | 0.97 | 0.84 | 3.4e-02 |
  | JPT1 | 0.95 | 0.81 | 2.5e-02 |
  | RPL13 | 0.95 | 1.00 | 4.6e-02 |
  | PTGES3 | 0.95 | 1.00 | 4.6e-02 |
  | RPS25 | 0.95 | 1.00 | 5.0e-03 |
  | KPNB1 | 0.93 | 0.99 | 1.5e-02 |
  | YWHAH | 0.93 | 1.00 | 2.0e-03 |
  | EEF1D | 0.93 | 1.00 | 9.0e-03 |
  | COTL1 | 0.93 | 0.52 | 8.4e-10 |
  | ENO2 | 0.90 | 0.72 | 1.4e-02 |
  | CALR | 0.89 | 0.99 | 2.2e-03 |
  | DCX | 0.89 | 0.63 | 4.1e-04 |
  | MAPK1 | 0.89 | 0.99 | 2.2e-03 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | TAGLN2 | 0.08 | 0.90 | 4.8e-40 |
  | CLIC1 | 0.16 | 0.94 | 3.8e-41 |
  | FLNA | 0.08 | 0.82 | 3.1e-30 |
  | MSI1 | 0.02 | 0.73 | 7.9e-28 |
  | CNN3 | 0.20 | 0.89 | 6.3e-29 |
  | HSPB1 | 0.20 | 0.88 | 1.4e-27 |
  | MAT2B | 0.15 | 0.78 | 4.9e-21 |
  | LGALS1 | 0.13 | 0.74 | 4.7e-19 |
  | ACAA2 | 0.18 | 0.78 | 2.3e-18 |
  | ALDH6A1 | 0.21 | 0.80 | 7.3e-18 |
  | ISYNA1 | 0.20 | 0.77 | 3.8e-17 |
  | HADH | 0.28 | 0.84 | 2.0e-17 |
  | TMSB4X | 0.31 | 0.86 | 2.0e-17 |
  | NEDD4L | 0.10 | 0.64 | 3.8e-15 |
  | LGALS3 | 0.11 | 0.65 | 1.8e-14 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | SELENBP1 | +3.38 | 2.4e-06 |
  | S100B | +2.55 | 1.7e-10 |
  | BLVRB | +2.48 | 3.8e-03 |
  | SLC2A1 | +2.41 | 2.0e-02 |
  | UGDH | +1.97 | 2.2e-23 |
  | SIRT2 | +1.46 | 6.8e-03 |
  | HBG2 | +1.42 | 3.9e-03 |
  | APPL1 | +1.21 | 8.4e-03 |
  | COTL1 | +1.19 | 4.1e-04 |
  | HBB | +1.15 | 8.1e-03 |
  | MAP2 | +1.09 | 8.6e-09 |
  | CTSD | +1.09 | 6.7e-09 |
  | ARL8B | +1.08 | 1.3e-02 |
  | MMP3 | +1.08 | 1.5e-02 |
  | SCARB2 | +1.02 | 1.4e-07 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | VIM | -3.93 | 3.1e-08 |
  | THY1 | -3.77 | 6.5e-03 |
  | VIM | -3.13 | 4.1e-04 |
  | CKB | -2.21 | 4.4e-23 |
  | SEPTIN5 | -1.68 | 2.8e-03 |
  | HDGFL3 | -1.63 | 2.4e-03 |
  | PDXP | -1.52 | 7.4e-03 |
  | TAGLN2 | -1.41 | 1.5e-02 |
  | SNRPG | -1.41 | 2.8e-02 |
  | KIF5C | -1.39 | 4.1e-02 |
  | ATAT1 | -1.38 | 1.9e-02 |
  | MAPT | -1.34 | 1.2e-04 |
  | EZR | -1.30 | 2.4e-19 |
  | PHF6 | -1.29 | 1.6e-02 |
  | GNA11 | -1.29 | 5.5e-03 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | S100B | +3.88 | 7.2e-07 |
  | BAIAP2 | +3.61 | 1.4e-07 |
  | SELENBP1 | +2.81 | 6.9e-06 |
  | G6PD | +2.28 | 4.0e-03 |
  | EEF1A2 | +2.27 | 1.7e-04 |
  | PITPNB | +2.23 | 1.6e-02 |
  | UGDH | +2.14 | 4.6e-25 |
  | SCARB2 | +2.04 | 2.4e-07 |
  | PGM3 | +1.85 | 2.9e-09 |
  | NME2P1 | +1.83 | 4.7e-02 |
  | CNRIP1 | +1.77 | 6.9e-06 |
  | QKI | +1.62 | 4.0e-11 |
  | CTSD | +1.58 | 1.4e-15 |
  | FHL1 | +1.57 | 1.5e-04 |
  | PSMC4 | +1.54 | 7.8e-03 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | VIM | -2.47 | 2.0e-07 |
  | VIM | -2.47 | 2.0e-07 |
  | CA2 | -2.32 | 9.4e-03 |
  | CKB | -2.03 | 1.3e-19 |
  | AHNAK | -1.84 | 3.9e-02 |
  | RAB21 | -1.79 | 4.4e-03 |
  | HMGB2 | -1.26 | 2.5e-11 |
  | HMGB2 | -1.26 | 2.5e-11 |
  | CDSN | -1.24 | 1.8e-02 |
  | SMARCC2 | -1.19 | 2.5e-02 |
  | SMARCC2 | -1.19 | 2.5e-02 |
  | FKBP1A | -1.19 | 2.6e-02 |
  | CROCC | -1.18 | 1.9e-04 |
  | H3-3B | -1.13 | 4.4e-02 |
  | LGALS1 | -1.10 | 4.3e-02 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_HEDGEHOG_SIGNALING | 0.114 | +0.026 |
  | HALLMARK_ALLOGRAFT_REJECTION | 0.048 | +0.009 |
  | HALLMARK_SPERMATOGENESIS | 0.015 | +0.008 |
  | HALLMARK_XENOBIOTIC_METABOLISM | 0.012 | +0.008 |
  | HALLMARK_HEME_METABOLISM | 0.022 | +0.007 |
  | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | 0.050 | +0.007 |
  | HALLMARK_ADIPOGENESIS | 0.013 | +0.006 |
  | HALLMARK_ESTROGEN_RESPONSE_EARLY | 0.029 | +0.004 |
  | HALLMARK_HYPOXIA | 0.049 | +0.004 |
  | HALLMARK_BILE_ACID_METABOLISM | 0.007 | +0.002 |
  | HALLMARK_DNA_REPAIR | 0.009 | +0.001 |
  | HALLMARK_PEROXISOME | 0.028 | +0.001 |
  | HALLMARK_ANDROGEN_RESPONSE | 0.018 | +0.001 |
  | HALLMARK_GLYCOLYSIS | 0.066 | +0.001 |
  | HALLMARK_UV_RESPONSE_UP | 0.015 | +0.001 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 3.07 | 15.38 | 1.09 | 0.1614 | 0.0000 |
  | TUBB3 | 1.00 | 0.96 | 2.36 | 15.38 | -0.30 | 0.3848 | 0.2983 |
  | DCX | 0.89 | 0.63 | 2.09 | 13.63 | 0.29 | 0.0004 | 0.4886 |
  | STMN1 | 1.00 | 1.00 | -3.62 | 15.72 | -0.28 | 1.0000 | 0.0264 |
  | STMN2 | 0.43 | 0.53 | -0.60 | 15.13 | -0.16 | 0.3650 | 0.2130 |
  | GAP43 | 0.05 | 0.08 | -0.61 | — | — | 0.7443 | — |
  | SYN1 | 0.02 | 0.02 | 0.25 | — | — | 1.0000 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.02 | 0.03 | -0.37 | — | — | 1.0000 | — |
  | NEFL | 0.02 | 0.03 | -0.49 | — | — | 0.9448 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.90 | 0.72 | 1.71 | 13.27 | 0.25 | 0.0137 | 0.7095 |
  | UCHL1 | 0.98 | 0.92 | 1.88 | 15.86 | 0.52 | 0.2434 | 0.0019 |
  | SOX2 | 0.08 | 0.11 | -0.42 | 12.46 | 0.10 | 0.8059 | 0.9769 |
  | NES | 0.43 | 0.76 | -2.11 | 14.24 | -0.36 | 0.0000 | 0.0506 |
  | VIM | 0.30 | 0.55 | -1.51 | 14.32 | -3.13 | 0.0023 | 0.0004 |
  | PAX6 | 0.02 | 0.05 | -1.18 | — | — | 0.6564 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 1.00 | 0.95 | 2.76 | 18.18 | -1.11 | 0.2984 | 0.0000 |
  | HOPX | 0.02 | 0.32 | -4.26 | — | — | 0.0000 | — |
  | EOMES | 0.00 | 0.09 | -3.61 | — | — | 0.0461 | — |
  | TBR1 | 0.03 | 0.11 | -1.56 | — | — | 0.2428 | — |
  | SATB2 | 0.07 | 0.01 | 2.63 | 16.70 | -0.02 | 0.0680 | 0.8379 |
  | BCL11B | 0.00 | 0.02 | -1.54 | — | — | 0.8704 | — |
  | CUX1 | 0.02 | 0.03 | -0.37 | — | — | 1.0000 | — |
  | CUX2 | 0.00 | 0.07 | -3.19 | — | — | 0.1189 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | -0.79 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.46 | 0.77 | -2.00 | 15.58 | 0.42 | 0.0000 | 0.1516 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.79 | 0.06 | 5.85 | 15.83 | 2.55 | 0.0000 | 0.0000 |
  | SLC1A3 | 0.00 | 0.03 | -1.77 | — | — | 0.6738 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.03 | 0.01 | 2.82 | — | — | 0.2221 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.46 | 0.11 | 2.76 | 12.45 | 0.44 | 0.0000 | 0.1390 |
  | AIF1 | 0.00 | 0.06 | -3.10 | — | — | 0.1632 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.07 | 0.11 | -0.70 | — | — | 0.6738 | — |
  | COL1A2 | 0.30 | 0.17 | 1.07 | 14.22 | -0.22 | 0.1057 | 0.7461 |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.04 | -2.21 | — | — | 0.5196 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.33 | 0.54 | -1.25 | 14.62 | 0.58 | 0.0176 | 0.0095 |
  | TOP2A | 0.56 | 0.60 | -0.26 | 13.82 | -0.52 | 0.8695 | 0.1193 |
  | HMGB2 | 0.59 | 0.84 | -1.89 | 12.98 | -1.01 | 0.0002 | 0.0000 |
  | CDK1 | 0.74 | 0.70 | 0.23 | 13.05 | -0.72 | 0.9023 | 0.0005 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C1

- **n_cells**: 75
- **median_n_detected**: 1202.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | CRMP1 | 0.67 | 3/3 |
  | DPYSL3 | 0.63 | 3/3 |
  | TUBB3 | 0.59 | 3/3 |
  | HDGF | 0.57 | 3/3 |
  | DCX | 0.55 | 3/3 |
  | YWHAG | 0.54 | 2/3 |
  | UCHL1 | 0.54 | 3/3 |
  | HP1BP3 | 0.50 | 2/3 |
  | PTMS | 0.50 | 2/3 |
  | VIM | 0.49 | 3/3 |
  | PAFAH1B3 | 0.49 | 3/3 |
  | PPP1R14B | 0.47 | 3/3 |
  | DDAH2 | 0.44 | 3/3 |
  | SUB1 | 0.44 | 2/3 |
  | H1-10 | 0.43 | 2/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | TOP2B | 1.00 | 0.82 | 4.4e-06 |
  | DPYSL3 | 1.00 | 0.86 | 8.9e-05 |
  | PGRMC1 | 1.00 | 0.94 | 4.7e-02 |
  | TRIM28 | 1.00 | 0.92 | 8.5e-03 |
  | SF3B1 | 1.00 | 0.78 | 1.7e-07 |
  | FKBP3 | 1.00 | 0.86 | 8.9e-05 |
  | CBX1 | 1.00 | 0.60 | 2.3e-15 |
  | SARNP | 1.00 | 0.90 | 1.5e-03 |
  | BASP1 | 1.00 | 0.72 | 4.1e-10 |
  | PGLS | 1.00 | 0.93 | 1.3e-02 |
  | DCX | 1.00 | 0.61 | 1.1e-14 |
  | TXNL1 | 1.00 | 0.85 | 3.3e-05 |
  | DHX15 | 1.00 | 0.84 | 1.2e-05 |
  | DDAH1 | 1.00 | 0.90 | 1.5e-03 |
  | UBE2I | 1.00 | 0.78 | 9.4e-08 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | MCM3 | 0.15 | 0.82 | 3.4e-31 |
  | MCM5 | 0.16 | 0.82 | 4.8e-30 |
  | MCM7 | 0.16 | 0.82 | 1.3e-29 |
  | MCM4 | 0.15 | 0.81 | 3.7e-29 |
  | MCM6 | 0.15 | 0.80 | 4.9e-29 |
  | MCM2 | 0.25 | 0.83 | 1.9e-23 |
  | MKI67 | 0.03 | 0.58 | 1.0e-21 |
  | DUT | 0.05 | 0.58 | 2.1e-19 |
  | FLNA | 0.29 | 0.82 | 2.1e-19 |
  | TOP2A | 0.13 | 0.64 | 4.0e-17 |
  | RRM1 | 0.11 | 0.62 | 1.6e-17 |
  | NES | 0.28 | 0.78 | 1.6e-17 |
  | PSAT1 | 0.35 | 0.85 | 3.4e-19 |
  | QKI | 0.24 | 0.74 | 1.4e-16 |
  | LINC01600 | 0.12 | 0.54 | 2.0e-12 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | FBL | +4.35 | 4.1e-19 |
  | CRMP1 | +3.40 | 9.7e-37 |
  | TUBB3 | +3.25 | 9.9e-37 |
  | DCX | +3.16 | 9.5e-33 |
  | DPYSL3 | +3.01 | 4.4e-35 |
  | SNW1 | +3.01 | 2.2e-07 |
  | SRM | +2.88 | 5.0e-25 |
  | LMNB2 | +2.76 | 4.3e-29 |
  | HNRNPAB | +2.75 | 3.0e-29 |
  | STMN1 | +2.68 | 1.2e-35 |
  | HDGF | +2.64 | 3.5e-37 |
  | PPP1R14B | +2.61 | 1.3e-32 |
  | DDAH2 | +2.60 | 7.3e-37 |
  | YWHAG | +2.58 | 7.3e-37 |
  | UCHL1 | +2.53 | 6.4e-35 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | VIM | -4.24 | 1.1e-21 |
  | SMNDC1 | -4.03 | 1.5e-08 |
  | IGF2BP1 | -4.03 | 4.0e-05 |
  | KCTD12 | -3.40 | 6.8e-06 |
  | ALB | -2.96 | 5.7e-12 |
  | PCNA | -2.90 | 1.7e-15 |
  | H2AC6 | -2.49 | 2.4e-05 |
  | PSAT1 | -2.42 | 3.2e-06 |
  | PRKAR1A | -1.88 | 3.5e-03 |
  | OSBP | -1.88 | 7.9e-03 |
  | HSPB1 | -1.72 | 3.2e-08 |
  | RAB35 | -1.55 | 7.7e-10 |
  | HSPH1 | -1.55 | 5.1e-03 |
  | XPOT | -1.51 | 1.2e-02 |
  | LGALS3 | -1.45 | 7.7e-05 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | MUCL1 | +3.92 | 2.8e-03 |
  | KRT6B | +3.02 | 3.7e-02 |
  | KRT6B | +3.02 | 3.7e-02 |
  | DPYSL3 | +2.94 | 3.0e-57 |
  | TUBB3 | +2.90 | 1.1e-47 |
  | CRMP1 | +2.87 | 1.5e-58 |
  | H1-0 | +2.79 | 2.4e-41 |
  | HPRT1 | +2.74 | 2.6e-02 |
  | UCHL1 | +2.63 | 3.2e-57 |
  | HDGF | +2.55 | 6.8e-59 |
  | DCX | +2.49 | 1.3e-38 |
  | YWHAG | +2.42 | 1.6e-60 |
  | SRM | +2.32 | 1.5e-34 |
  | HP1BP3 | +2.25 | 9.6e-64 |
  | NFIX | +2.24 | 2.5e-08 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | TERF2IP | -5.92 | 4.9e-20 |
  | AP2A2 | -4.95 | 1.4e-07 |
  | VIM | -3.69 | 4.2e-33 |
  | VIM | -3.69 | 4.2e-33 |
  | MYLK | -2.79 | 9.6e-03 |
  | HELZ2 | -2.61 | 4.9e-03 |
  | PCNA | -2.58 | 7.6e-30 |
  | RAB21 | -2.06 | 3.0e-07 |
  | DNHD1 | -1.90 | 1.7e-02 |
  | C11orf54 | -1.55 | 2.5e-05 |
  | RBPMS | -1.38 | 5.4e-08 |
  | OCLN | -1.33 | 2.0e-06 |
  | AHNAK | -1.32 | 3.5e-02 |
  | IAH1 | -1.28 | 2.6e-03 |
  | HSPB1 | -1.28 | 3.0e-09 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_HEDGEHOG_SIGNALING | 0.136 | +0.050 |
  | HALLMARK_PANCREAS_BETA_CELLS | 0.043 | +0.040 |
  | HALLMARK_G2M_CHECKPOINT | 0.108 | +0.022 |
  | HALLMARK_E2F_TARGETS | 0.094 | +0.017 |
  | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | 0.058 | +0.016 |
  | HALLMARK_SPERMATOGENESIS | 0.018 | +0.011 |
  | HALLMARK_MYC_TARGETS_V1 | 0.188 | +0.011 |
  | HALLMARK_ADIPOGENESIS | 0.014 | +0.007 |
  | HALLMARK_ANDROGEN_RESPONSE | 0.023 | +0.006 |
  | HALLMARK_HEME_METABOLISM | 0.019 | +0.005 |
  | HALLMARK_WNT_BETA_CATENIN_SIGNALING | 0.004 | +0.004 |
  | HALLMARK_XENOBIOTIC_METABOLISM | 0.008 | +0.004 |
  | HALLMARK_ANGIOGENESIS | 0.000 | -0.000 |
  | HALLMARK_NOTCH_SIGNALING | 0.000 | -0.000 |
  | HALLMARK_IL2_STAT5_SIGNALING | 0.001 | -0.001 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 3.39 | 16.49 | 2.26 | 0.0311 | 0.0000 |
  | TUBB3 | 1.00 | 0.96 | 2.69 | 18.81 | 3.25 | 0.1584 | 0.0000 |
  | DCX | 1.00 | 0.61 | 6.57 | 16.35 | 3.16 | 0.0000 | 0.0000 |
  | STMN1 | 1.00 | 1.00 | -3.30 | 18.52 | 2.68 | 1.0000 | 0.0000 |
  | STMN2 | 0.89 | 0.49 | 3.07 | 17.53 | 2.49 | 0.0000 | 0.0000 |
  | GAP43 | 0.47 | 0.04 | 4.31 | 12.60 | -0.15 | 0.0000 | 0.5279 |
  | SYN1 | 0.08 | 0.01 | 2.71 | 11.42 | -0.17 | 0.0039 | 0.4344 |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.17 | 0.01 | 3.78 | 12.35 | 0.47 | 0.0000 | 0.3074 |
  | NEFL | 0.24 | 0.01 | 4.80 | 12.38 | 0.40 | 0.0000 | 0.9018 |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 1.00 | 0.71 | 5.95 | 15.15 | 2.27 | 0.0000 | 0.0000 |
  | UCHL1 | 1.00 | 0.91 | 3.84 | 17.81 | 2.53 | 0.0056 | 0.0000 |
  | SOX2 | 0.01 | 0.12 | -2.80 | — | — | 0.0037 | — |
  | NES | 0.28 | 0.78 | -3.21 | 14.00 | -0.60 | 0.0000 | 0.1597 |
  | VIM | 0.72 | 0.51 | 1.29 | 13.32 | -4.24 | 0.0014 | 0.0000 |
  | PAX6 | 0.00 | 0.06 | -3.16 | — | — | 0.0466 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 1.00 | 0.95 | 3.08 | 19.06 | -0.18 | 0.0705 | 0.0509 |
  | HOPX | 0.05 | 0.32 | -2.93 | — | — | 0.0000 | — |
  | EOMES | 0.01 | 0.09 | -2.31 | — | — | 0.0276 | — |
  | TBR1 | 0.73 | 0.04 | 5.99 | 13.49 | -0.30 | 0.0000 | 0.5649 |
  | SATB2 | 0.07 | 0.01 | 2.76 | 13.32 | -0.45 | 0.0082 | 0.9673 |
  | BCL11B | 0.07 | 0.02 | 2.19 | 12.34 | 0.76 | 0.0268 | 0.2105 |
  | CUX1 | 0.00 | 0.03 | -2.36 | — | — | 0.2346 | — |
  | CUX2 | 0.48 | 0.02 | 5.35 | 12.93 | -0.23 | 0.0000 | 0.1910 |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.05 | 0.01 | 2.84 | — | — | 0.0170 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.87 | 0.74 | 1.15 | 15.21 | 0.03 | 0.0223 | 0.3155 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.00 | 0.12 | -4.42 | — | — | 0.0004 | — |
  | SLC1A3 | 0.00 | 0.03 | -2.10 | — | — | 0.3440 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | -0.41 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.05 | 0.15 | -1.44 | — | — | 0.0563 | — |
  | AIF1 | 0.00 | 0.07 | -3.43 | — | — | 0.0319 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.11 | 0.11 | 0.04 | 12.37 | -1.46 | 1.0000 | 0.0681 |
  | COL1A2 | 0.07 | 0.19 | -1.57 | 14.68 | 0.32 | 0.0128 | 0.8681 |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.04 | -2.53 | — | — | 0.2409 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.03 | 0.58 | -5.31 | — | — | 0.0000 | — |
  | TOP2A | 0.13 | 0.64 | -3.49 | 14.23 | -0.09 | 0.0000 | 0.7823 |
  | HMGB2 | 1.00 | 0.81 | 5.20 | 15.32 | 1.50 | 0.0000 | 0.0000 |
  | CDK1 | 0.79 | 0.70 | 0.66 | 13.25 | -0.55 | 0.1729 | 0.0000 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C10

- **n_cells**: 31
- **median_n_detected**: 853.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | RRM1 | 0.46 | 3/3 |
  | EZR | 0.42 | 2/3 |
  | PHGDH | 0.39 | 3/3 |
  | KRT84 | 0.37 | 2/3 |
  | AKAP12 | 0.34 | 3/3 |
  | MAT2B | 0.29 | 3/3 |
  | SCRN1 | 0.28 | 3/3 |
  | VIM | 0.24 | 3/3 |
  | RBBP7 | 0.24 | 3/3 |
  | PAICS | 0.23 | 3/3 |
  | SYNE2 | 0.21 | 3/3 |
  | SMC4 | 0.21 | 3/3 |
  | PSAT1 | 0.21 | 3/3 |
  | CKB | 0.21 | 2/3 |
  | TERF2IP | 0.20 | 1/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | MSI1 | 1.00 | 0.66 | 9.5e-04 |
  | TMSB15A | 1.00 | 0.78 | 3.1e-02 |
  | MCM3 | 1.00 | 0.75 | 1.4e-02 |
  | GART | 1.00 | 0.72 | 4.7e-03 |
  | RRM1 | 1.00 | 0.55 | 1.6e-05 |
  | SRM | 1.00 | 0.73 | 8.2e-03 |
  | FLNA | 1.00 | 0.76 | 1.4e-02 |
  | MAT2B | 1.00 | 0.72 | 6.6e-03 |
  | SAE1 | 1.00 | 0.80 | 4.4e-02 |
  | RUVBL2 | 1.00 | 0.80 | 4.5e-02 |
  | SUPT16H | 1.00 | 0.79 | 3.1e-02 |
  | HADH | 1.00 | 0.79 | 3.1e-02 |
  | MCM7 | 1.00 | 0.75 | 1.4e-02 |
  | MCM5 | 1.00 | 0.75 | 1.4e-02 |
  | MCM4 | 1.00 | 0.73 | 8.2e-03 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | CALM1 | 0.13 | 0.49 | 5.2e-03 |
  | RBMX | 0.16 | 0.50 | 8.2e-03 |
  | EEF1A1 | 0.16 | 0.50 | 8.2e-03 |
  | UBE2N | 0.16 | 0.50 | 8.2e-03 |
  | UBE2V1 | 0.16 | 0.50 | 8.2e-03 |
  | KHDRBS1 | 0.16 | 0.50 | 8.2e-03 |
  | RAB5A | 0.16 | 0.50 | 8.2e-03 |
  | KRT13 | 0.16 | 0.49 | 1.4e-02 |
  | MAGOH | 0.16 | 0.48 | 1.5e-02 |
  | SRSF4 | 0.16 | 0.47 | 2.3e-02 |
  | RAB35 | 0.10 | 0.40 | 1.8e-02 |
  | DYNLL1 | 0.32 | 0.62 | 3.4e-02 |
  | H3-7 | 0.16 | 0.46 | 3.5e-02 |
  | TUBB2A | 0.16 | 0.45 | 3.5e-02 |
  | HSPA6 | 0.16 | 0.45 | 3.7e-02 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | VIM | +1.67 | 6.7e-04 |
  | HNRNPAB | +1.27 | 1.7e-03 |
  | RRM2 | +1.02 | 7.0e-03 |
  | PHGDH | +0.96 | 8.3e-10 |
  | EZR | +0.95 | 5.9e-10 |
  | TOP2B | +0.93 | 1.5e-02 |
  | CSTB | +0.92 | 2.0e-02 |
  | TPM4 | +0.92 | 3.7e-02 |
  | LMNB2 | +0.87 | 3.4e-02 |
  | RAB5A | +0.84 | 1.1e-02 |
  | SCRN1 | +0.80 | 2.4e-08 |
  | FBL | +0.78 | 2.8e-02 |
  | CROCC | +0.78 | 2.8e-03 |
  | STMN1 | +0.77 | 2.6e-03 |
  | RAB6A | +0.77 | 5.9e-05 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | ALB | -2.43 | 2.9e-02 |
  | SSBP2 | -1.58 | 1.6e-02 |
  | SMARCD3 | -1.55 | 1.9e-03 |
  | PSAP | -1.51 | 3.6e-03 |
  | SUMO1 | -1.44 | 4.7e-02 |
  | RAP1B | -1.22 | 6.7e-03 |
  | KRT5 | -1.18 | 1.4e-03 |
  | EFHD2 | -1.17 | 2.4e-02 |
  | ARHGAP1 | -1.15 | 8.5e-03 |
  | YLPM1 | -1.14 | 1.4e-02 |
  | LSM14B | -1.12 | 2.1e-02 |
  | LAMP1 | -1.11 | 3.4e-02 |
  | AK1 | -1.07 | 1.7e-02 |
  | ELAVL3 | -1.04 | 4.4e-03 |
  | ANXA2 | -0.99 | 3.1e-03 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | KRT86 | +4.04 | 7.4e-04 |
  | PSAT1 | +1.28 | 1.3e-11 |
  | FHL1 | +1.12 | 3.3e-02 |
  | BLVRA | +1.09 | 1.3e-03 |
  | TUBB2A | +1.09 | 2.1e-02 |
  | TUBB2A | +1.09 | 2.1e-02 |
  | PEA15 | +1.06 | 1.8e-10 |
  | PPA2 | +1.06 | 1.0e-03 |
  | SUCLG1 | +1.04 | 4.5e-02 |
  | NIT1 | +0.99 | 7.3e-03 |
  | ALDH6A1 | +0.98 | 5.3e-08 |
  | UBE2V1 | +0.97 | 3.0e-03 |
  | TCERG1 | +0.94 | 3.3e-02 |
  | RAB2A | +0.94 | 4.6e-02 |
  | GSTM2 | +0.94 | 4.2e-02 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | TERF2IP | -7.21 | 9.0e-04 |
  | KRT84 | -4.91 | 4.8e-08 |
  | VCL | -2.09 | 7.9e-03 |
  | WASHC2A | -2.00 | 4.5e-02 |
  | IAH1 | -1.96 | 4.6e-02 |
  | A2ML1 | -1.92 | 1.7e-02 |
  | FKBP1A | -1.78 | 4.7e-02 |
  | SMARCC2 | -1.27 | 4.7e-02 |
  | SMARCC2 | -1.27 | 4.7e-02 |
  | CD207 | -1.27 | 2.2e-02 |
  | POLR2G | -0.76 | 4.6e-02 |
  | ELAVL3 | -0.66 | 3.0e-03 |
  | ATL1 | -0.66 | 1.4e-03 |
  | ACTA1 | -0.63 | 3.0e-02 |
  | COTL1 | -0.51 | 4.3e-02 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_APICAL_SURFACE | 0.028 | +0.014 |
  | HALLMARK_MYC_TARGETS_V1 | 0.190 | +0.012 |
  | HALLMARK_E2F_TARGETS | 0.087 | +0.009 |
  | HALLMARK_G2M_CHECKPOINT | 0.095 | +0.006 |
  | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | 0.048 | +0.005 |
  | HALLMARK_DNA_REPAIR | 0.012 | +0.005 |
  | HALLMARK_MYC_TARGETS_V2 | 0.026 | +0.004 |
  | HALLMARK_GLYCOLYSIS | 0.069 | +0.004 |
  | HALLMARK_P53_PATHWAY | 0.027 | +0.003 |
  | HALLMARK_ESTROGEN_RESPONSE_LATE | 0.036 | +0.002 |
  | HALLMARK_MTORC1_SIGNALING | 0.073 | +0.002 |
  | HALLMARK_APOPTOSIS | 0.016 | +0.002 |
  | HALLMARK_ANDROGEN_RESPONSE | 0.019 | +0.001 |
  | HALLMARK_PI3K_AKT_MTOR_SIGNALING | 0.042 | +0.001 |
  | HALLMARK_MITOTIC_SPINDLE | 0.036 | +0.001 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 2.04 | 14.53 | 0.22 | 0.7812 | 0.2553 |
  | TUBB3 | 1.00 | 0.96 | 1.34 | 15.58 | -0.10 | 1.0000 | 0.8076 |
  | DCX | 0.77 | 0.64 | 0.85 | 13.19 | -0.20 | 0.6551 | 0.3564 |
  | STMN1 | 1.00 | 1.00 | -4.65 | 16.70 | 0.77 | 1.0000 | 0.0026 |
  | STMN2 | 0.84 | 0.51 | 2.21 | 15.46 | 0.20 | 0.0136 | 0.5179 |
  | GAP43 | 0.00 | 0.08 | -2.54 | — | — | 0.6317 | — |
  | SYN1 | 0.00 | 0.02 | -0.43 | — | — | 1.0000 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.00 | 0.03 | -1.01 | — | — | 1.0000 | — |
  | NEFL | 0.00 | 0.03 | -1.13 | — | — | 1.0000 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.97 | 0.73 | 2.93 | 12.93 | -0.12 | 0.0352 | 0.4180 |
  | UCHL1 | 1.00 | 0.92 | 2.49 | 15.58 | 0.22 | 0.6317 | 0.6347 |
  | SOX2 | 0.03 | 0.12 | -1.42 | — | — | 0.7812 | — |
  | NES | 0.97 | 0.73 | 2.92 | 14.38 | -0.22 | 0.0352 | 0.2209 |
  | VIM | 0.84 | 0.52 | 2.17 | 17.62 | 0.26 | 0.0139 | 0.2707 |
  | PAX6 | 0.03 | 0.05 | -0.14 | — | — | 1.0000 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 1.00 | 0.95 | 1.73 | 19.65 | 0.46 | 0.9504 | 0.0014 |
  | HOPX | 0.19 | 0.30 | -0.77 | 12.69 | -0.77 | 0.7812 | 0.2840 |
  | EOMES | 0.32 | 0.07 | 2.61 | 12.67 | 0.49 | 0.0065 | 0.2827 |
  | TBR1 | 0.00 | 0.11 | -2.94 | — | — | 0.4221 | — |
  | SATB2 | 0.00 | 0.02 | -0.13 | 17.09 | 0.41 | 1.0000 | 0.1278 |
  | BCL11B | 0.00 | 0.02 | -0.52 | — | — | 1.0000 | — |
  | CUX1 | 0.00 | 0.03 | -1.01 | — | — | 1.0000 | — |
  | CUX2 | 0.00 | 0.07 | -2.17 | — | — | 0.7812 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | 0.23 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.97 | 0.74 | 2.83 | 15.28 | 0.10 | 0.0536 | 0.2362 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.03 | 0.12 | -1.42 | — | — | 0.7812 | — |
  | SLC1A3 | 0.00 | 0.03 | -0.75 | — | — | 1.0000 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | 0.93 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.13 | 0.14 | 0.04 | — | — | 1.0000 | — |
  | AIF1 | 0.00 | 0.06 | -2.08 | — | — | 0.7812 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.10 | 0.11 | -0.00 | — | — | 1.0000 | — |
  | COL1A2 | 0.10 | 0.18 | -0.84 | — | — | 0.9434 | — |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.03 | -1.19 | — | — | 1.0000 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.97 | 0.51 | 4.30 | 14.18 | 0.13 | 0.0000 | 0.4627 |
  | TOP2A | 0.97 | 0.58 | 3.87 | 14.80 | 0.53 | 0.0005 | 0.0599 |
  | HMGB2 | 1.00 | 0.82 | 3.83 | 14.34 | 0.42 | 0.0656 | 0.0505 |
  | CDK1 | 0.87 | 0.70 | 1.40 | 14.38 | 0.72 | 0.3296 | 0.0000 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C11

- **n_cells**: 48
- **median_n_detected**: 601.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | TUFM | 0.36 | 3/3 |
  | DDX58 | 0.35 | 2/3 |
  | FBL | 0.33 | 2/3 |
  | COTL1 | 0.28 | 2/3 |
  | TOP2A | 0.26 | 3/3 |
  | RAB11B | 0.25 | 3/3 |
  | PCNA | 0.25 | 3/3 |
  | ARPC4 | 0.22 | 3/3 |
  | AK1 | 0.22 | 3/3 |
  | COX6C | 0.22 | 3/3 |
  | NPEPPS | 0.21 | 3/3 |
  | ATP5MG | 0.21 | 3/3 |
  | BOLA2B | 0.21 | 3/3 |
  | PAK2 | 0.20 | 2/3 |
  | MAP2 | 0.20 | 3/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | TAGLN2 | 1.00 | 0.83 | 2.6e-03 |
  | PCNA | 1.00 | 0.84 | 5.8e-03 |
  | PSAT1 | 0.98 | 0.79 | 5.4e-03 |
  | ACTBL2 | 0.94 | 0.77 | 4.2e-02 |
  | RAB11B | 0.92 | 0.99 | 6.3e-03 |
  | MCM4 | 0.92 | 0.73 | 2.6e-02 |
  | CCT4 | 0.90 | 0.98 | 2.1e-02 |
  | QKI | 0.90 | 0.68 | 1.0e-02 |
  | EEF1B2 | 0.90 | 0.98 | 4.0e-02 |
  | MSI1 | 0.90 | 0.66 | 4.4e-03 |
  | PGD | 0.90 | 1.00 | 5.5e-04 |
  | DSC1 | 0.85 | 0.53 | 1.2e-04 |
  | MDH1 | 0.85 | 0.99 | 4.0e-04 |
  | PPA2 | 0.85 | 0.47 | 9.2e-06 |
  | LMNB2 | 0.85 | 0.98 | 3.9e-03 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | COX6C | 0.19 | 0.71 | 9.6e-10 |
  | ATP5MG | 0.29 | 0.79 | 9.6e-10 |
  | COTL1 | 0.08 | 0.58 | 1.5e-09 |
  | PSMA7 | 0.21 | 0.70 | 6.3e-09 |
  | ARPC5 | 0.21 | 0.69 | 1.2e-08 |
  | ARPC4 | 0.33 | 0.82 | 1.5e-09 |
  | RAB18 | 0.12 | 0.60 | 2.5e-08 |
  | NAGK | 0.21 | 0.68 | 3.1e-08 |
  | ATIC | 0.33 | 0.79 | 1.6e-08 |
  | TUFM | 0.48 | 0.94 | 1.9e-12 |
  | SRM | 0.31 | 0.77 | 2.9e-08 |
  | COX5B | 0.23 | 0.68 | 8.9e-08 |
  | NPEPPS | 0.42 | 0.87 | 1.5e-09 |
  | DLST | 0.29 | 0.74 | 9.8e-08 |
  | OLA1 | 0.27 | 0.71 | 2.0e-07 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | ALB | +2.95 | 1.7e-04 |
  | AK1 | +2.56 | 3.7e-02 |
  | NCALD | +1.76 | 2.4e-02 |
  | BOLA2B | +1.73 | 5.4e-04 |
  | DYNC1LI1 | +1.62 | 1.1e-03 |
  | LBR | +1.60 | 1.8e-02 |
  | SNRPC | +1.57 | 4.8e-02 |
  | CHAMP1 | +1.47 | 7.1e-03 |
  | LETM1 | +1.43 | 2.0e-03 |
  | PLEC | +1.34 | 1.1e-02 |
  | CGGBP1 | +1.26 | 3.7e-02 |
  | SEC16A | +1.22 | 1.1e-02 |
  | RPL22 | +1.18 | 4.6e-02 |
  | PPIH | +1.18 | 3.6e-02 |
  | CZIB | +1.13 | 4.6e-02 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | FBL | -4.35 | 1.2e-06 |
  | COX7A2 | -1.71 | 2.4e-07 |
  | KPNA2 | -1.66 | 5.9e-03 |
  | HNRNPAB | -1.52 | 9.8e-06 |
  | RAB15 | -1.47 | 2.6e-03 |
  | TOP2A | -1.47 | 5.5e-10 |
  | JPT1 | -1.47 | 2.8e-06 |
  | CSTB | -1.46 | 2.2e-07 |
  | RPS14 | -1.30 | 3.5e-02 |
  | LAP3 | -1.29 | 2.9e-02 |
  | RAB11B | -1.27 | 5.6e-12 |
  | TUBA4B | -1.22 | 2.7e-03 |
  | MAP2 | -1.18 | 1.6e-11 |
  | ARPC4 | -1.18 | 1.4e-02 |
  | LAMP1 | -1.17 | 8.8e-03 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | DDX58 | +5.65 | 5.2e-08 |
  | NFIX | +3.37 | 2.0e-03 |
  | LSM4 | +3.30 | 4.5e-02 |
  | MAPT | +2.76 | 3.2e-02 |
  | CDIPT | +2.55 | 6.4e-03 |
  | EML1 | +2.42 | 2.7e-03 |
  | UQCRFS1 | +2.35 | 4.2e-05 |
  | SEPTIN5 | +2.33 | 5.1e-03 |
  | CMPK1 | +2.32 | 8.1e-05 |
  | BCL7A | +2.29 | 2.3e-02 |
  | AASDHPPT | +2.28 | 6.2e-03 |
  | ARPC3 | +2.26 | 3.3e-06 |
  | AK1 | +2.25 | 1.6e-03 |
  | SHTN1 | +2.22 | 2.2e-02 |
  | ENTPD8 | +2.19 | 5.8e-04 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | TERF2IP | -5.59 | 1.0e-02 |
  | AP2A2 | -3.81 | 2.2e-02 |
  | RBPMS | -2.28 | 1.9e-07 |
  | TF | -2.08 | 1.3e-05 |
  | RAB21 | -1.79 | 4.7e-02 |
  | SERPINH1 | -1.56 | 4.4e-02 |
  | SMARCC2 | -1.35 | 1.7e-02 |
  | SMARCC2 | -1.35 | 1.7e-02 |
  | MKI67 | -1.30 | 9.0e-05 |
  | PCNA | -1.29 | 6.4e-15 |
  | TOP2A | -1.24 | 3.1e-09 |
  | TOP2A | -1.24 | 3.1e-09 |
  | KPNA2 | -1.16 | 1.9e-03 |
  | STMN1 | -1.08 | 5.8e-08 |
  | JPT1 | -1.07 | 5.5e-04 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_ESTROGEN_RESPONSE_EARLY | 0.034 | +0.009 |
  | HALLMARK_MYC_TARGETS_V1 | 0.183 | +0.005 |
  | HALLMARK_UNFOLDED_PROTEIN_RESPONSE | 0.055 | +0.005 |
  | HALLMARK_UV_RESPONSE_DN | 0.014 | +0.005 |
  | HALLMARK_P53_PATHWAY | 0.026 | +0.003 |
  | HALLMARK_APICAL_SURFACE | 0.018 | +0.003 |
  | HALLMARK_ESTROGEN_RESPONSE_LATE | 0.036 | +0.002 |
  | HALLMARK_UV_RESPONSE_UP | 0.017 | +0.002 |
  | HALLMARK_APOPTOSIS | 0.016 | +0.002 |
  | HALLMARK_E2F_TARGETS | 0.080 | +0.002 |
  | HALLMARK_CHOLESTEROL_HOMEOSTASIS | 0.033 | +0.001 |
  | HALLMARK_IL2_STAT5_SIGNALING | 0.002 | +0.001 |
  | HALLMARK_G2M_CHECKPOINT | 0.089 | +0.001 |
  | HALLMARK_WNT_BETA_CATENIN_SIGNALING | 0.001 | +0.001 |
  | HALLMARK_PEROXISOME | 0.027 | +0.000 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 2.70 | 13.24 | -1.18 | 0.3592 | 0.0000 |
  | TUBB3 | 1.00 | 0.96 | 2.00 | 15.07 | -0.65 | 0.6205 | 0.0000 |
  | DCX | 0.38 | 0.67 | -1.72 | 12.61 | -0.80 | 0.0012 | 0.0252 |
  | STMN1 | 1.00 | 1.00 | -3.99 | 15.07 | -0.96 | 1.0000 | 0.0000 |
  | STMN2 | 0.29 | 0.54 | -1.46 | 14.58 | -0.73 | 0.0093 | 0.1797 |
  | GAP43 | 0.04 | 0.08 | -0.76 | — | — | 0.7893 | — |
  | SYN1 | 0.04 | 0.02 | 1.49 | — | — | 0.6111 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.02 | 0.03 | 0.01 | — | — | 1.0000 | — |
  | NEFL | 0.00 | 0.03 | -1.79 | — | — | 0.7635 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.42 | 0.76 | -2.11 | 12.06 | -0.99 | 0.0001 | 0.0027 |
  | UCHL1 | 0.98 | 0.92 | 1.50 | 14.81 | -0.67 | 0.4798 | 0.0000 |
  | SOX2 | 0.08 | 0.11 | -0.36 | — | — | 0.9611 | — |
  | NES | 0.79 | 0.74 | 0.40 | 14.51 | -0.09 | 0.8659 | 0.6659 |
  | VIM | 0.21 | 0.55 | -2.16 | 15.73 | -1.69 | 0.0001 | 0.2074 |
  | PAX6 | 0.00 | 0.05 | -2.46 | — | — | 0.4766 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 0.98 | 0.95 | 0.73 | 19.20 | -0.01 | 1.0000 | 0.8869 |
  | HOPX | 0.38 | 0.29 | 0.54 | 12.86 | -0.63 | 0.6205 | 0.0128 |
  | EOMES | 0.04 | 0.09 | -0.81 | — | — | 0.7893 | — |
  | TBR1 | 0.08 | 0.11 | -0.22 | — | — | 1.0000 | — |
  | SATB2 | 0.00 | 0.02 | -0.79 | 16.73 | 0.03 | 1.0000 | 0.7491 |
  | BCL11B | 0.04 | 0.02 | 1.39 | — | — | 0.6205 | — |
  | CUX1 | 0.00 | 0.03 | -1.67 | — | — | 0.7618 | — |
  | CUX2 | 0.06 | 0.06 | 0.17 | — | — | 1.0000 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | -0.42 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.77 | 0.75 | 0.14 | 14.26 | -0.96 | 1.0000 | 0.0002 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.00 | 0.12 | -3.72 | — | — | 0.0303 | — |
  | SLC1A3 | 0.00 | 0.03 | -1.40 | — | — | 0.9542 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.02 | 0.01 | 2.14 | — | — | 0.6757 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.06 | 0.14 | -1.10 | — | — | 0.5194 | — |
  | AIF1 | 0.00 | 0.06 | -2.73 | — | — | 0.3629 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.04 | 0.11 | -1.25 | — | — | 0.4612 | — |
  | COL1A2 | 0.17 | 0.18 | -0.05 | 14.06 | -0.36 | 1.0000 | 0.6715 |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.02 | 0.03 | -0.17 | — | — | 1.0000 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.33 | 0.54 | -1.19 | 12.94 | -1.13 | 0.0460 | 0.0003 |
  | TOP2A | 0.69 | 0.59 | 0.58 | 12.93 | -1.47 | 0.5848 | 0.0000 |
  | HMGB2 | 0.94 | 0.82 | 1.56 | 13.38 | -0.61 | 0.1530 | 0.0001 |
  | CDK1 | 0.71 | 0.70 | 0.00 | 13.08 | -0.69 | 1.0000 | 0.0001 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C12

- **n_cells**: 44
- **median_n_detected**: 662.5
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | STMN1 | 0.67 | 2/3 |
  | PEA15 | 0.57 | 2/3 |
  | STMN2 | 0.53 | 2/3 |
  | PCNA | 0.52 | 2/3 |
  | RCC2 | 0.39 | 3/3 |
  | HSPA5 | 0.39 | 2/3 |
  | TMSB15A | 0.38 | 3/3 |
  | NASP | 0.37 | 3/3 |
  | JPT1 | 0.35 | 3/3 |
  | TNC | 0.34 | 2/3 |
  | RANBP1 | 0.32 | 3/3 |
  | CDK1 | 0.30 | 3/3 |
  | PSAT1 | 0.30 | 3/3 |
  | LMNB1 | 0.30 | 2/3 |
  | H1-4 | 0.29 | 3/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | HSPB1 | 1.00 | 0.82 | 5.2e-03 |
  | FLNA | 1.00 | 0.76 | 3.2e-04 |
  | DARS1 | 1.00 | 0.84 | 1.1e-02 |
  | PSAT1 | 1.00 | 0.79 | 2.2e-03 |
  | TAGLN2 | 1.00 | 0.83 | 7.4e-03 |
  | SEPTIN9 | 1.00 | 0.84 | 1.7e-02 |
  | RCN1 | 1.00 | 0.78 | 8.9e-04 |
  | CALU | 1.00 | 0.84 | 1.1e-02 |
  | QKI | 1.00 | 0.68 | 7.8e-06 |
  | HMGCS1 | 1.00 | 0.73 | 8.1e-05 |
  | ALDH6A1 | 1.00 | 0.74 | 2.1e-04 |
  | NDRG2 | 0.98 | 0.82 | 4.4e-02 |
  | LGALS3 | 0.95 | 0.59 | 8.3e-06 |
  | TMED10 | 0.95 | 0.73 | 7.4e-03 |
  | GSTP1 | 0.93 | 0.99 | 3.5e-02 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | IMPDH2 | 0.05 | 0.66 | 5.6e-14 |
  | FEN1 | 0.11 | 0.70 | 3.5e-12 |
  | COTL1 | 0.02 | 0.58 | 4.5e-12 |
  | ATIC | 0.30 | 0.79 | 4.0e-09 |
  | DUT | 0.07 | 0.56 | 1.0e-08 |
  | DPYSL4 | 0.18 | 0.67 | 4.4e-08 |
  | SRM | 0.30 | 0.77 | 4.4e-08 |
  | AKR1B1 | 0.30 | 0.76 | 6.5e-08 |
  | CLIC4 | 0.02 | 0.49 | 6.2e-09 |
  | HSPA4 | 0.23 | 0.69 | 2.4e-07 |
  | SARS1 | 0.14 | 0.59 | 3.3e-07 |
  | NUDT5 | 0.16 | 0.61 | 7.6e-07 |
  | AKR1A1 | 0.05 | 0.49 | 1.2e-07 |
  | SSRP1 | 0.30 | 0.74 | 7.6e-07 |
  | PPP1CA | 0.30 | 0.73 | 1.0e-06 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | NELFCD | +1.99 | 6.0e-03 |
  | PEA15 | +1.52 | 3.3e-11 |
  | MFN1 | +1.18 | 2.5e-02 |
  | LPCAT1 | +1.00 | 1.9e-04 |
  | SATB2 | +0.94 | 4.0e-02 |
  | HSPA5 | +0.82 | 3.2e-10 |
  | PSAT1 | +0.76 | 2.1e-05 |
  | HMGCS1 | +0.75 | 2.6e-04 |
  | GLUD1 | +0.74 | 3.9e-06 |
  | CALR | +0.71 | 2.5e-06 |
  | GRM3 | +0.67 | 1.1e-02 |
  | FABP7 | +0.65 | 5.4e-09 |
  | HADHB | +0.62 | 2.3e-02 |
  | QKI | +0.60 | 1.4e-06 |
  | CALU | +0.59 | 5.8e-05 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | HDGFL3 | -2.26 | 1.3e-03 |
  | ELAVL3 | -1.95 | 1.1e-04 |
  | H3-4 | -1.84 | 1.4e-04 |
  | STMN2 | -1.78 | 2.5e-10 |
  | PCNA | -1.73 | 1.5e-09 |
  | VIM | -1.64 | 1.4e-02 |
  | ZNF512 | -1.57 | 2.9e-04 |
  | CDK1 | -1.52 | 4.5e-07 |
  | PSMA6 | -1.43 | 2.1e-04 |
  | JPT1 | -1.39 | 2.9e-09 |
  | OPA1 | -1.38 | 1.6e-06 |
  | ARPC2 | -1.35 | 9.7e-04 |
  | STMN1 | -1.33 | 3.1e-14 |
  | CAPZA2 | -1.33 | 7.7e-06 |
  | UBE2L3 | -1.31 | 2.6e-03 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | PLAAT3 | +4.49 | 2.0e-02 |
  | NELFCD | +2.44 | 2.1e-02 |
  | PEA15 | +1.91 | 6.3e-28 |
  | OGT | +1.84 | 4.6e-03 |
  | PSAT1 | +1.56 | 5.3e-23 |
  | QKI | +1.53 | 2.1e-10 |
  | HMGCS1 | +1.44 | 2.0e-13 |
  | HSPA5 | +1.40 | 9.9e-34 |
  | PTPRZ1 | +1.34 | 2.1e-05 |
  | NDRG2 | +1.23 | 4.3e-17 |
  | GLUD1 | +1.23 | 5.7e-12 |
  | TMED9 | +1.21 | 4.0e-14 |
  | PSAP | +1.20 | 3.7e-03 |
  | TMED7 | +1.17 | 2.5e-03 |
  | CUL3 | +1.17 | 3.1e-02 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | VCL | -2.31 | 1.7e-02 |
  | STMN2 | -1.99 | 1.3e-20 |
  | STMN1 | -1.92 | 2.6e-33 |
  | BAIAP2 | -1.87 | 1.1e-07 |
  | P4HA1 | -1.82 | 2.5e-02 |
  | NIBAN2 | -1.81 | 4.7e-03 |
  | PCNA | -1.81 | 2.6e-25 |
  | FKBP9 | -1.79 | 4.0e-05 |
  | JPT1 | -1.69 | 2.4e-15 |
  | SUPT5H | -1.59 | 1.8e-05 |
  | CDK1 | -1.51 | 1.4e-12 |
  | CDK1 | -1.51 | 1.4e-12 |
  | POU3F2 | -1.51 | 2.3e-04 |
  | C11orf54 | -1.50 | 3.8e-03 |
  | CDSN | -1.43 | 2.8e-02 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_UNFOLDED_PROTEIN_RESPONSE | 0.079 | +0.031 |
  | HALLMARK_APICAL_SURFACE | 0.025 | +0.010 |
  | HALLMARK_MYC_TARGETS_V2 | 0.030 | +0.008 |
  | HALLMARK_APOPTOSIS | 0.022 | +0.007 |
  | HALLMARK_MTORC1_SIGNALING | 0.078 | +0.007 |
  | HALLMARK_INTERFERON_ALPHA_RESPONSE | 0.014 | +0.005 |
  | HALLMARK_ANDROGEN_RESPONSE | 0.022 | +0.004 |
  | HALLMARK_ESTROGEN_RESPONSE_LATE | 0.037 | +0.003 |
  | HALLMARK_INTERFERON_GAMMA_RESPONSE | 0.009 | +0.003 |
  | HALLMARK_HEDGEHOG_SIGNALING | 0.093 | +0.003 |
  | HALLMARK_CHOLESTEROL_HOMEOSTASIS | 0.034 | +0.002 |
  | HALLMARK_BILE_ACID_METABOLISM | 0.008 | +0.002 |
  | HALLMARK_MYC_TARGETS_V1 | 0.180 | +0.002 |
  | HALLMARK_PI3K_AKT_MTOR_SIGNALING | 0.043 | +0.002 |
  | HALLMARK_NOTCH_SIGNALING | 0.001 | +0.001 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 2.57 | 13.84 | -0.55 | 0.3139 | 0.0001 |
  | TUBB3 | 1.00 | 0.96 | 1.86 | 15.32 | -0.39 | 0.7024 | 0.0022 |
  | DCX | 0.52 | 0.66 | -0.81 | 12.49 | -0.95 | 0.2518 | 0.0000 |
  | STMN1 | 1.00 | 1.00 | -4.12 | 14.72 | -1.33 | 1.0000 | 0.0000 |
  | STMN2 | 0.52 | 0.52 | -0.00 | 13.56 | -1.78 | 1.0000 | 0.0000 |
  | GAP43 | 0.00 | 0.09 | -3.06 | — | — | 0.1771 | — |
  | SYN1 | 0.00 | 0.02 | -0.95 | — | — | 1.0000 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.00 | 0.03 | -1.54 | — | — | 0.9140 | — |
  | NEFL | 0.00 | 0.03 | -1.65 | — | — | 0.7024 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.52 | 0.75 | -1.44 | 12.06 | -1.02 | 0.0232 | 0.0000 |
  | UCHL1 | 1.00 | 0.92 | 3.01 | 15.64 | 0.28 | 0.1771 | 0.3289 |
  | SOX2 | 0.05 | 0.12 | -1.17 | — | — | 0.5067 | — |
  | NES | 0.77 | 0.74 | 0.23 | 14.44 | -0.17 | 1.0000 | 0.1191 |
  | VIM | 0.75 | 0.52 | 1.45 | 17.46 | 0.08 | 0.0272 | 0.6300 |
  | PAX6 | 0.00 | 0.05 | -2.33 | — | — | 0.4161 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 1.00 | 0.95 | 2.25 | 19.83 | 0.65 | 0.5388 | 0.0000 |
  | HOPX | 0.77 | 0.27 | 3.13 | 13.43 | -0.01 | 0.0000 | 0.9899 |
  | EOMES | 0.00 | 0.09 | -3.11 | — | — | 0.1775 | — |
  | TBR1 | 0.02 | 0.11 | -1.83 | — | — | 0.2518 | — |
  | SATB2 | 0.00 | 0.02 | -0.66 | 17.62 | 0.94 | 1.0000 | 0.0400 |
  | BCL11B | 0.00 | 0.02 | -1.04 | — | — | 1.0000 | — |
  | CUX1 | 0.00 | 0.03 | -1.54 | — | — | 0.9140 | — |
  | CUX2 | 0.00 | 0.07 | -2.69 | — | — | 0.3139 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | -0.29 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.86 | 0.74 | 1.04 | 15.49 | 0.34 | 0.2520 | 0.1712 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.00 | 0.12 | -3.59 | — | — | 0.0674 | — |
  | SLC1A3 | 0.05 | 0.02 | 1.27 | — | — | 0.5890 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | 0.41 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.09 | 0.14 | -0.55 | — | — | 0.8262 | — |
  | AIF1 | 0.00 | 0.06 | -2.60 | — | — | 0.3139 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.11 | 0.11 | 0.19 | 12.90 | -0.86 | 1.0000 | 0.5476 |
  | COL1A2 | 0.07 | 0.18 | -1.42 | — | — | 0.2319 | — |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.03 | -1.71 | — | — | 0.7024 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.14 | 0.55 | -2.84 | 13.70 | -0.38 | 0.0000 | 0.2534 |
  | TOP2A | 0.30 | 0.61 | -1.89 | 13.06 | -1.28 | 0.0012 | 0.0097 |
  | HMGB2 | 0.80 | 0.82 | -0.33 | 13.29 | -0.71 | 0.9704 | 0.0000 |
  | CDK1 | 0.39 | 0.72 | -2.03 | 12.24 | -1.52 | 0.0003 | 0.0000 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C13

- **n_cells**: 34
- **median_n_detected**: 874.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | PEA15 | 0.67 | 2/3 |
  | PSAT1 | 0.46 | 3/3 |
  | FABP7 | 0.41 | 3/3 |
  | TMED9 | 0.40 | 3/3 |
  | HSPA5 | 0.38 | 2/3 |
  | KHSRP | 0.36 | 2/3 |
  | MSI1 | 0.35 | 3/3 |
  | PHIP | 0.35 | 3/3 |
  | IGF2BP2 | 0.35 | 3/3 |
  | RPS18 | 0.34 | 2/3 |
  | PDIA4 | 0.33 | 3/3 |
  | PSIP1 | 0.33 | 3/3 |
  | NDRG2 | 0.32 | 3/3 |
  | ALDH6A1 | 0.32 | 3/3 |
  | RPS3 | 0.32 | 2/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | MSI1 | 1.00 | 0.66 | 3.1e-04 |
  | NFIB | 1.00 | 0.82 | 3.2e-02 |
  | HSPB1 | 1.00 | 0.82 | 3.2e-02 |
  | SF3B1 | 1.00 | 0.79 | 1.6e-02 |
  | CYB5R3 | 1.00 | 0.83 | 4.6e-02 |
  | UGDH | 1.00 | 0.64 | 2.1e-04 |
  | MAP4 | 1.00 | 0.74 | 4.6e-03 |
  | RPA1 | 1.00 | 0.75 | 4.6e-03 |
  | RCC1 | 1.00 | 0.78 | 1.1e-02 |
  | TAGLN2 | 1.00 | 0.83 | 4.6e-02 |
  | TPR | 1.00 | 0.71 | 2.1e-03 |
  | PRMT1 | 1.00 | 0.82 | 4.7e-02 |
  | TMED9 | 1.00 | 0.75 | 4.6e-03 |
  | PSPC1 | 1.00 | 0.74 | 3.2e-03 |
  | RPL23 | 1.00 | 0.76 | 6.9e-03 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | RBP1 | 0.06 | 0.40 | 1.3e-03 |
  | KIF5B | 0.21 | 0.53 | 9.2e-03 |
  | CAPZA2 | 0.06 | 0.36 | 3.4e-03 |
  | COTL1 | 0.26 | 0.56 | 1.6e-02 |
  | RPS27A | 0.32 | 0.61 | 2.3e-02 |
  | TPM4 | 0.32 | 0.61 | 3.2e-02 |
  | KLC1 | 0.06 | 0.31 | 1.9e-02 |
  | ELAVL2 | 0.15 | 0.39 | 4.7e-02 |
  | ALOX12B | 0.12 | 0.36 | 4.3e-02 |
  | CCDC183 | 0.03 | 0.24 | 2.8e-02 |
  | HELZ2 | 0.00 | 0.19 | 3.2e-02 |
  | VAT1 | 0.00 | 0.18 | 3.2e-02 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | ALB | +2.78 | 3.0e-02 |
  | PEA15 | +1.36 | 3.8e-08 |
  | TMED4 | +1.36 | 5.7e-03 |
  | PSAT1 | +1.08 | 2.6e-07 |
  | RPS2 | +0.88 | 1.3e-05 |
  | FUBP3 | +0.83 | 3.1e-03 |
  | FABP7 | +0.83 | 7.1e-09 |
  | RPS8 | +0.72 | 3.0e-07 |
  | CALU | +0.72 | 6.1e-06 |
  | PDIA4 | +0.71 | 3.2e-07 |
  | SPTBN2 | +0.70 | 2.9e-02 |
  | SEPTIN11 | +0.70 | 2.1e-04 |
  | RPL5 | +0.68 | 2.1e-08 |
  | CTBP2 | +0.66 | 9.2e-07 |
  | ALDH6A1 | +0.66 | 6.6e-05 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | SERPINB12 | -2.26 | 3.5e-02 |
  | S100A13 | -1.71 | 1.4e-02 |
  | SAMHD1 | -1.58 | 9.9e-03 |
  | NUDT3 | -1.53 | 5.1e-04 |
  | H1-1 | -1.45 | 3.4e-02 |
  | ELAVL3 | -1.30 | 9.1e-03 |
  | MTA2 | -1.16 | 1.0e-02 |
  | AK1 | -1.08 | 2.0e-02 |
  | CENPV | -0.99 | 3.2e-02 |
  | RAC2 | -0.98 | 3.4e-02 |
  | ARPC2 | -0.94 | 4.2e-02 |
  | RAP1B | -0.94 | 2.4e-02 |
  | GTF2A2 | -0.93 | 2.9e-02 |
  | COX7A2 | -0.91 | 5.0e-04 |
  | SRM | -0.89 | 3.4e-04 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | IGHG1 | +3.18 | 2.8e-05 |
  | RNMT | +2.67 | 3.8e-02 |
  | PITPNB | +2.00 | 2.5e-03 |
  | PEA15 | +1.94 | 1.5e-25 |
  | PICK1 | +1.91 | 3.4e-04 |
  | SRSF11 | +1.79 | 3.2e-04 |
  | H1-2 | +1.70 | 4.6e-02 |
  | PSMC4 | +1.70 | 5.9e-03 |
  | H1-0 | +1.61 | 1.3e-13 |
  | ADD2 | +1.56 | 4.0e-06 |
  | PSAT1 | +1.53 | 8.3e-19 |
  | QKI | +1.51 | 1.8e-09 |
  | DPY19L2 | +1.48 | 1.2e-02 |
  | TMED9 | +1.47 | 6.5e-21 |
  | NDRG2 | +1.44 | 2.0e-17 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | TERF2IP | -5.15 | 1.9e-02 |
  | RAB21 | -2.23 | 6.2e-06 |
  | C3 | -1.95 | 3.2e-02 |
  | CA2 | -1.95 | 2.0e-02 |
  | VCL | -1.87 | 4.9e-03 |
  | FKBP9 | -1.47 | 1.6e-04 |
  | C11orf54 | -1.28 | 3.4e-02 |
  | INSR | -1.07 | 3.7e-02 |
  | SMARCC2 | -1.01 | 9.1e-04 |
  | SMARCC2 | -1.01 | 9.1e-04 |
  | GPN2 | -0.87 | 2.1e-02 |
  | S100A8 | -0.82 | 3.5e-02 |
  | TF | -0.79 | 1.5e-03 |
  | JPT1 | -0.69 | 8.3e-04 |
  | MAD2L1 | -0.65 | 5.4e-03 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_MYC_TARGETS_V1 | 0.197 | +0.020 |
  | HALLMARK_UNFOLDED_PROTEIN_RESPONSE | 0.066 | +0.017 |
  | HALLMARK_G2M_CHECKPOINT | 0.095 | +0.007 |
  | HALLMARK_MYC_TARGETS_V2 | 0.029 | +0.007 |
  | HALLMARK_ESTROGEN_RESPONSE_EARLY | 0.030 | +0.005 |
  | HALLMARK_APOPTOSIS | 0.019 | +0.005 |
  | HALLMARK_E2F_TARGETS | 0.083 | +0.005 |
  | HALLMARK_APICAL_SURFACE | 0.019 | +0.004 |
  | HALLMARK_P53_PATHWAY | 0.027 | +0.004 |
  | HALLMARK_ESTROGEN_RESPONSE_LATE | 0.037 | +0.004 |
  | HALLMARK_MTORC1_SIGNALING | 0.075 | +0.004 |
  | HALLMARK_UV_RESPONSE_UP | 0.018 | +0.003 |
  | HALLMARK_UV_RESPONSE_DN | 0.012 | +0.002 |
  | HALLMARK_ANDROGEN_RESPONSE | 0.019 | +0.002 |
  | HALLMARK_IL2_STAT5_SIGNALING | 0.002 | +0.001 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 2.18 | 13.85 | -0.53 | 0.7387 | 0.0146 |
  | TUBB3 | 1.00 | 0.96 | 1.48 | 15.62 | -0.06 | 1.0000 | 0.5983 |
  | DCX | 0.79 | 0.64 | 1.02 | 12.93 | -0.49 | 0.4570 | 0.0086 |
  | STMN1 | 1.00 | 1.00 | -4.51 | 15.56 | -0.45 | 1.0000 | 0.1119 |
  | STMN2 | 0.38 | 0.53 | -0.84 | 15.05 | -0.24 | 0.4761 | 0.7905 |
  | GAP43 | 0.03 | 0.08 | -1.02 | — | — | 1.0000 | — |
  | SYN1 | 0.03 | 0.02 | 1.15 | — | — | 1.0000 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.03 | 0.03 | 0.54 | — | — | 1.0000 | — |
  | NEFL | 0.00 | 0.03 | -1.27 | — | — | 1.0000 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.88 | 0.73 | 1.33 | 12.54 | -0.55 | 0.2867 | 0.0046 |
  | UCHL1 | 1.00 | 0.92 | 2.62 | 15.66 | 0.31 | 0.4617 | 0.3041 |
  | SOX2 | 0.29 | 0.10 | 1.87 | 12.66 | 0.38 | 0.0377 | 0.0777 |
  | NES | 0.85 | 0.73 | 0.96 | 14.57 | -0.02 | 0.5865 | 0.7804 |
  | VIM | 0.32 | 0.54 | -1.25 | 17.46 | 0.08 | 0.1605 | 0.9901 |
  | PAX6 | 0.21 | 0.04 | 2.57 | 11.70 | -0.29 | 0.0193 | 0.9948 |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 1.00 | 0.95 | 1.87 | 20.01 | 0.83 | 0.9393 | 0.0000 |
  | HOPX | 0.79 | 0.28 | 3.25 | 13.62 | 0.29 | 0.0000 | 0.4998 |
  | EOMES | 0.00 | 0.09 | -2.72 | — | — | 0.4617 | — |
  | TBR1 | 0.03 | 0.11 | -1.43 | — | — | 0.7386 | — |
  | SATB2 | 0.03 | 0.02 | 1.47 | 16.64 | -0.10 | 0.9548 | 0.9724 |
  | BCL11B | 0.00 | 0.02 | -0.65 | — | — | 1.0000 | — |
  | CUX1 | 0.00 | 0.03 | -1.15 | — | — | 1.0000 | — |
  | CUX2 | 0.03 | 0.06 | -0.65 | — | — | 1.0000 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | 0.10 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.85 | 0.75 | 0.88 | 14.68 | -0.53 | 0.7078 | 0.3836 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.03 | 0.12 | -1.56 | — | — | 0.5865 | — |
  | SLC1A3 | 0.03 | 0.02 | 0.82 | — | — | 1.0000 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.03 | 0.01 | 2.67 | — | — | 0.7078 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.09 | 0.14 | -0.54 | — | — | 1.0000 | — |
  | AIF1 | 0.00 | 0.06 | -2.21 | — | — | 0.7387 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.06 | 0.11 | -0.70 | — | — | 1.0000 | — |
  | COL1A2 | 0.29 | 0.17 | 1.04 | 14.78 | 0.43 | 0.4617 | 0.5890 |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.03 | 0.03 | 0.36 | — | — | 1.0000 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.74 | 0.52 | 1.33 | 13.98 | -0.09 | 0.1161 | 0.7426 |
  | TOP2A | 0.88 | 0.58 | 2.27 | 14.39 | 0.09 | 0.0086 | 0.8376 |
  | HMGB2 | 1.00 | 0.82 | 3.97 | 13.87 | -0.07 | 0.0323 | 0.8385 |
  | CDK1 | 0.91 | 0.70 | 1.98 | 13.73 | 0.03 | 0.0671 | 0.8385 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C14

- **n_cells**: 31
- **median_n_detected**: 704.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | VIM | 0.64 | 3/3 |
  | CASP14 | 0.34 | 3/3 |
  | FBL | 0.34 | 3/3 |
  | RAB2A | 0.33 | 1/3 |
  | RPS27A | 0.31 | 2/3 |
  | H2AC21 | 0.30 | 3/3 |
  | SSR4 | 0.30 | 1/3 |
  | UBB | 0.30 | 3/3 |
  | TPM4 | 0.30 | 3/3 |
  | SMARCC2 | 0.30 | 3/3 |
  | GSTM2 | 0.29 | 3/3 |
  | RAB21 | 0.29 | 1/3 |
  | TBC1D30 | 0.28 | 3/3 |
  | PPA2 | 0.28 | 3/3 |
  | FTO | 0.28 | 2/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | RANBP1 | 1.00 | 0.80 | 2.9e-02 |
  | MCM4 | 1.00 | 0.73 | 4.7e-03 |
  | ACTBL2 | 1.00 | 0.77 | 1.2e-02 |
  | MCM2 | 1.00 | 0.77 | 1.2e-02 |
  | HMGB2 | 1.00 | 0.82 | 4.4e-02 |
  | SUPT16H | 1.00 | 0.79 | 1.9e-02 |
  | NES | 0.97 | 0.73 | 2.3e-02 |
  | MCM6 | 0.97 | 0.73 | 2.3e-02 |
  | PPA2 | 0.97 | 0.48 | 1.4e-06 |
  | CASP14 | 0.97 | 0.43 | 1.2e-07 |
  | TOP2A | 0.97 | 0.58 | 1.1e-04 |
  | CDK1 | 0.97 | 0.69 | 9.1e-03 |
  | BANF1 | 0.97 | 0.73 | 2.3e-02 |
  | MKI67 | 0.94 | 0.51 | 5.6e-05 |
  | UBE2V1 | 0.94 | 0.47 | 7.6e-06 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | RPS27A | 0.06 | 0.62 | 1.2e-07 |
  | TPM4 | 0.06 | 0.62 | 1.3e-07 |
  | RAB2A | 0.03 | 0.58 | 1.2e-07 |
  | PPP2CA | 0.06 | 0.56 | 2.4e-06 |
  | QKI | 0.23 | 0.71 | 6.7e-06 |
  | GATAD2A | 0.06 | 0.55 | 3.7e-06 |
  | VIM | 0.06 | 0.55 | 3.9e-06 |
  | DECR1 | 0.13 | 0.61 | 1.0e-05 |
  | GSTM2 | 0.03 | 0.50 | 3.9e-06 |
  | UBE2N | 0.06 | 0.53 | 7.6e-06 |
  | EEF1A1 | 0.06 | 0.53 | 7.6e-06 |
  | UBE2V2 | 0.06 | 0.53 | 7.6e-06 |
  | CALM1 | 0.06 | 0.53 | 7.6e-06 |
  | RBMX | 0.06 | 0.53 | 7.6e-06 |
  | KHDRBS1 | 0.06 | 0.52 | 8.6e-06 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | ALB | +2.91 | 1.4e-04 |
  | HSPH1 | +1.93 | 6.7e-03 |
  | TUBB4B | +1.72 | 1.5e-03 |
  | KRT5 | +1.18 | 9.8e-04 |
  | KRT2 | +0.75 | 7.6e-04 |
  | RPS2 | +0.70 | 4.3e-04 |
  | PSMA6 | +0.66 | 1.1e-03 |
  | EIF4G2 | +0.62 | 9.2e-03 |
  | KRT10 | +0.52 | 9.0e-03 |
  | CSTA | +0.50 | 3.4e-04 |
  | STMN1 | +0.50 | 2.6e-02 |
  | H4-16 | +0.50 | 1.3e-04 |
  | KRT1 | +0.47 | 6.8e-03 |
  | H1-4 | +0.46 | 2.1e-04 |
  | KRT9 | +0.45 | 3.7e-03 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | FBL | -4.31 | 8.3e-05 |
  | SRSF6 | -3.46 | 4.1e-03 |
  | HMGA2 | -1.98 | 4.8e-04 |
  | H1-1 | -1.57 | 3.7e-03 |
  | GFAP | -1.48 | 3.1e-06 |
  | HNRNPAB | -1.31 | 3.6e-03 |
  | CSTB | -1.28 | 5.4e-05 |
  | ERP29 | -1.22 | 1.4e-04 |
  | CSRP2 | -1.19 | 3.4e-02 |
  | ALDH2 | -1.18 | 1.8e-03 |
  | PPIB | -1.18 | 8.8e-10 |
  | MTA2 | -1.14 | 3.4e-04 |
  | NCBP2L | -1.12 | 4.0e-02 |
  | COX4I1 | -1.12 | 3.2e-03 |
  | DSTN | -1.11 | 4.1e-03 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | H1-2 | +2.23 | 3.6e-02 |
  | COL1A1 | +2.20 | 2.4e-02 |
  | FHL1 | +1.98 | 1.8e-02 |
  | NME2P1 | +1.96 | 1.9e-02 |
  | CARHSP1 | +1.80 | 4.0e-02 |
  | ECHS1 | +1.54 | 1.6e-02 |
  | CUL3 | +1.48 | 1.0e-02 |
  | SCARB2 | +1.47 | 3.9e-02 |
  | FTO | +1.42 | 1.9e-03 |
  | ETFB | +1.39 | 3.6e-02 |
  | UBAP2L | +1.37 | 6.6e-03 |
  | GSTO1 | +1.36 | 6.6e-03 |
  | USP14 | +1.35 | 4.1e-03 |
  | AP1B1 | +1.35 | 2.3e-02 |
  | CENPV | +1.34 | 1.5e-02 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | VIM | -2.97 | 1.6e-02 |
  | VIM | -2.97 | 1.6e-02 |
  | CA2 | -2.49 | 3.9e-02 |
  | RAB21 | -2.40 | 1.2e-02 |
  | PDLIM1 | -1.76 | 1.2e-02 |
  | SMARCC2 | -1.51 | 7.4e-03 |
  | SMARCC2 | -1.51 | 7.4e-03 |
  | DCUN1D5 | -1.02 | 1.7e-02 |
  | CCDC158 | -0.91 | 6.6e-03 |
  | SMARCC1 | -0.90 | 2.4e-02 |
  | KRT25 | -0.78 | 2.0e-02 |
  | CD1D | -0.71 | 3.2e-02 |
  | MYBPC1 | -0.54 | 1.0e-02 |
  | ERP29 | -0.49 | 2.0e-02 |
  | GFAP | -0.47 | 2.9e-02 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_E2F_TARGETS | 0.096 | +0.018 |
  | HALLMARK_G2M_CHECKPOINT | 0.101 | +0.013 |
  | HALLMARK_ESTROGEN_RESPONSE_EARLY | 0.036 | +0.011 |
  | HALLMARK_PEROXISOME | 0.036 | +0.010 |
  | HALLMARK_P53_PATHWAY | 0.032 | +0.009 |
  | HALLMARK_ESTROGEN_RESPONSE_LATE | 0.041 | +0.008 |
  | HALLMARK_MYC_TARGETS_V1 | 0.183 | +0.005 |
  | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | 0.047 | +0.004 |
  | HALLMARK_UV_RESPONSE_UP | 0.018 | +0.003 |
  | HALLMARK_MITOTIC_SPINDLE | 0.038 | +0.003 |
  | HALLMARK_CHOLESTEROL_HOMEOSTASIS | 0.034 | +0.002 |
  | HALLMARK_DNA_REPAIR | 0.009 | +0.001 |
  | HALLMARK_IL2_STAT5_SIGNALING | 0.002 | +0.001 |
  | HALLMARK_ADIPOGENESIS | 0.008 | +0.000 |
  | HALLMARK_TGF_BETA_SIGNALING | 0.002 | +0.000 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 2.04 | 14.23 | -0.12 | 0.7058 | 0.3691 |
  | TUBB3 | 1.00 | 0.96 | 1.34 | 15.51 | -0.17 | 1.0000 | 0.5364 |
  | DCX | 0.84 | 0.64 | 1.43 | 12.99 | -0.40 | 0.2136 | 0.1175 |
  | STMN1 | 1.00 | 1.00 | -4.65 | 16.43 | 0.50 | 1.0000 | 0.0263 |
  | STMN2 | 0.13 | 0.54 | -2.83 | — | — | 0.0002 | — |
  | GAP43 | 0.03 | 0.08 | -0.88 | — | — | 0.9699 | — |
  | SYN1 | 0.03 | 0.02 | 1.30 | — | — | 0.9469 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.00 | 0.03 | -1.01 | — | — | 1.0000 | — |
  | NEFL | 0.00 | 0.03 | -1.13 | — | — | 1.0000 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.77 | 0.73 | 0.24 | 12.80 | -0.25 | 1.0000 | 0.3151 |
  | UCHL1 | 1.00 | 0.92 | 2.49 | 14.90 | -0.55 | 0.5587 | 0.0004 |
  | SOX2 | 0.10 | 0.11 | -0.06 | — | — | 1.0000 | — |
  | NES | 0.97 | 0.73 | 2.92 | 14.37 | -0.24 | 0.0226 | 0.0537 |
  | VIM | 0.06 | 0.55 | -3.83 | 17.65 | -0.55 | 0.0000 | 0.3221 |
  | PAX6 | 0.13 | 0.05 | 1.71 | — | — | 0.3334 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 1.00 | 0.95 | 1.73 | 19.11 | -0.12 | 0.8761 | 0.2350 |
  | HOPX | 0.00 | 0.31 | -4.83 | — | — | 0.0008 | — |
  | EOMES | 0.19 | 0.08 | 1.57 | 12.03 | -0.39 | 0.2307 | 0.2809 |
  | TBR1 | 0.03 | 0.11 | -1.29 | — | — | 0.7031 | — |
  | SATB2 | 0.00 | 0.02 | -0.13 | — | — | 1.0000 | — |
  | BCL11B | 0.03 | 0.02 | 1.20 | — | — | 0.9699 | — |
  | CUX1 | 0.00 | 0.03 | -1.01 | — | — | 1.0000 | — |
  | CUX2 | 0.03 | 0.06 | -0.51 | — | — | 1.0000 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | 0.23 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.74 | 0.75 | -0.12 | 13.74 | -1.48 | 1.0000 | 0.0000 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.00 | 0.12 | -3.06 | — | — | 0.2453 | — |
  | SLC1A3 | 0.00 | 0.03 | -0.75 | — | — | 1.0000 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | 0.93 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.06 | 0.14 | -0.94 | — | — | 0.7675 | — |
  | AIF1 | 0.00 | 0.06 | -2.08 | — | — | 0.7058 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.06 | 0.11 | -0.56 | — | — | 1.0000 | — |
  | COL1A2 | 0.23 | 0.18 | 0.52 | 14.56 | 0.23 | 0.9617 | 0.8554 |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.03 | -1.19 | — | — | 1.0000 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.94 | 0.51 | 3.51 | 14.35 | 0.29 | 0.0001 | 0.1618 |
  | TOP2A | 0.97 | 0.58 | 3.87 | 14.77 | 0.52 | 0.0001 | 0.1623 |
  | HMGB2 | 1.00 | 0.82 | 3.83 | 13.93 | -0.01 | 0.0439 | 0.9695 |
  | CDK1 | 0.97 | 0.69 | 3.16 | 13.85 | 0.15 | 0.0091 | 0.3137 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C15

- **n_cells**: 17
- **median_n_detected**: 636.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | CRMP1 | 0.60 | 2/3 |
  | TUBB3 | 0.58 | 2/3 |
  | VIM | 0.56 | 2/3 |
  | YWHAG | 0.48 | 2/3 |
  | DPYSL3 | 0.44 | 3/3 |
  | DCX | 0.43 | 3/3 |
  | TAGLN3 | 0.36 | 3/3 |
  | UCHL1 | 0.36 | 2/3 |
  | EEF1D | 0.33 | 3/3 |
  | PTMS | 0.31 | 2/3 |
  | PDIA4 | 0.30 | 3/3 |
  | PPP1R14B | 0.27 | 3/3 |
  | ANXA5 | 0.26 | 3/3 |
  | QKI | 0.25 | 1/3 |
  | MAPK1 | 0.24 | 3/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | DCX | 1.00 | 0.64 | 3.3e-02 |
  | IRGQ | 1.00 | 0.51 | 1.3e-03 |
  | CEP170 | 1.00 | 0.34 | 5.2e-06 |
  | TAGLN3 | 1.00 | 0.22 | 3.9e-08 |
  | ELAVL3 | 1.00 | 0.45 | 2.1e-04 |
  | ELAVL2 | 1.00 | 0.37 | 1.8e-05 |
  | AK1 | 0.94 | 0.47 | 4.4e-03 |
  | DCLK1 | 0.94 | 0.38 | 3.2e-04 |
  | MACROH2A2 | 0.94 | 0.46 | 4.4e-03 |
  | SH3BGRL | 0.82 | 0.99 | 1.4e-02 |
  | EEF1D | 0.82 | 1.00 | 1.0e-02 |
  | ANXA5 | 0.76 | 0.99 | 4.5e-03 |
  | MAPK1 | 0.76 | 0.98 | 1.2e-02 |
  | RBP1 | 0.76 | 0.37 | 4.4e-02 |
  | RPS8 | 0.71 | 0.96 | 2.4e-02 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | GLUD1 | 0.18 | 0.89 | 1.1e-07 |
  | QKI | 0.00 | 0.71 | 1.0e-06 |
  | PCNA | 0.18 | 0.87 | 6.6e-07 |
  | FASN | 0.18 | 0.85 | 1.9e-06 |
  | UGDH | 0.00 | 0.67 | 5.2e-06 |
  | SHMT2 | 0.12 | 0.78 | 6.0e-06 |
  | PSME1 | 0.18 | 0.83 | 5.2e-06 |
  | MCM7 | 0.18 | 0.77 | 1.1e-04 |
  | MCM4 | 0.18 | 0.76 | 2.0e-04 |
  | MCM6 | 0.18 | 0.76 | 2.0e-04 |
  | PRDX4 | 0.06 | 0.64 | 2.1e-04 |
  | NES | 0.18 | 0.75 | 2.2e-04 |
  | RCN1 | 0.24 | 0.80 | 1.9e-04 |
  | SYNE1 | 0.29 | 0.86 | 9.2e-05 |
  | LGALS3 | 0.06 | 0.62 | 3.2e-04 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | TUBB3 | +2.24 | 9.4e-06 |
  | CRMP1 | +2.21 | 1.8e-05 |
  | DCX | +2.20 | 1.3e-04 |
  | YWHAG | +1.68 | 1.8e-05 |
  | DPYSL3 | +1.60 | 1.1e-04 |
  | PPP1R14B | +1.57 | 1.3e-03 |
  | DDAH2 | +1.52 | 8.5e-05 |
  | STMN1 | +1.42 | 9.1e-04 |
  | UCHL1 | +1.28 | 3.3e-04 |
  | JPT1 | +1.22 | 5.8e-04 |
  | PAFAH1B3 | +1.20 | 3.5e-04 |
  | PTMS | +1.10 | 1.8e-05 |
  | H1-10 | +1.08 | 2.8e-04 |
  | ENO2 | +1.04 | 9.6e-04 |
  | MAP2 | +1.02 | 1.4e-03 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | VIM | -4.02 | 5.8e-04 |
  | VIM | -3.28 | 2.1e-02 |
  | RPS14 | -1.71 | 2.1e-02 |
  | DBI | -1.66 | 4.1e-05 |
  | SEPTIN5 | -1.65 | 4.1e-02 |
  | CROCC | -1.61 | 4.6e-02 |
  | HSPB1 | -1.59 | 1.4e-04 |
  | PEA15 | -1.54 | 1.8e-05 |
  | NDRG2 | -1.50 | 1.6e-04 |
  | PDIA4 | -1.43 | 1.1e-04 |
  | ZBTB18 | -1.39 | 3.6e-02 |
  | RPN1 | -1.38 | 1.9e-02 |
  | FABP7 | -1.37 | 4.7e-05 |
  | ANXA5 | -1.37 | 3.4e-05 |
  | EEF1D | -1.31 | 2.3e-06 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | DDX58 | +4.19 | 5.1e-06 |
  | CAT | +3.61 | 1.4e-05 |
  | SELENBP1 | +2.91 | 7.0e-04 |
  | TUBB3 | +2.34 | 8.7e-16 |
  | CRMP1 | +2.27 | 2.3e-18 |
  | DPYSL3 | +2.24 | 1.2e-15 |
  | DCX | +1.96 | 2.2e-11 |
  | UCHL1 | +1.85 | 6.6e-17 |
  | YWHAG | +1.82 | 7.7e-19 |
  | SH3BGRL3 | +1.71 | 5.1e-05 |
  | PFN2 | +1.70 | 2.1e-06 |
  | RTN1 | +1.68 | 9.8e-04 |
  | H1-0 | +1.67 | 8.1e-10 |
  | RBP1 | +1.60 | 1.5e-06 |
  | PPP1R14B | +1.56 | 8.9e-11 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | VIM | -3.75 | 4.9e-08 |
  | VIM | -3.75 | 4.9e-08 |
  | PCNA | -2.60 | 8.7e-03 |
  | NDRG1 | -2.15 | 2.0e-04 |
  | C1QBP | -1.95 | 7.4e-05 |
  | JPT2 | -1.78 | 1.6e-03 |
  | FABP7 | -1.71 | 1.6e-03 |
  | NAP1L1 | -1.66 | 3.6e-03 |
  | HSPB1 | -1.63 | 3.1e-06 |
  | CDK1 | -1.58 | 3.2e-04 |
  | CDK1 | -1.58 | 3.2e-04 |
  | MCM6 | -1.50 | 5.6e-06 |
  | PHGDH | -1.31 | 1.3e-05 |
  | EEF1D | -1.21 | 6.4e-08 |
  | FABP5 | -1.19 | 4.9e-09 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_HEDGEHOG_SIGNALING | 0.134 | +0.044 |
  | HALLMARK_PANCREAS_BETA_CELLS | 0.040 | +0.034 |
  | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | 0.057 | +0.014 |
  | HALLMARK_ADIPOGENESIS | 0.015 | +0.007 |
  | HALLMARK_HEME_METABOLISM | 0.021 | +0.006 |
  | HALLMARK_SPERMATOGENESIS | 0.013 | +0.005 |
  | HALLMARK_XENOBIOTIC_METABOLISM | 0.010 | +0.005 |
  | HALLMARK_BILE_ACID_METABOLISM | 0.009 | +0.004 |
  | HALLMARK_G2M_CHECKPOINT | 0.092 | +0.004 |
  | HALLMARK_ESTROGEN_RESPONSE_EARLY | 0.029 | +0.003 |
  | HALLMARK_TGF_BETA_SIGNALING | 0.004 | +0.003 |
  | HALLMARK_IL2_STAT5_SIGNALING | 0.002 | +0.001 |
  | HALLMARK_ANGIOGENESIS | 0.000 | -0.000 |
  | HALLMARK_NOTCH_SIGNALING | 0.000 | -0.000 |
  | HALLMARK_WNT_BETA_CATENIN_SIGNALING | 0.000 | -0.000 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 1.17 | 15.32 | 1.02 | 1.0000 | 0.0014 |
  | TUBB3 | 1.00 | 0.96 | 0.47 | 17.88 | 2.24 | 1.0000 | 0.0000 |
  | DCX | 1.00 | 0.64 | 4.29 | 15.54 | 2.20 | 0.0333 | 0.0001 |
  | STMN1 | 1.00 | 1.00 | -5.52 | 17.37 | 1.42 | 1.0000 | 0.0009 |
  | STMN2 | 0.59 | 0.52 | 0.36 | 16.14 | 0.88 | 1.0000 | 0.1053 |
  | GAP43 | 0.00 | 0.08 | -1.66 | — | — | 1.0000 | — |
  | SYN1 | 0.00 | 0.02 | 0.45 | — | — | 1.0000 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.06 | 0.03 | 1.59 | — | — | 1.0000 | — |
  | NEFL | 0.00 | 0.03 | -0.26 | — | — | 1.0000 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 1.00 | 0.73 | 3.69 | 14.04 | 1.04 | 0.1470 | 0.0010 |
  | UCHL1 | 1.00 | 0.92 | 1.61 | 16.63 | 1.28 | 1.0000 | 0.0003 |
  | SOX2 | 0.00 | 0.11 | -2.19 | — | — | 0.9598 | — |
  | NES | 0.18 | 0.75 | -3.64 | — | — | 0.0002 | — |
  | VIM | 0.47 | 0.53 | -0.33 | 13.41 | -4.02 | 1.0000 | 0.0006 |
  | PAX6 | 0.00 | 0.05 | -0.93 | — | — | 1.0000 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 1.00 | 0.95 | 0.86 | 17.87 | -1.37 | 1.0000 | 0.0000 |
  | HOPX | 0.00 | 0.31 | -3.95 | — | — | 0.0586 | — |
  | EOMES | 0.00 | 0.08 | -1.71 | — | — | 1.0000 | — |
  | TBR1 | 0.35 | 0.10 | 2.36 | 12.15 | -1.53 | 0.0955 | 0.1320 |
  | SATB2 | 0.18 | 0.01 | 4.19 | 16.60 | -0.11 | 0.0484 | 0.9683 |
  | BCL11B | 0.06 | 0.02 | 2.12 | — | — | 1.0000 | — |
  | CUX1 | 0.00 | 0.03 | -0.14 | — | — | 1.0000 | — |
  | CUX2 | 0.00 | 0.06 | -1.29 | — | — | 1.0000 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | 1.11 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.47 | 0.76 | -1.79 | 15.33 | 0.14 | 0.2289 | 0.6463 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.00 | 0.11 | -2.19 | — | — | 0.9598 | — |
  | SLC1A3 | 0.00 | 0.02 | 0.13 | — | — | 1.0000 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.06 | 0.01 | 3.72 | — | — | 0.6914 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.06 | 0.14 | -0.83 | — | — | 1.0000 | — |
  | AIF1 | 0.00 | 0.06 | -1.20 | — | — | 1.0000 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.18 | 0.11 | 1.00 | — | — | 1.0000 | — |
  | COL1A2 | 0.18 | 0.18 | 0.16 | — | — | 1.0000 | — |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.03 | -0.31 | — | — | 1.0000 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.00 | 0.54 | -5.34 | — | — | 0.0003 | — |
  | TOP2A | 0.41 | 0.60 | -1.08 | 13.65 | -0.69 | 0.7680 | 0.2343 |
  | HMGB2 | 1.00 | 0.82 | 2.95 | 14.20 | 0.27 | 0.4396 | 0.9365 |
  | CDK1 | 0.24 | 0.71 | -2.91 | — | — | 0.0044 | — |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C16

- **n_cells**: 39
- **median_n_detected**: 688.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | VIM | 0.46 | 3/3 |
  | UBB | 0.40 | 3/3 |
  | TPM4 | 0.36 | 3/3 |
  | ARHGDIA | 0.35 | 3/3 |
  | GNA11 | 0.35 | 2/3 |
  | TF | 0.33 | 2/3 |
  | FBL | 0.33 | 2/3 |
  | PPP2CA | 0.32 | 3/3 |
  | RPS27A | 0.32 | 2/3 |
  | H2AC21 | 0.32 | 3/3 |
  | NDRG2 | 0.30 | 3/3 |
  | ALB | 0.30 | 2/3 |
  | RAB11B | 0.29 | 3/3 |
  | RAB2A | 0.29 | 3/3 |
  | GSTM2 | 0.29 | 3/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | PCNA | 1.00 | 0.84 | 2.4e-02 |
  | CALU | 1.00 | 0.84 | 2.4e-02 |
  | SMARCA5 | 1.00 | 0.83 | 1.6e-02 |
  | TAGLN2 | 1.00 | 0.83 | 1.6e-02 |
  | RCN1 | 1.00 | 0.78 | 2.7e-03 |
  | SUPT16H | 1.00 | 0.79 | 2.7e-03 |
  | MCM4 | 0.97 | 0.73 | 3.6e-03 |
  | FLNA | 0.97 | 0.76 | 9.0e-03 |
  | MCM7 | 0.97 | 0.75 | 5.7e-03 |
  | MCM6 | 0.97 | 0.73 | 3.6e-03 |
  | PSAT1 | 0.97 | 0.80 | 3.0e-02 |
  | MCM2 | 0.97 | 0.77 | 1.4e-02 |
  | ACTBL2 | 0.97 | 0.77 | 1.4e-02 |
  | MCM5 | 0.97 | 0.75 | 5.8e-03 |
  | MCM3 | 0.97 | 0.75 | 5.8e-03 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | RPS27A | 0.08 | 0.62 | 2.7e-09 |
  | TPM4 | 0.08 | 0.62 | 2.8e-09 |
  | COX4I1 | 0.10 | 0.62 | 2.5e-08 |
  | ARPC5 | 0.18 | 0.69 | 4.2e-08 |
  | RAB2A | 0.08 | 0.58 | 3.5e-08 |
  | COTL1 | 0.08 | 0.57 | 4.2e-08 |
  | PPP2CA | 0.08 | 0.57 | 5.8e-08 |
  | ARHGDIA | 0.38 | 0.87 | 5.3e-09 |
  | KIF5B | 0.05 | 0.54 | 5.1e-08 |
  | GATAD2A | 0.08 | 0.55 | 1.5e-07 |
  | VIM | 0.08 | 0.55 | 1.6e-07 |
  | UBE2V2 | 0.08 | 0.53 | 4.3e-07 |
  | EEF1A1 | 0.08 | 0.53 | 4.3e-07 |
  | UBE2N | 0.08 | 0.53 | 4.3e-07 |
  | RBMX | 0.08 | 0.53 | 4.3e-07 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | ALB | +2.99 | 1.2e-05 |
  | EIF3D | +1.81 | 4.6e-02 |
  | TUBB4B | +1.73 | 2.5e-03 |
  | DCPS | +1.26 | 4.4e-03 |
  | KRT5 | +1.00 | 3.7e-04 |
  | PDCD6 | +0.91 | 3.7e-02 |
  | RPS2 | +0.76 | 2.6e-05 |
  | KRT2 | +0.72 | 1.0e-04 |
  | HK1 | +0.71 | 4.0e-02 |
  | EIF4G2 | +0.66 | 6.0e-03 |
  | PSME1 | +0.64 | 2.1e-04 |
  | EIF3C | +0.56 | 4.4e-02 |
  | PSMA6 | +0.53 | 4.6e-03 |
  | MCM3 | +0.52 | 5.0e-07 |
  | IPO7 | +0.51 | 3.3e-02 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | FBL | -4.24 | 1.3e-04 |
  | SRSF6 | -3.35 | 7.0e-03 |
  | H1-1 | -1.87 | 7.2e-04 |
  | MTA2 | -1.69 | 3.8e-06 |
  | SNRPG | -1.62 | 2.7e-02 |
  | LAP3 | -1.54 | 1.7e-03 |
  | TMSB10 | -1.53 | 1.4e-05 |
  | GFAP | -1.53 | 1.8e-07 |
  | HNRNPAB | -1.45 | 2.1e-04 |
  | ARPC4 | -1.41 | 1.1e-05 |
  | ARHGDIA | -1.39 | 3.4e-04 |
  | NCBP2L | -1.39 | 1.3e-03 |
  | RTN1 | -1.36 | 8.8e-03 |
  | PPP1R14B | -1.36 | 1.1e-04 |
  | TUFM | -1.35 | 1.1e-08 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | GNA11 | +2.66 | 2.3e-05 |
  | PTPRZ1 | +2.42 | 3.8e-02 |
  | UQCRFS1 | +1.99 | 1.2e-03 |
  | PITPNB | +1.63 | 8.0e-03 |
  | FLNB | +1.62 | 2.0e-03 |
  | CARHSP1 | +1.51 | 4.7e-04 |
  | ATP6V1A | +1.50 | 5.4e-03 |
  | UBQLN1 | +1.48 | 1.1e-02 |
  | CUL3 | +1.47 | 2.2e-02 |
  | SEPTIN3 | +1.38 | 2.1e-02 |
  | ATP5PB | +1.38 | 1.6e-02 |
  | SUPT4H1 | +1.35 | 3.8e-02 |
  | EIF3J | +1.31 | 1.9e-02 |
  | TRA2B | +1.30 | 7.7e-03 |
  | PHF5A | +1.28 | 3.5e-02 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | VIM | -2.42 | 1.4e-02 |
  | VIM | -2.42 | 1.4e-02 |
  | TF | -2.29 | 3.3e-06 |
  | VCL | -1.88 | 2.1e-02 |
  | SOGA1 | -1.67 | 3.2e-04 |
  | SERPINH1 | -1.66 | 9.8e-03 |
  | SMARCC2 | -1.50 | 5.7e-03 |
  | SMARCC2 | -1.50 | 5.7e-03 |
  | CTPS1 | -1.02 | 1.5e-03 |
  | STMN1 | -0.90 | 1.7e-08 |
  | JPT1 | -0.87 | 3.4e-05 |
  | MKI67 | -0.82 | 1.8e-04 |
  | CALM1 | -0.80 | 3.5e-02 |
  | CALM1 | -0.80 | 3.5e-02 |
  | H3-3B | -0.75 | 3.5e-02 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_MYC_TARGETS_V1 | 0.194 | +0.016 |
  | HALLMARK_E2F_TARGETS | 0.091 | +0.013 |
  | HALLMARK_G2M_CHECKPOINT | 0.098 | +0.010 |
  | HALLMARK_ESTROGEN_RESPONSE_EARLY | 0.034 | +0.009 |
  | HALLMARK_P53_PATHWAY | 0.030 | +0.007 |
  | HALLMARK_UV_RESPONSE_DN | 0.013 | +0.003 |
  | HALLMARK_ESTROGEN_RESPONSE_LATE | 0.036 | +0.003 |
  | HALLMARK_UV_RESPONSE_UP | 0.017 | +0.002 |
  | HALLMARK_UNFOLDED_PROTEIN_RESPONSE | 0.052 | +0.002 |
  | HALLMARK_PEROXISOME | 0.028 | +0.001 |
  | HALLMARK_ALLOGRAFT_REJECTION | 0.041 | +0.001 |
  | HALLMARK_DNA_REPAIR | 0.009 | +0.001 |
  | HALLMARK_CHOLESTEROL_HOMEOSTASIS | 0.033 | +0.000 |
  | HALLMARK_IL2_STAT5_SIGNALING | 0.002 | +0.000 |
  | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITI | 0.028 | -0.000 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 2.39 | 13.38 | -1.03 | 0.4721 | 0.0000 |
  | TUBB3 | 1.00 | 0.96 | 1.68 | 15.07 | -0.64 | 0.7757 | 0.0002 |
  | DCX | 0.56 | 0.65 | -0.55 | 12.31 | -1.11 | 0.6798 | 0.0000 |
  | STMN1 | 1.00 | 1.00 | -4.30 | 15.27 | -0.77 | 1.0000 | 0.0000 |
  | STMN2 | 0.10 | 0.54 | -3.23 | — | — | 0.0000 | — |
  | GAP43 | 0.03 | 0.08 | -1.23 | — | — | 0.7316 | — |
  | SYN1 | 0.03 | 0.02 | 0.94 | — | — | 0.9299 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.00 | 0.03 | -1.36 | — | — | 0.9705 | — |
  | NEFL | 0.00 | 0.03 | -1.47 | — | — | 0.9705 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.56 | 0.74 | -1.18 | 12.22 | -0.85 | 0.1011 | 0.0001 |
  | UCHL1 | 1.00 | 0.92 | 2.83 | 14.93 | -0.55 | 0.2686 | 0.0000 |
  | SOX2 | 0.28 | 0.10 | 1.79 | 12.55 | 0.25 | 0.0242 | 0.1756 |
  | NES | 0.87 | 0.73 | 1.20 | 14.44 | -0.17 | 0.2578 | 0.0993 |
  | VIM | 0.08 | 0.55 | -3.68 | 17.81 | -0.33 | 0.0000 | 0.3352 |
  | PAX6 | 0.15 | 0.04 | 2.03 | 12.24 | 0.39 | 0.0735 | 0.7360 |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 1.00 | 0.95 | 2.07 | 19.40 | 0.21 | 0.6036 | 0.0829 |
  | HOPX | 0.41 | 0.29 | 0.75 | 13.07 | -0.38 | 0.4563 | 0.2026 |
  | EOMES | 0.03 | 0.09 | -1.28 | — | — | 0.6036 | — |
  | TBR1 | 0.03 | 0.11 | -1.64 | — | — | 0.4742 | — |
  | SATB2 | 0.00 | 0.02 | -0.47 | — | — | 1.0000 | — |
  | BCL11B | 0.00 | 0.02 | -0.86 | — | — | 1.0000 | — |
  | CUX1 | 0.00 | 0.03 | -1.36 | — | — | 0.9705 | — |
  | CUX2 | 0.03 | 0.07 | -0.86 | — | — | 0.8915 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | -0.11 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.79 | 0.75 | 0.33 | 13.70 | -1.53 | 0.9633 | 0.0000 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.00 | 0.12 | -3.41 | — | — | 0.1035 | — |
  | SLC1A3 | 0.03 | 0.02 | 0.61 | — | — | 1.0000 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | 0.59 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.00 | 0.14 | -3.73 | — | — | 0.0523 | — |
  | AIF1 | 0.00 | 0.06 | -2.42 | — | — | 0.4721 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.03 | 0.11 | -1.71 | — | — | 0.3581 | — |
  | COL1A2 | 0.28 | 0.17 | 0.96 | 14.58 | 0.25 | 0.3125 | 0.5346 |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.03 | -1.53 | — | — | 0.9708 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.74 | 0.51 | 1.41 | 13.60 | -0.48 | 0.0434 | 0.0015 |
  | TOP2A | 0.95 | 0.58 | 3.44 | 13.46 | -0.93 | 0.0000 | 0.0005 |
  | HMGB2 | 0.95 | 0.82 | 1.76 | 13.63 | -0.35 | 0.1651 | 0.0251 |
  | CDK1 | 0.87 | 0.70 | 1.45 | 13.32 | -0.45 | 0.1157 | 0.0066 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C17

- **n_cells**: 40
- **median_n_detected**: 558.5
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | FBL | 0.35 | 3/3 |
  | RAB1A | 0.34 | 3/3 |
  | TAF15 | 0.34 | 3/3 |
  | CFL1 | 0.34 | 3/3 |
  | TAGLN2 | 0.33 | 3/3 |
  | H1-4 | 0.33 | 2/3 |
  | IDH2 | 0.32 | 3/3 |
  | PCBP2 | 0.32 | 3/3 |
  | CCT6A | 0.32 | 3/3 |
  | HS3ST3A1 | 0.32 | 2/3 |
  | PABPC1 | 0.32 | 3/3 |
  | H1-3 | 0.31 | 2/3 |
  | KRT1 | 0.31 | 3/3 |
  | TUBA4A | 0.31 | 3/3 |
  | RAB6A | 0.31 | 2/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | KRT5 | 1.00 | 0.79 | 2.1e-03 |
  | NES | 1.00 | 0.73 | 8.5e-05 |
  | HSPB1 | 1.00 | 0.82 | 5.3e-03 |
  | MCM3 | 0.97 | 0.75 | 2.9e-03 |
  | MCM6 | 0.97 | 0.73 | 1.1e-03 |
  | ACTBL2 | 0.97 | 0.77 | 7.2e-03 |
  | MCM7 | 0.97 | 0.75 | 1.8e-03 |
  | MCM2 | 0.97 | 0.77 | 4.6e-03 |
  | RALY | 0.95 | 1.00 | 3.2e-02 |
  | MARCKS | 0.95 | 0.79 | 4.0e-02 |
  | HADH | 0.95 | 0.79 | 4.0e-02 |
  | MCM4 | 0.95 | 0.73 | 7.4e-03 |
  | RCN1 | 0.95 | 0.78 | 3.9e-02 |
  | ILF3 | 0.95 | 1.00 | 3.2e-02 |
  | UBA1 | 0.95 | 1.00 | 3.2e-02 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | HDAC2 | 0.07 | 0.90 | 1.8e-28 |
  | H1-4 | 0.15 | 0.95 | 5.1e-31 |
  | TAF15 | 0.15 | 0.95 | 5.1e-31 |
  | PTBP1 | 0.15 | 0.94 | 1.8e-29 |
  | KRT1 | 0.17 | 0.95 | 4.0e-30 |
  | CFL1 | 0.17 | 0.95 | 4.0e-30 |
  | SRSF2 | 0.15 | 0.93 | 1.7e-27 |
  | PABPC1 | 0.17 | 0.95 | 6.4e-30 |
  | IDH2 | 0.17 | 0.95 | 1.0e-29 |
  | SLC25A6 | 0.17 | 0.95 | 1.0e-29 |
  | RAB1A | 0.17 | 0.95 | 2.7e-29 |
  | PCBP2 | 0.17 | 0.95 | 2.7e-29 |
  | RPS6 | 0.10 | 0.87 | 7.2e-24 |
  | RAB6A | 0.17 | 0.93 | 6.3e-26 |
  | SMARCA5 | 0.15 | 0.87 | 2.9e-21 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | USP7 | +3.90 | 1.2e-10 |
  | TGM1 | +3.28 | 2.5e-03 |
  | ALB | +2.65 | 1.5e-02 |
  | HSPA4L | +2.48 | 1.5e-04 |
  | PDCD6IP | +2.44 | 2.8e-03 |
  | HMGA2 | +2.17 | 1.1e-07 |
  | DLD | +2.03 | 4.2e-05 |
  | NSF | +1.73 | 4.6e-05 |
  | AZGP1 | +1.51 | 5.4e-04 |
  | RAB6A | +1.46 | 3.4e-02 |
  | CMBL | +1.42 | 1.7e-04 |
  | COL11A2 | +1.40 | 3.5e-08 |
  | DENND1A | +1.31 | 4.7e-14 |
  | HARS1 | +1.23 | 1.0e-03 |
  | TBCB | +1.14 | 1.1e-06 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | FBL | -4.77 | 6.2e-12 |
  | CCDC183 | -3.22 | 4.4e-04 |
  | H2BC10 | -2.43 | 2.7e-15 |
  | PSAT1 | -2.25 | 3.6e-08 |
  | DLST | -2.10 | 2.6e-07 |
  | YWHAZ | -2.10 | 8.9e-13 |
  | HNRNPAB | -1.98 | 1.2e-09 |
  | MKI67 | -1.87 | 9.4e-05 |
  | HADHB | -1.81 | 2.4e-04 |
  | CSTB | -1.64 | 2.2e-06 |
  | RCC1 | -1.61 | 1.6e-07 |
  | JPT1 | -1.60 | 9.8e-04 |
  | FABP7 | -1.58 | 1.2e-07 |
  | MYEF2 | -1.53 | 9.3e-07 |
  | TUBB2A | -1.49 | 3.6e-06 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | HBB | +4.11 | 7.5e-03 |
  | SEPTIN5 | +2.76 | 1.4e-02 |
  | FHL1 | +2.71 | 1.4e-03 |
  | SEC11A | +2.41 | 2.6e-02 |
  | UQCRFS1 | +2.40 | 1.2e-06 |
  | NUCKS1 | +2.40 | 4.5e-02 |
  | GNA13 | +2.33 | 1.3e-08 |
  | GNA11 | +2.33 | 6.6e-07 |
  | TAGLN2 | +2.27 | 1.1e-22 |
  | EML5 | +2.27 | 9.5e-07 |
  | IMPA1 | +2.24 | 1.1e-03 |
  | ITGB6 | +2.24 | 4.7e-02 |
  | PITPNB | +2.22 | 3.3e-02 |
  | RHEB | +2.16 | 4.1e-02 |
  | RUFY3 | +2.12 | 3.5e-02 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | AP2A2 | -4.47 | 1.2e-02 |
  | HELZ2 | -2.55 | 3.0e-05 |
  | TOP2A | -2.26 | 3.8e-05 |
  | TOP2A | -2.26 | 3.8e-05 |
  | C11orf54 | -2.15 | 5.5e-10 |
  | PIP | -1.61 | 2.5e-02 |
  | HEXB | -1.52 | 5.6e-04 |
  | KPNA2 | -1.41 | 2.5e-02 |
  | VCL | -1.36 | 3.0e-02 |
  | ARHGAP17 | -1.33 | 3.1e-03 |
  | GTF3C5 | -1.26 | 7.3e-03 |
  | RAB21 | -1.24 | 3.8e-03 |
  | MACF1 | -1.06 | 1.0e-02 |
  | MACF1 | -1.06 | 1.0e-02 |
  | GCNT3 | -1.05 | 1.3e-02 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_INTERFERON_ALPHA_RESPONSE | 0.023 | +0.015 |
  | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITI | 0.036 | +0.009 |
  | HALLMARK_INTERFERON_GAMMA_RESPONSE | 0.015 | +0.008 |
  | HALLMARK_PEROXISOME | 0.035 | +0.008 |
  | HALLMARK_FATTY_ACID_METABOLISM | 0.051 | +0.007 |
  | HALLMARK_UV_RESPONSE_UP | 0.020 | +0.005 |
  | HALLMARK_UV_RESPONSE_DN | 0.015 | +0.005 |
  | HALLMARK_SPERMATOGENESIS | 0.012 | +0.005 |
  | HALLMARK_MYC_TARGETS_V2 | 0.026 | +0.004 |
  | HALLMARK_GLYCOLYSIS | 0.069 | +0.004 |
  | HALLMARK_MTORC1_SIGNALING | 0.075 | +0.003 |
  | HALLMARK_XENOBIOTIC_METABOLISM | 0.008 | +0.003 |
  | HALLMARK_ALLOGRAFT_REJECTION | 0.042 | +0.002 |
  | HALLMARK_APICAL_JUNCTION | 0.037 | +0.002 |
  | HALLMARK_HYPOXIA | 0.048 | +0.002 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 0.97 | 0.94 | 0.77 | 13.31 | -1.09 | 1.0000 | 0.0000 |
  | TUBB3 | 1.00 | 0.96 | 1.72 | 15.90 | 0.25 | 0.7346 | 0.4712 |
  | DCX | 0.28 | 0.67 | -2.37 | 12.88 | -0.52 | 0.0000 | 0.2585 |
  | STMN1 | 1.00 | 1.00 | -4.27 | 15.63 | -0.37 | 1.0000 | 0.0250 |
  | STMN2 | 0.47 | 0.53 | -0.29 | 14.72 | -0.58 | 0.9218 | 0.1143 |
  | GAP43 | 0.05 | 0.08 | -0.47 | — | — | 1.0000 | — |
  | SYN1 | 0.00 | 0.02 | -0.81 | — | — | 1.0000 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.00 | 0.03 | -1.39 | — | — | 0.9218 | — |
  | NEFL | 0.00 | 0.03 | -1.51 | — | — | 0.9218 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.50 | 0.75 | -1.57 | 12.47 | -0.57 | 0.0080 | 0.0157 |
  | UCHL1 | 0.97 | 0.92 | 1.22 | 14.94 | -0.53 | 0.6994 | 0.0001 |
  | SOX2 | 0.30 | 0.10 | 1.93 | 11.68 | -0.73 | 0.0051 | 0.0311 |
  | NES | 1.00 | 0.73 | 4.94 | 15.50 | 0.96 | 0.0001 | 0.0000 |
  | VIM | 0.93 | 0.51 | 3.37 | 18.14 | 0.84 | 0.0000 | 0.0003 |
  | PAX6 | 0.00 | 0.05 | -2.19 | — | — | 0.5559 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 0.97 | 0.95 | 0.45 | 17.67 | -1.58 | 1.0000 | 0.0000 |
  | HOPX | 0.12 | 0.31 | -1.52 | 12.06 | -1.38 | 0.0559 | 0.2030 |
  | EOMES | 0.00 | 0.09 | -2.96 | — | — | 0.2169 | — |
  | TBR1 | 0.03 | 0.11 | -1.68 | — | — | 0.3120 | — |
  | SATB2 | 0.00 | 0.02 | -0.51 | — | — | 1.0000 | — |
  | BCL11B | 0.00 | 0.02 | -0.89 | — | — | 1.0000 | — |
  | CUX1 | 0.35 | 0.01 | 5.33 | 13.41 | 0.50 | 0.0000 | 0.4702 |
  | CUX2 | 0.03 | 0.07 | -0.90 | — | — | 0.8475 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | -0.14 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.25 | 0.78 | -3.32 | 15.17 | -0.02 | 0.0000 | 0.6530 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.05 | 0.12 | -1.02 | — | — | 0.6266 | — |
  | SLC1A3 | 0.00 | 0.03 | -1.13 | — | — | 0.9218 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | 0.55 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.05 | 0.14 | -1.35 | — | — | 0.3988 | — |
  | AIF1 | 0.00 | 0.06 | -2.46 | — | — | 0.4172 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.25 | 0.10 | 1.60 | 13.09 | -0.80 | 0.0353 | 0.1539 |
  | COL1A2 | 0.38 | 0.17 | 1.60 | 14.64 | 0.34 | 0.0121 | 0.3624 |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.03 | 0.03 | 0.11 | — | — | 1.0000 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.30 | 0.54 | -1.40 | 12.21 | -1.87 | 0.0260 | 0.0001 |
  | TOP2A | 0.07 | 0.62 | -4.15 | 13.43 | 0.04 | 0.0000 | 0.8136 |
  | HMGB2 | 0.15 | 0.86 | -4.99 | 14.00 | 0.06 | 0.0000 | 0.5777 |
  | CDK1 | 0.15 | 0.73 | -3.86 | 13.59 | -0.12 | 0.0000 | 0.6331 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C18

- **n_cells**: 23
- **median_n_detected**: 831.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | RRM1 | 0.35 | 3/3 |
  | HLA-B | 0.34 | 2/3 |
  | DUT | 0.34 | 3/3 |
  | POTEJ | 0.34 | 3/3 |
  | ALB | 0.34 | 2/3 |
  | H2AC21 | 0.33 | 3/3 |
  | NEDD4L | 0.31 | 3/3 |
  | TPM4 | 0.30 | 3/3 |
  | PPP2CA | 0.29 | 3/3 |
  | SEMA3C | 0.28 | 3/3 |
  | VIM | 0.28 | 3/3 |
  | NFU1 | 0.27 | 3/3 |
  | MBNL1 | 0.27 | 3/3 |
  | MYH10 | 0.27 | 3/3 |
  | HSPA1A | 0.26 | 3/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | DUT | 1.00 | 0.52 | 5.4e-05 |
  | MCM7 | 1.00 | 0.75 | 4.0e-02 |
  | MCM3 | 1.00 | 0.76 | 4.0e-02 |
  | RRM1 | 1.00 | 0.56 | 3.8e-04 |
  | QDPR | 1.00 | 0.70 | 1.5e-02 |
  | ALDH7A1 | 1.00 | 0.72 | 2.4e-02 |
  | ALDH6A1 | 1.00 | 0.74 | 4.0e-02 |
  | NES | 1.00 | 0.73 | 2.4e-02 |
  | MCM6 | 1.00 | 0.74 | 2.5e-02 |
  | ISYNA1 | 1.00 | 0.72 | 2.4e-02 |
  | MCM4 | 1.00 | 0.74 | 4.0e-02 |
  | MCM5 | 1.00 | 0.76 | 4.0e-02 |
  | MAP4 | 1.00 | 0.75 | 4.0e-02 |
  | TOP2A | 1.00 | 0.59 | 6.9e-04 |
  | JPT2 | 1.00 | 0.58 | 4.3e-04 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | RPS27A | 0.09 | 0.61 | 5.3e-05 |
  | TPM4 | 0.09 | 0.61 | 5.4e-05 |
  | GATAD2A | 0.04 | 0.55 | 1.2e-04 |
  | RAB2A | 0.09 | 0.57 | 3.8e-04 |
  | PPP2CA | 0.09 | 0.56 | 4.3e-04 |
  | VIM | 0.09 | 0.54 | 9.9e-04 |
  | UBE2N | 0.09 | 0.52 | 9.9e-04 |
  | EEF1A1 | 0.09 | 0.52 | 9.9e-04 |
  | UBE2V2 | 0.09 | 0.52 | 9.9e-04 |
  | RBMX | 0.09 | 0.52 | 1.0e-03 |
  | CALM1 | 0.09 | 0.52 | 1.0e-03 |
  | KHDRBS1 | 0.09 | 0.52 | 1.1e-03 |
  | RAB5C | 0.09 | 0.52 | 1.1e-03 |
  | HSPA1A | 0.09 | 0.52 | 1.1e-03 |
  | MAGOHB | 0.09 | 0.51 | 2.2e-03 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | ALB | +2.90 | 1.6e-03 |
  | TUBB4B | +2.43 | 8.6e-03 |
  | DYNLT1 | +1.50 | 4.4e-02 |
  | SLC25A5 | +1.30 | 4.2e-02 |
  | KRT5 | +1.18 | 1.7e-02 |
  | KRT2 | +0.99 | 3.3e-03 |
  | PSMA6 | +0.98 | 7.8e-04 |
  | EIF4G2 | +0.87 | 4.3e-03 |
  | PSME1 | +0.83 | 9.7e-03 |
  | H2AC21 | +0.81 | 9.9e-03 |
  | KRT1 | +0.79 | 1.6e-02 |
  | STMN1 | +0.79 | 1.5e-02 |
  | RPS2 | +0.78 | 9.1e-04 |
  | CTPS1 | +0.74 | 4.2e-02 |
  | EIF3C | +0.71 | 3.4e-02 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | FBL | -4.08 | 4.9e-02 |
  | SRSF6 | -2.85 | 1.2e-02 |
  | PFKM | -2.84 | 2.1e-02 |
  | NDUFB4 | -1.86 | 4.2e-02 |
  | PTMA | -1.49 | 2.5e-02 |
  | PDXP | -1.49 | 3.7e-02 |
  | H1-1 | -1.46 | 9.7e-03 |
  | ARHGDIA | -1.42 | 3.3e-03 |
  | ERP29 | -1.40 | 2.0e-02 |
  | GFAP | -1.28 | 2.3e-02 |
  | ARPC2 | -1.27 | 3.9e-02 |
  | PSMA2 | -1.24 | 3.3e-02 |
  | CORO1A | -1.22 | 1.1e-02 |
  | VARS1 | -1.16 | 1.5e-02 |
  | CSTB | -1.13 | 1.5e-02 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | HLA-B | +3.35 | 5.8e-06 |
  | DLD | +2.65 | 3.6e-02 |
  | HOPX | +2.00 | 1.3e-03 |
  | PITPNB | +1.40 | 1.3e-02 |
  | UBAP2L | +1.38 | 6.3e-03 |
  | COL1A1 | +1.37 | 2.8e-02 |
  | LSM7 | +1.33 | 6.3e-03 |
  | UBA6 | +1.32 | 2.1e-02 |
  | RBBP9 | +1.29 | 9.1e-03 |
  | UQCRC2 | +1.25 | 2.0e-03 |
  | BAZ1B | +1.23 | 4.7e-02 |
  | LZIC | +1.23 | 2.1e-03 |
  | LSM6 | +1.21 | 6.3e-03 |
  | BLVRA | +1.17 | 3.2e-04 |
  | TXNRD1 | +1.15 | 4.0e-07 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | VCL | -2.05 | 3.9e-03 |
  | CA2 | -2.03 | 5.5e-03 |
  | MFN1 | -1.34 | 2.2e-02 |
  | CAMSAP3 | -1.23 | 4.6e-02 |
  | PDLIM1 | -1.23 | 1.9e-02 |
  | CDSN | -1.19 | 2.9e-02 |
  | SMARCC2 | -1.18 | 3.8e-03 |
  | SMARCC2 | -1.18 | 3.8e-03 |
  | IGHG1 | -0.90 | 4.3e-02 |
  | LSM5 | -0.75 | 3.9e-02 |
  | AZGP1 | -0.58 | 2.6e-02 |
  | ERP29 | -0.56 | 4.3e-02 |
  | IMMT | -0.49 | 1.9e-02 |
  | PIP | -0.46 | 1.8e-02 |
  | MCM7 | -0.17 | 4.8e-02 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_PEROXISOME | 0.040 | +0.013 |
  | HALLMARK_ESTROGEN_RESPONSE_EARLY | 0.034 | +0.008 |
  | HALLMARK_E2F_TARGETS | 0.086 | +0.008 |
  | HALLMARK_APICAL_SURFACE | 0.022 | +0.007 |
  | HALLMARK_ESTROGEN_RESPONSE_LATE | 0.040 | +0.007 |
  | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | 0.049 | +0.006 |
  | HALLMARK_UV_RESPONSE_DN | 0.015 | +0.005 |
  | HALLMARK_G2M_CHECKPOINT | 0.093 | +0.005 |
  | HALLMARK_P53_PATHWAY | 0.028 | +0.004 |
  | HALLMARK_GLYCOLYSIS | 0.070 | +0.004 |
  | HALLMARK_FATTY_ACID_METABOLISM | 0.048 | +0.004 |
  | HALLMARK_OXIDATIVE_PHOSPHORYLATION | 0.028 | +0.003 |
  | HALLMARK_MTORC1_SIGNALING | 0.074 | +0.002 |
  | HALLMARK_MITOTIC_SPINDLE | 0.037 | +0.002 |
  | HALLMARK_CHOLESTEROL_HOMEOSTASIS | 0.034 | +0.002 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 1.61 | 14.55 | 0.24 | 0.9170 | 0.4387 |
  | TUBB3 | 1.00 | 0.96 | 0.90 | 15.90 | 0.22 | 1.0000 | 0.7638 |
  | DCX | 0.91 | 0.64 | 2.26 | 13.29 | -0.10 | 0.0894 | 0.5899 |
  | STMN1 | 1.00 | 1.00 | -5.08 | 16.72 | 0.79 | 1.0000 | 0.0148 |
  | STMN2 | 0.13 | 0.53 | -2.75 | — | — | 0.0053 | — |
  | GAP43 | 0.00 | 0.08 | -2.10 | — | — | 0.7724 | — |
  | SYN1 | 0.00 | 0.02 | 0.01 | — | — | 1.0000 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.00 | 0.03 | -0.58 | — | — | 1.0000 | — |
  | NEFL | 0.00 | 0.03 | -0.69 | — | — | 1.0000 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.96 | 0.73 | 2.47 | 12.86 | -0.17 | 0.1475 | 0.8123 |
  | UCHL1 | 1.00 | 0.92 | 2.05 | 15.33 | -0.06 | 0.7724 | 0.7329 |
  | SOX2 | 0.22 | 0.11 | 1.27 | 12.26 | -0.11 | 0.6610 | 0.9252 |
  | NES | 1.00 | 0.73 | 4.11 | 14.70 | 0.11 | 0.0243 | 0.5918 |
  | VIM | 0.09 | 0.54 | -3.34 | 18.60 | 0.63 | 0.0010 | 0.1018 |
  | PAX6 | 0.13 | 0.05 | 1.75 | — | — | 0.5102 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 1.00 | 0.95 | 1.29 | 19.42 | 0.23 | 1.0000 | 0.2252 |
  | HOPX | 0.09 | 0.31 | -1.92 | — | — | 0.2041 | — |
  | EOMES | 0.35 | 0.08 | 2.74 | 11.90 | -0.55 | 0.0081 | 0.2112 |
  | TBR1 | 0.00 | 0.11 | -2.50 | — | — | 0.6456 | — |
  | SATB2 | 0.00 | 0.02 | 0.30 | — | — | 1.0000 | — |
  | BCL11B | 0.00 | 0.02 | -0.08 | — | — | 1.0000 | — |
  | CUX1 | 0.00 | 0.03 | -0.58 | — | — | 1.0000 | — |
  | CUX2 | 0.00 | 0.07 | -1.73 | — | — | 0.9170 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | 0.67 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.65 | 0.75 | -0.74 | 13.93 | -1.28 | 0.9127 | 0.0232 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.00 | 0.12 | -2.63 | — | — | 0.5065 | — |
  | SLC1A3 | 0.04 | 0.02 | 1.41 | — | — | 0.9789 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.04 | 0.01 | 3.26 | — | — | 0.6456 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.13 | 0.14 | 0.10 | — | — | 1.0000 | — |
  | AIF1 | 0.00 | 0.06 | -1.64 | — | — | 0.9170 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.17 | 0.11 | 0.94 | — | — | 0.8536 | — |
  | COL1A2 | 0.48 | 0.17 | 2.18 | 14.63 | 0.31 | 0.0163 | 0.3805 |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.03 | -0.75 | — | — | 1.0000 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.96 | 0.51 | 3.83 | 14.44 | 0.39 | 0.0006 | 0.1305 |
  | TOP2A | 1.00 | 0.59 | 5.06 | 14.75 | 0.45 | 0.0007 | 0.0823 |
  | HMGB2 | 1.00 | 0.82 | 3.39 | 14.01 | 0.08 | 0.2041 | 1.0000 |
  | CDK1 | 0.96 | 0.70 | 2.71 | 14.06 | 0.39 | 0.0661 | 0.0157 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C19

- **n_cells**: 18
- **median_n_detected**: 883.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | KRT5 | 0.43 | 2/3 |
  | KRT78 | 0.39 | 3/3 |
  | TBCA | 0.38 | 2/3 |
  | SRSF6 | 0.37 | 2/3 |
  | RPL22 | 0.33 | 2/3 |
  | KRT14 | 0.31 | 2/3 |
  | AP2B1 | 0.31 | 3/3 |
  | UBB | 0.29 | 3/3 |
  | S100A9 | 0.29 | 3/3 |
  | HSPA1A | 0.28 | 3/3 |
  | TUBA4A | 0.28 | 3/3 |
  | VIM | 0.27 | 2/3 |
  | RAB2A | 0.27 | 3/3 |
  | ALB | 0.27 | 2/3 |
  | CAT | 0.27 | 2/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | CLTA | 1.00 | 0.68 | 3.4e-02 |
  | ECI1 | 1.00 | 0.63 | 1.1e-02 |
  | AZGP1 | 1.00 | 0.54 | 2.6e-03 |
  | CASP14 | 0.94 | 0.43 | 1.2e-03 |
  | S100A9 | 0.94 | 0.36 | 2.6e-04 |
  | EEF1A1 | 0.94 | 0.48 | 2.8e-03 |
  | EIF5A | 0.94 | 0.41 | 6.9e-04 |
  | CALM1 | 0.94 | 0.47 | 2.8e-03 |
  | KRT13 | 0.94 | 0.47 | 2.8e-03 |
  | RAB5A | 0.94 | 0.47 | 2.8e-03 |
  | RBMX | 0.94 | 0.48 | 2.8e-03 |
  | UBE2N | 0.94 | 0.48 | 2.8e-03 |
  | MAGOH | 0.94 | 0.45 | 2.6e-03 |
  | RAB18 | 0.94 | 0.56 | 2.7e-02 |
  | TOP2A | 0.94 | 0.59 | 3.3e-02 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | RPS27A | 0.06 | 0.61 | 4.0e-04 |
  | TPM4 | 0.06 | 0.61 | 4.0e-04 |
  | RAB2A | 0.06 | 0.57 | 1.2e-03 |
  | PPP2CA | 0.06 | 0.55 | 1.7e-03 |
  | GATAD2A | 0.06 | 0.54 | 2.7e-03 |
  | VIM | 0.06 | 0.54 | 2.7e-03 |
  | STMN2 | 0.06 | 0.53 | 2.8e-03 |
  | EEF1A1 | 0.06 | 0.52 | 2.8e-03 |
  | UBE2V2 | 0.06 | 0.52 | 2.8e-03 |
  | UBE2N | 0.06 | 0.52 | 2.8e-03 |
  | CALM1 | 0.06 | 0.52 | 2.8e-03 |
  | RBMX | 0.06 | 0.52 | 2.8e-03 |
  | MBNL1 | 0.00 | 0.46 | 2.6e-03 |
  | KHDRBS1 | 0.06 | 0.52 | 2.9e-03 |
  | RAB5C | 0.06 | 0.52 | 2.9e-03 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | KRT78 | +4.27 | 1.7e-02 |
  | ALB | +2.85 | 4.6e-03 |
  | KRT5 | +2.27 | 4.3e-04 |
  | HSPH1 | +1.69 | 2.5e-02 |
  | SLC25A5 | +1.64 | 3.8e-02 |
  | KRT14 | +1.45 | 2.7e-04 |
  | KRT2 | +1.32 | 4.3e-04 |
  | KRT1 | +1.20 | 4.9e-04 |
  | HADHA | +1.19 | 3.6e-02 |
  | KRT3 | +1.17 | 1.4e-02 |
  | FLG2 | +1.13 | 1.3e-02 |
  | STMN1 | +1.10 | 9.9e-03 |
  | CDSN | +0.99 | 4.1e-02 |
  | EIF3C | +0.96 | 1.3e-02 |
  | RAB35 | +0.87 | 1.0e-02 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | SRSF6 | -3.13 | 2.5e-03 |
  | PFKM | -2.67 | 1.9e-02 |
  | NDUFB4 | -1.88 | 1.4e-02 |
  | TAGLN3 | -1.78 | 2.4e-02 |
  | ERP29 | -1.54 | 2.5e-03 |
  | ACTR2 | -1.53 | 2.7e-02 |
  | RPS14 | -1.52 | 3.1e-02 |
  | ARPC2 | -1.40 | 2.7e-02 |
  | ZNF428 | -1.40 | 3.2e-02 |
  | LAP3 | -1.39 | 2.6e-02 |
  | GFAP | -1.31 | 2.7e-03 |
  | ARPC3 | -1.27 | 3.1e-02 |
  | AGO1 | -1.27 | 2.7e-02 |
  | NUDT3 | -1.24 | 4.3e-02 |
  | NDUFS8 | -1.20 | 4.2e-02 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | PSMC4 | +2.99 | 3.5e-03 |
  | CAT | +2.69 | 7.4e-03 |
  | LAMP1 | +2.44 | 9.3e-03 |
  | ARG1 | +2.34 | 6.8e-03 |
  | KRT3 | +2.31 | 3.2e-03 |
  | KRT3 | +2.31 | 3.2e-03 |
  | MUCL1 | +2.13 | 4.2e-02 |
  | GPX1 | +1.87 | 2.1e-02 |
  | PITPNB | +1.86 | 8.0e-03 |
  | ARCN1 | +1.72 | 4.8e-03 |
  | KRT5 | +1.72 | 2.0e-04 |
  | IGHG1 | +1.55 | 1.1e-02 |
  | PGM2 | +1.52 | 2.2e-02 |
  | KRT14 | +1.51 | 2.6e-04 |
  | FHL1 | +1.48 | 2.3e-02 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | CTSV | -2.44 | 2.8e-02 |
  | KRT77 | -2.35 | 1.6e-02 |
  | HDLBP | -2.12 | 2.0e-02 |
  | VCL | -1.87 | 9.3e-03 |
  | SERPINH1 | -1.81 | 2.5e-02 |
  | CAMSAP3 | -1.34 | 1.2e-06 |
  | ANXA11 | -1.16 | 2.6e-03 |
  | SMARCC2 | -0.95 | 3.8e-02 |
  | SMARCC2 | -0.95 | 3.8e-02 |
  | C11orf96 | -0.74 | 8.7e-03 |
  | ERP29 | -0.63 | 1.2e-04 |
  | PCNA | -0.32 | 1.3e-02 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | 0.057 | +0.014 |
  | HALLMARK_PEROXISOME | 0.040 | +0.013 |
  | HALLMARK_ESTROGEN_RESPONSE_LATE | 0.043 | +0.009 |
  | HALLMARK_P53_PATHWAY | 0.030 | +0.007 |
  | HALLMARK_ESTROGEN_RESPONSE_EARLY | 0.032 | +0.007 |
  | HALLMARK_APICAL_SURFACE | 0.021 | +0.006 |
  | HALLMARK_E2F_TARGETS | 0.084 | +0.006 |
  | HALLMARK_G2M_CHECKPOINT | 0.093 | +0.005 |
  | HALLMARK_BILE_ACID_METABOLISM | 0.009 | +0.004 |
  | HALLMARK_HEME_METABOLISM | 0.017 | +0.003 |
  | HALLMARK_KRAS_SIGNALING_UP | 0.004 | +0.003 |
  | HALLMARK_IL2_STAT5_SIGNALING | 0.003 | +0.002 |
  | HALLMARK_ADIPOGENESIS | 0.009 | +0.001 |
  | HALLMARK_UV_RESPONSE_DN | 0.011 | +0.001 |
  | HALLMARK_GLYCOLYSIS | 0.066 | +0.001 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 0.94 | 0.94 | -0.45 | 14.76 | 0.45 | 1.0000 | 0.2478 |
  | TUBB3 | 0.94 | 0.96 | -1.17 | 16.17 | 0.52 | 1.0000 | 0.0848 |
  | DCX | 0.94 | 0.64 | 2.70 | 13.54 | 0.16 | 0.0843 | 0.6942 |
  | STMN1 | 1.00 | 1.00 | -5.44 | 17.05 | 1.10 | 1.0000 | 0.0099 |
  | STMN2 | 0.06 | 0.53 | -3.74 | — | — | 0.0028 | — |
  | GAP43 | 0.00 | 0.08 | -1.74 | — | — | 1.0000 | — |
  | SYN1 | 0.00 | 0.02 | 0.36 | — | — | 1.0000 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.06 | 0.03 | 1.51 | — | — | 1.0000 | — |
  | NEFL | 0.00 | 0.03 | -0.34 | — | — | 1.0000 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.94 | 0.73 | 2.10 | 13.07 | 0.04 | 0.3864 | 0.6666 |
  | UCHL1 | 0.94 | 0.92 | 0.00 | 15.48 | 0.11 | 1.0000 | 0.7596 |
  | SOX2 | 0.17 | 0.11 | 0.85 | — | — | 1.0000 | — |
  | NES | 0.94 | 0.73 | 2.08 | 14.60 | 0.01 | 0.3864 | 0.9865 |
  | VIM | 0.06 | 0.54 | -3.77 | 18.34 | 0.31 | 0.0027 | 0.4760 |
  | PAX6 | 0.00 | 0.05 | -1.01 | — | — | 1.0000 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 0.94 | 0.95 | -0.77 | 19.34 | 0.14 | 1.0000 | 0.3765 |
  | HOPX | 0.11 | 0.30 | -1.52 | — | — | 0.6129 | — |
  | EOMES | 0.33 | 0.08 | 2.62 | 11.67 | -0.75 | 0.0488 | 0.1644 |
  | TBR1 | 0.00 | 0.11 | -2.14 | — | — | 0.8498 | — |
  | SATB2 | 0.00 | 0.02 | 0.66 | — | — | 1.0000 | — |
  | BCL11B | 0.06 | 0.02 | 2.03 | — | — | 1.0000 | — |
  | CUX1 | 0.00 | 0.03 | -0.22 | — | — | 1.0000 | — |
  | CUX2 | 0.00 | 0.06 | -1.37 | — | — | 1.0000 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | 1.03 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.78 | 0.75 | 0.11 | 13.90 | -1.31 | 1.0000 | 0.0027 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.00 | 0.11 | -2.27 | — | — | 0.8498 | — |
  | SLC1A3 | 0.00 | 0.02 | 0.04 | — | — | 1.0000 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | 1.73 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.00 | 0.14 | -2.59 | — | — | 0.6905 | — |
  | AIF1 | 0.06 | 0.06 | 0.41 | — | — | 1.0000 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.17 | 0.11 | 0.90 | — | — | 1.0000 | — |
  | COL1A2 | 0.11 | 0.18 | -0.52 | — | — | 1.0000 | — |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.03 | -0.39 | — | — | 1.0000 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.83 | 0.52 | 2.04 | 14.89 | 0.84 | 0.1082 | 0.0060 |
  | TOP2A | 0.94 | 0.59 | 3.02 | 15.03 | 0.78 | 0.0325 | 0.0112 |
  | HMGB2 | 1.00 | 0.82 | 3.04 | 13.85 | -0.10 | 0.3864 | 0.4460 |
  | CDK1 | 0.94 | 0.70 | 2.33 | 14.41 | 0.73 | 0.2779 | 0.0027 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C2

- **n_cells**: 47
- **median_n_detected**: 712.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | GSN | 0.71 | 3/3 |
  | COTL1 | 0.64 | 3/3 |
  | CTSB | 0.64 | 3/3 |
  | FABP5 | 0.56 | 2/3 |
  | CTSD | 0.55 | 3/3 |
  | PFN1 | 0.46 | 2/3 |
  | TUBB2B | 0.43 | 3/3 |
  | GPX1 | 0.43 | 3/3 |
  | CORO1A | 0.35 | 3/3 |
  | RBBP4 | 0.34 | 2/3 |
  | AIF1 | 0.33 | 1/3 |
  | RNASET2 | 0.33 | 1/3 |
  | DPYSL2 | 0.33 | 3/3 |
  | NASP | 0.31 | 3/3 |
  | FABP7 | 0.30 | 3/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | NDUFA4 | 1.00 | 0.60 | 1.6e-09 |
  | RNASET2 | 1.00 | 0.00 | 2.2e-71 |
  | WDR1 | 1.00 | 0.89 | 2.1e-02 |
  | TLN1 | 1.00 | 0.31 | 8.7e-22 |
  | NAGK | 1.00 | 0.63 | 9.9e-09 |
  | SAMHD1 | 1.00 | 0.10 | 3.8e-41 |
  | ARPC5 | 1.00 | 0.64 | 2.4e-08 |
  | GSN | 1.00 | 0.35 | 1.6e-19 |
  | PGM2 | 1.00 | 0.39 | 6.2e-18 |
  | CNDP2 | 1.00 | 0.84 | 3.7e-03 |
  | COTL1 | 1.00 | 0.52 | 2.8e-12 |
  | SH3BGRL3 | 1.00 | 0.45 | 2.7e-15 |
  | CTSZ | 1.00 | 0.06 | 5.0e-48 |
  | S100A11 | 1.00 | 0.05 | 1.7e-51 |
  | SERPINB9 | 1.00 | 0.04 | 2.5e-53 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | MAP2 | 0.09 | 0.99 | 7.8e-59 |
  | PHGDH | 0.06 | 0.96 | 4.3e-47 |
  | DDAH1 | 0.06 | 0.96 | 2.4e-46 |
  | PSIP1 | 0.09 | 0.97 | 2.7e-48 |
  | SEPTIN11 | 0.09 | 0.96 | 1.7e-45 |
  | PPA1 | 0.13 | 0.99 | 1.1e-55 |
  | PFN2 | 0.13 | 0.98 | 2.2e-47 |
  | CNN3 | 0.04 | 0.89 | 7.6e-35 |
  | DPYSL5 | 0.11 | 0.94 | 3.5e-38 |
  | NFIB | 0.04 | 0.88 | 2.5e-33 |
  | RCN1 | 0.02 | 0.84 | 3.4e-31 |
  | PSAT1 | 0.04 | 0.85 | 1.7e-30 |
  | DDX6 | 0.11 | 0.90 | 2.4e-31 |
  | LANCL1 | 0.04 | 0.83 | 4.6e-28 |
  | ARPC1A | 0.04 | 0.82 | 1.4e-27 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | CTSB | +4.30 | 7.3e-23 |
  | GSN | +4.06 | 4.6e-25 |
  | COTL1 | +3.50 | 7.0e-26 |
  | GPX1 | +3.41 | 1.9e-24 |
  | CORO1A | +2.84 | 2.5e-21 |
  | CTSD | +2.84 | 7.0e-26 |
  | CTSZ | +2.72 | 2.2e-13 |
  | KCTD12 | +2.54 | 1.5e-08 |
  | LCP1 | +2.23 | 3.8e-14 |
  | PYCARD | +2.21 | 1.7e-06 |
  | RAP1B | +2.18 | 1.8e-11 |
  | SAMHD1 | +2.17 | 7.2e-16 |
  | CSTB | +2.13 | 3.1e-14 |
  | FSCN1 | +2.01 | 4.6e-20 |
  | PFN1 | +2.00 | 7.5e-26 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | HMOX1 | -6.41 | 5.5e-05 |
  | FABP7 | -4.66 | 1.0e-10 |
  | PLS1 | -4.23 | 6.2e-03 |
  | TUBB2B | -3.46 | 7.7e-20 |
  | FABP5 | -3.24 | 7.5e-26 |
  | TMEM109 | -2.73 | 7.8e-06 |
  | HSPB1 | -2.71 | 9.4e-06 |
  | BLVRB | -2.58 | 2.1e-04 |
  | DPYSL2 | -2.41 | 3.3e-25 |
  | ALB | -2.36 | 4.2e-03 |
  | PTMS | -2.30 | 1.1e-21 |
  | PCNA | -2.25 | 7.9e-11 |
  | TUBB2A | -2.24 | 3.5e-03 |
  | TUBA1A | -2.18 | 8.1e-25 |
  | DBI | -1.97 | 5.2e-20 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | CTSB | +4.93 | 9.8e-24 |
  | GSN | +4.21 | 2.3e-38 |
  | GPX1 | +3.95 | 1.4e-13 |
  | COTL1 | +3.66 | 2.3e-42 |
  | CTSD | +3.39 | 7.5e-44 |
  | HPRT1 | +2.91 | 1.9e-03 |
  | PSAP | +2.85 | 1.8e-09 |
  | CORO1A | +2.77 | 2.6e-10 |
  | RAP1B | +2.66 | 2.5e-12 |
  | GNAI2 | +2.39 | 2.3e-21 |
  | LAMP1 | +2.38 | 7.4e-06 |
  | LGALS3BP | +2.36 | 1.2e-04 |
  | COL1A1 | +2.35 | 2.2e-02 |
  | CNDP2 | +2.29 | 6.7e-24 |
  | GNA13 | +2.28 | 6.2e-14 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | FABP7 | -3.47 | 5.8e-15 |
  | FABP5 | -3.18 | 1.4e-43 |
  | TUBB2B | -3.07 | 7.0e-27 |
  | STMN1 | -2.41 | 2.0e-27 |
  | PCNA | -2.36 | 1.7e-27 |
  | HSPB1 | -2.20 | 2.4e-11 |
  | VCL | -2.18 | 3.0e-03 |
  | DPYSL2 | -2.17 | 2.5e-29 |
  | PTMS | -2.10 | 2.1e-27 |
  | RBBP4 | -2.00 | 3.1e-42 |
  | HMGB3 | -2.00 | 1.0e-05 |
  | JPT1 | -1.99 | 3.3e-11 |
  | VIM | -1.99 | 1.3e-02 |
  | VIM | -1.99 | 1.3e-02 |
  | NASP | -1.99 | 2.8e-41 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_COMPLEMENT | 0.077 | +0.055 |
  | HALLMARK_COAGULATION | 0.043 | +0.039 |
  | HALLMARK_OXIDATIVE_PHOSPHORYLATION | 0.060 | +0.036 |
  | HALLMARK_APICAL_JUNCTION | 0.066 | +0.032 |
  | HALLMARK_PI3K_AKT_MTOR_SIGNALING | 0.068 | +0.029 |
  | HALLMARK_MTORC1_SIGNALING | 0.098 | +0.027 |
  | HALLMARK_MYOGENESIS | 0.037 | +0.020 |
  | HALLMARK_HYPOXIA | 0.064 | +0.020 |
  | HALLMARK_PEROXISOME | 0.042 | +0.016 |
  | HALLMARK_APOPTOSIS | 0.029 | +0.015 |
  | HALLMARK_KRAS_SIGNALING_UP | 0.016 | +0.015 |
  | HALLMARK_MITOTIC_SPINDLE | 0.048 | +0.014 |
  | HALLMARK_CHOLESTEROL_HOMEOSTASIS | 0.043 | +0.011 |
  | HALLMARK_INTERFERON_GAMMA_RESPONSE | 0.016 | +0.010 |
  | HALLMARK_BILE_ACID_METABOLISM | 0.015 | +0.009 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 0.09 | 0.99 | -10.40 | — | — | 0.0000 | — |
  | TUBB3 | 0.40 | 1.00 | -8.82 | 14.20 | -1.49 | 0.0000 | 0.0000 |
  | DCX | 0.09 | 0.68 | -4.39 | — | — | 0.0000 | — |
  | STMN1 | 1.00 | 1.00 | -4.02 | 14.23 | -1.82 | 1.0000 | 0.0000 |
  | STMN2 | 0.38 | 0.53 | -0.85 | 13.69 | -1.64 | 0.1229 | 0.0000 |
  | GAP43 | 0.00 | 0.09 | -3.16 | — | — | 0.0700 | — |
  | SYN1 | 0.00 | 0.02 | -1.05 | — | — | 1.0000 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.00 | 0.03 | -1.64 | — | — | 0.6148 | — |
  | NEFL | 0.00 | 0.03 | -1.75 | — | — | 0.6148 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.09 | 0.78 | -5.06 | — | — | 0.0000 | — |
  | UCHL1 | 0.23 | 0.96 | -6.37 | 13.71 | -1.68 | 0.0000 | 0.0000 |
  | SOX2 | 0.02 | 0.12 | -2.06 | — | — | 0.1240 | — |
  | NES | 0.19 | 0.77 | -3.78 | 14.80 | 0.21 | 0.0000 | 0.8670 |
  | VIM | 0.11 | 0.55 | -3.26 | 14.60 | -2.81 | 0.0000 | 0.1762 |
  | PAX6 | 0.00 | 0.05 | -2.43 | — | — | 0.3152 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 0.47 | 0.98 | -5.89 | 14.58 | -4.66 | 0.0000 | 0.0000 |
  | HOPX | 0.02 | 0.32 | -3.84 | — | — | 0.0000 | — |
  | EOMES | 0.00 | 0.09 | -3.21 | — | — | 0.0700 | — |
  | TBR1 | 0.00 | 0.11 | -3.56 | — | — | 0.0326 | — |
  | SATB2 | 0.00 | 0.02 | -0.76 | 17.43 | 0.75 | 1.0000 | 0.0279 |
  | BCL11B | 0.00 | 0.02 | -1.14 | — | — | 0.8340 | — |
  | CUX1 | 0.02 | 0.03 | 0.04 | — | — | 1.0000 | — |
  | CUX2 | 0.00 | 0.07 | -2.79 | — | — | 0.1526 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | -0.39 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.70 | 0.75 | -0.39 | 15.45 | 0.29 | 0.7128 | 0.2530 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.02 | 0.12 | -2.06 | — | — | 0.1240 | — |
  | SLC1A3 | 0.00 | 0.03 | -1.37 | — | — | 0.8340 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | 0.31 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.04 | 0.14 | -1.60 | — | — | 0.1210 | — |
  | AIF1 | 1.00 | 0.00 | 14.84 | — | — | 0.0000 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.04 | 0.11 | -1.22 | — | — | 0.4097 | — |
  | COL1A2 | 0.02 | 0.19 | -2.83 | — | — | 0.0053 | — |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.03 | -1.81 | — | — | 0.6148 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.04 | 0.55 | -4.50 | — | — | 0.0000 | — |
  | TOP2A | 0.06 | 0.63 | -4.43 | — | — | 0.0000 | — |
  | HMGB2 | 0.51 | 0.84 | -2.35 | 12.74 | -1.25 | 0.0000 | 0.0000 |
  | CDK1 | 0.45 | 0.72 | -1.66 | 12.24 | -1.51 | 0.0009 | 0.0000 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C20

- **n_cells**: 40
- **median_n_detected**: 655.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | CNN3 | 0.44 | 3/3 |
  | PCBP2 | 0.42 | 3/3 |
  | PABPC1 | 0.38 | 3/3 |
  | PEA15 | 0.36 | 2/3 |
  | PSAT1 | 0.35 | 3/3 |
  | FBL | 0.35 | 3/3 |
  | PTBP1 | 0.35 | 3/3 |
  | RAB6A | 0.35 | 3/3 |
  | PRDX6 | 0.34 | 2/3 |
  | IDH2 | 0.34 | 3/3 |
  | TAF15 | 0.34 | 3/3 |
  | HMGB2 | 0.33 | 3/3 |
  | CFL1 | 0.33 | 3/3 |
  | RPS27 | 0.33 | 2/3 |
  | SLC25A5 | 0.32 | 3/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | RCN1 | 1.00 | 0.78 | 1.1e-03 |
  | DPYSL3 | 1.00 | 0.86 | 4.4e-02 |
  | FASN | 1.00 | 0.83 | 7.9e-03 |
  | ARHGDIA | 1.00 | 0.84 | 2.0e-02 |
  | GLUD1 | 1.00 | 0.87 | 4.4e-02 |
  | TAGLN2 | 1.00 | 0.83 | 1.2e-02 |
  | ACTBL2 | 1.00 | 0.77 | 1.1e-03 |
  | CNN3 | 1.00 | 0.83 | 1.2e-02 |
  | MCM5 | 0.97 | 0.75 | 4.2e-03 |
  | NES | 0.97 | 0.73 | 1.5e-03 |
  | MCM3 | 0.97 | 0.75 | 4.2e-03 |
  | HADH | 0.97 | 0.79 | 1.6e-02 |
  | FLNA | 0.97 | 0.76 | 4.3e-03 |
  | ARL3 | 0.97 | 0.62 | 6.1e-06 |
  | MTHFD1 | 0.97 | 0.81 | 3.8e-02 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | CFL1 | 0.23 | 0.95 | 2.0e-26 |
  | KRT1 | 0.23 | 0.95 | 2.0e-26 |
  | H1-4 | 0.23 | 0.95 | 2.9e-26 |
  | PABPC1 | 0.23 | 0.95 | 2.9e-26 |
  | TAF15 | 0.23 | 0.95 | 2.9e-26 |
  | SLC25A6 | 0.23 | 0.95 | 4.5e-26 |
  | IDH2 | 0.23 | 0.95 | 4.5e-26 |
  | RAB1A | 0.23 | 0.94 | 1.2e-25 |
  | PCBP2 | 0.23 | 0.94 | 1.2e-25 |
  | PTBP1 | 0.23 | 0.94 | 3.7e-24 |
  | RAB6A | 0.23 | 0.92 | 1.4e-22 |
  | SRSF2 | 0.23 | 0.92 | 2.0e-22 |
  | ELMOD1 | 0.12 | 0.79 | 1.6e-16 |
  | HDAC2 | 0.23 | 0.89 | 1.1e-18 |
  | RPS27 | 0.23 | 0.88 | 5.5e-18 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | ALB | +2.73 | 8.5e-03 |
  | USP7 | +2.42 | 1.4e-05 |
  | KRT80 | +2.11 | 2.9e-04 |
  | EWSR1 | +1.84 | 7.3e-07 |
  | HMGA2 | +1.60 | 1.0e-03 |
  | VIM | +1.24 | 1.2e-07 |
  | TBCB | +1.18 | 4.2e-09 |
  | TNIK | +1.16 | 1.6e-02 |
  | ADA | +1.14 | 2.9e-03 |
  | CACYBP | +1.11 | 2.1e-06 |
  | PGRMC1 | +1.07 | 4.6e-07 |
  | RAB6A | +1.07 | 5.2e-03 |
  | CNN3 | +1.07 | 1.2e-11 |
  | DENND1A | +1.04 | 2.5e-08 |
  | FHL1 | +1.03 | 9.6e-03 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | FBL | -4.63 | 1.2e-07 |
  | CALD1 | -3.26 | 3.1e-02 |
  | DNM2 | -2.51 | 1.5e-04 |
  | H2BC10 | -2.08 | 3.0e-10 |
  | GOT1 | -1.98 | 5.9e-03 |
  | HADHB | -1.81 | 9.6e-06 |
  | UBA6 | -1.79 | 2.6e-03 |
  | YWHAZ | -1.70 | 3.4e-08 |
  | HNRNPAB | -1.62 | 7.8e-05 |
  | DLST | -1.62 | 1.4e-05 |
  | GPX1 | -1.61 | 2.0e-05 |
  | MYEF2 | -1.43 | 7.3e-07 |
  | CSTB | -1.39 | 9.9e-05 |
  | DYNC1LI1 | -1.35 | 1.2e-02 |
  | RCC1 | -1.19 | 1.6e-05 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | SNX27 | +2.84 | 9.3e-05 |
  | PSAT1 | +2.34 | 1.5e-28 |
  | TAGLN2 | +2.22 | 5.3e-26 |
  | FHL1 | +2.16 | 4.3e-05 |
  | PEA15 | +2.11 | 1.6e-30 |
  | HUWE1 | +2.06 | 2.6e-02 |
  | CNN3 | +2.01 | 6.5e-30 |
  | SMAD4 | +1.96 | 2.5e-03 |
  | RPS28 | +1.94 | 6.6e-14 |
  | PTPRZ1 | +1.87 | 3.6e-03 |
  | RHEB | +1.85 | 9.7e-03 |
  | NDRG2 | +1.84 | 9.2e-17 |
  | NCAM1 | +1.77 | 2.0e-04 |
  | DDAH1 | +1.70 | 5.0e-27 |
  | PRDX6 | +1.63 | 5.9e-40 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | HELZ2 | -2.71 | 2.8e-04 |
  | C11orf54 | -2.47 | 1.2e-12 |
  | KRT77 | -2.43 | 2.0e-02 |
  | CNP | -2.10 | 1.7e-03 |
  | S100A9 | -1.88 | 1.7e-02 |
  | SNCB | -1.70 | 2.6e-02 |
  | VCL | -1.55 | 4.6e-03 |
  | HARS1 | -1.35 | 3.6e-02 |
  | FKBP9 | -1.26 | 3.0e-02 |
  | ALOX12B | -1.10 | 8.0e-03 |
  | S100A8 | -1.02 | 4.9e-02 |
  | SMARCC2 | -0.96 | 4.9e-03 |
  | SMARCC2 | -0.96 | 4.9e-03 |
  | PLEC | -0.92 | 2.2e-02 |
  | MAD2L1 | -0.71 | 8.6e-04 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_INTERFERON_ALPHA_RESPONSE | 0.023 | +0.014 |
  | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITI | 0.039 | +0.011 |
  | HALLMARK_MYC_TARGETS_V2 | 0.031 | +0.009 |
  | HALLMARK_MTORC1_SIGNALING | 0.080 | +0.009 |
  | HALLMARK_SPERMATOGENESIS | 0.015 | +0.008 |
  | HALLMARK_INTERFERON_GAMMA_RESPONSE | 0.014 | +0.008 |
  | HALLMARK_UV_RESPONSE_UP | 0.021 | +0.007 |
  | HALLMARK_PEROXISOME | 0.034 | +0.007 |
  | HALLMARK_UV_RESPONSE_DN | 0.016 | +0.007 |
  | HALLMARK_PI3K_AKT_MTOR_SIGNALING | 0.047 | +0.006 |
  | HALLMARK_FATTY_ACID_METABOLISM | 0.048 | +0.004 |
  | HALLMARK_GLYCOLYSIS | 0.069 | +0.004 |
  | HALLMARK_XENOBIOTIC_METABOLISM | 0.009 | +0.004 |
  | HALLMARK_MITOTIC_SPINDLE | 0.038 | +0.003 |
  | HALLMARK_HYPOXIA | 0.048 | +0.003 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 2.43 | 13.80 | -0.58 | 0.4498 | 0.0051 |
  | TUBB3 | 1.00 | 0.96 | 1.72 | 15.88 | 0.23 | 0.7450 | 0.3776 |
  | DCX | 0.23 | 0.67 | -2.76 | 13.05 | -0.34 | 0.0000 | 0.7599 |
  | STMN1 | 1.00 | 1.00 | -4.27 | 15.80 | -0.20 | 1.0000 | 0.4623 |
  | STMN2 | 0.55 | 0.52 | 0.16 | 14.61 | -0.71 | 1.0000 | 0.0779 |
  | GAP43 | 0.00 | 0.08 | -2.92 | — | — | 0.2574 | — |
  | SYN1 | 0.00 | 0.02 | -0.81 | — | — | 1.0000 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.00 | 0.03 | -1.39 | — | — | 0.9531 | — |
  | NEFL | 0.00 | 0.03 | -1.51 | — | — | 0.9531 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.75 | 0.74 | 0.06 | 12.58 | -0.47 | 1.0000 | 0.0332 |
  | UCHL1 | 1.00 | 0.92 | 2.87 | 15.30 | -0.08 | 0.2560 | 0.3945 |
  | SOX2 | 0.30 | 0.10 | 1.93 | 11.85 | -0.54 | 0.0074 | 0.2252 |
  | NES | 0.97 | 0.73 | 3.31 | 15.38 | 0.82 | 0.0015 | 0.0001 |
  | VIM | 0.90 | 0.51 | 2.96 | 18.50 | 1.24 | 0.0000 | 0.0000 |
  | PAX6 | 0.05 | 0.05 | 0.28 | — | — | 1.0000 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 1.00 | 0.95 | 2.11 | 18.53 | -0.70 | 0.5766 | 0.1149 |
  | HOPX | 0.35 | 0.30 | 0.37 | 13.15 | -0.30 | 0.8393 | 0.8886 |
  | EOMES | 0.00 | 0.09 | -2.96 | — | — | 0.2628 | — |
  | TBR1 | 0.00 | 0.11 | -3.32 | — | — | 0.0936 | — |
  | SATB2 | 0.00 | 0.02 | -0.51 | — | — | 1.0000 | — |
  | BCL11B | 0.00 | 0.02 | -0.89 | — | — | 1.0000 | — |
  | CUX1 | 0.20 | 0.02 | 3.60 | 12.91 | -0.50 | 0.0001 | 0.4727 |
  | CUX2 | 0.00 | 0.07 | -2.55 | — | — | 0.4528 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | -0.14 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.33 | 0.77 | -2.78 | 14.69 | -0.50 | 0.0000 | 0.4345 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.07 | 0.11 | -0.47 | — | — | 0.9531 | — |
  | SLC1A3 | 0.03 | 0.02 | 0.57 | — | — | 1.0000 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | 0.55 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.12 | 0.14 | -0.05 | 11.25 | -0.87 | 1.0000 | 0.1956 |
  | AIF1 | 0.00 | 0.06 | -2.46 | — | — | 0.4510 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.25 | 0.10 | 1.60 | 12.77 | -1.13 | 0.0512 | 0.2604 |
  | COL1A2 | 0.45 | 0.16 | 2.07 | 14.21 | -0.16 | 0.0005 | 0.9858 |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.03 | -1.57 | — | — | 0.9531 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.45 | 0.53 | -0.45 | 13.27 | -0.80 | 0.7024 | 0.0180 |
  | TOP2A | 0.20 | 0.62 | -2.62 | 14.23 | -0.08 | 0.0000 | 0.8548 |
  | HMGB2 | 0.23 | 0.85 | -4.27 | 14.13 | 0.20 | 0.0000 | 0.5512 |
  | CDK1 | 0.23 | 0.73 | -3.16 | 13.85 | 0.14 | 0.0000 | 0.2296 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C21

- **n_cells**: 28
- **median_n_detected**: 999.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | PEA15 | 0.65 | 2/3 |
  | PSAT1 | 0.60 | 3/3 |
  | MAPK1 | 0.47 | 2/3 |
  | SRI | 0.46 | 2/3 |
  | MAT2A | 0.46 | 3/3 |
  | ALDH9A1 | 0.46 | 2/3 |
  | PPA1 | 0.45 | 3/3 |
  | TBCA | 0.44 | 3/3 |
  | JPT2 | 0.44 | 3/3 |
  | PFDN1 | 0.43 | 3/3 |
  | FABP7 | 0.43 | 2/3 |
  | GRHPR | 0.43 | 3/3 |
  | HMGB1 | 0.43 | 3/3 |
  | SEPTIN11 | 0.42 | 3/3 |
  | EZR | 0.41 | 2/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | JPT2 | 1.00 | 0.58 | 4.2e-05 |
  | NEDD4L | 1.00 | 0.59 | 7.7e-05 |
  | GGCT | 1.00 | 0.66 | 5.0e-04 |
  | SF3B1 | 1.00 | 0.79 | 2.3e-02 |
  | EDF1 | 1.00 | 0.70 | 1.6e-03 |
  | TOP2A | 1.00 | 0.58 | 4.9e-05 |
  | CHD4 | 1.00 | 0.65 | 5.0e-04 |
  | LANCL1 | 1.00 | 0.77 | 1.5e-02 |
  | SUPT16H | 1.00 | 0.79 | 2.3e-02 |
  | VDAC3 | 1.00 | 0.76 | 9.3e-03 |
  | PFDN1 | 1.00 | 0.57 | 4.2e-05 |
  | RUVBL2 | 1.00 | 0.80 | 3.4e-02 |
  | PSAT1 | 1.00 | 0.80 | 2.3e-02 |
  | RCN1 | 1.00 | 0.79 | 2.4e-02 |
  | DBN1 | 1.00 | 0.74 | 5.9e-03 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | CAPZA2 | 0.07 | 0.36 | 9.3e-03 |
  | OPA1 | 0.00 | 0.24 | 9.3e-03 |
  | C11orf96 | 0.04 | 0.25 | 3.8e-02 |
  | SNRPD1 | 0.64 | 0.86 | 3.1e-02 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | EIF3D | +1.98 | 1.7e-02 |
  | PEA15 | +1.33 | 8.4e-07 |
  | BCL7C | +1.27 | 4.4e-03 |
  | PSAT1 | +1.18 | 1.5e-06 |
  | RAB35 | +1.12 | 2.0e-04 |
  | TAGLN2 | +1.06 | 1.0e-04 |
  | VIM | +1.06 | 4.0e-03 |
  | KHDRBS1 | +1.05 | 5.2e-06 |
  | STMN1 | +1.01 | 5.0e-05 |
  | EIF1AX | +1.00 | 1.2e-04 |
  | RPL22 | +1.00 | 3.9e-05 |
  | MAT2A | +0.97 | 1.9e-09 |
  | HMGB1 | +0.97 | 8.5e-04 |
  | SOD1 | +0.96 | 9.1e-04 |
  | FABP7 | +0.93 | 5.1e-07 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | CD109 | -2.66 | 4.8e-02 |
  | PFKM | -2.26 | 1.8e-02 |
  | CTSZ | -2.10 | 5.1e-03 |
  | PDXP | -1.82 | 2.9e-03 |
  | ATAT1 | -1.76 | 6.5e-03 |
  | CSDE1 | -1.66 | 1.0e-02 |
  | CUL3 | -1.58 | 2.7e-03 |
  | PDS5B | -1.58 | 4.4e-03 |
  | KIF5C | -1.56 | 9.9e-03 |
  | MECP2 | -1.51 | 1.5e-02 |
  | TAGLN3 | -1.45 | 4.9e-04 |
  | NDUFB4 | -1.42 | 5.8e-03 |
  | HSDL1 | -1.36 | 1.7e-03 |
  | NOVA2 | -1.33 | 4.7e-04 |
  | GOT1 | -1.32 | 3.4e-02 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | KRT16 | +2.75 | 2.1e-02 |
  | NME2P1 | +2.57 | 9.3e-04 |
  | EEF1A2 | +2.33 | 3.6e-02 |
  | KRT86 | +2.20 | 2.1e-02 |
  | PSMC4 | +1.94 | 1.5e-03 |
  | ANXA1 | +1.91 | 3.6e-02 |
  | PEA15 | +1.82 | 5.2e-21 |
  | PITPNB | +1.82 | 8.4e-04 |
  | FHL1 | +1.68 | 4.4e-04 |
  | PSAT1 | +1.66 | 1.0e-17 |
  | SOD1 | +1.63 | 5.9e-08 |
  | ARG1 | +1.63 | 2.4e-02 |
  | VAPA | +1.63 | 1.1e-06 |
  | PSMF1 | +1.59 | 9.6e-03 |
  | DYNLL1 | +1.56 | 4.7e-11 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | CNP | -2.70 | 1.2e-02 |
  | SNCB | -2.62 | 1.5e-02 |
  | PDLIM1 | -2.19 | 1.1e-03 |
  | RAB21 | -2.01 | 1.9e-02 |
  | VCL | -1.89 | 1.3e-03 |
  | ZNF263 | -1.56 | 2.5e-03 |
  | RNASEH2C | -1.35 | 6.8e-04 |
  | FAM43B | -1.31 | 8.8e-03 |
  | CD207 | -1.17 | 6.4e-03 |
  | EXOSC8 | -1.07 | 3.9e-02 |
  | FKBP9 | -1.01 | 3.7e-03 |
  | S100A8 | -0.87 | 4.3e-02 |
  | SMARCC2 | -0.81 | 1.8e-02 |
  | SMARCC2 | -0.81 | 1.8e-02 |
  | WBP2 | -0.76 | 3.5e-02 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_MYC_TARGETS_V1 | 0.190 | +0.012 |
  | HALLMARK_UNFOLDED_PROTEIN_RESPONSE | 0.059 | +0.010 |
  | HALLMARK_G2M_CHECKPOINT | 0.096 | +0.008 |
  | HALLMARK_APICAL_SURFACE | 0.022 | +0.008 |
  | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | 0.049 | +0.005 |
  | HALLMARK_E2F_TARGETS | 0.083 | +0.005 |
  | HALLMARK_MTORC1_SIGNALING | 0.076 | +0.004 |
  | HALLMARK_MYC_TARGETS_V2 | 0.026 | +0.004 |
  | HALLMARK_ESTROGEN_RESPONSE_LATE | 0.036 | +0.003 |
  | HALLMARK_APOPTOSIS | 0.017 | +0.003 |
  | HALLMARK_PI3K_AKT_MTOR_SIGNALING | 0.044 | +0.002 |
  | HALLMARK_UV_RESPONSE_DN | 0.012 | +0.002 |
  | HALLMARK_GLYCOLYSIS | 0.067 | +0.001 |
  | HALLMARK_HEDGEHOG_SIGNALING | 0.092 | +0.001 |
  | HALLMARK_SPERMATOGENESIS | 0.009 | +0.001 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 1.89 | 14.92 | 0.62 | 0.8039 | 0.0094 |
  | TUBB3 | 1.00 | 0.96 | 1.19 | 16.03 | 0.39 | 1.0000 | 0.0331 |
  | DCX | 0.89 | 0.64 | 2.03 | 13.73 | 0.36 | 0.0310 | 0.3999 |
  | STMN1 | 1.00 | 1.00 | -4.80 | 16.93 | 1.01 | 1.0000 | 0.0000 |
  | STMN2 | 0.50 | 0.52 | -0.13 | 15.56 | 0.30 | 1.0000 | 0.5064 |
  | GAP43 | 0.00 | 0.08 | -2.39 | — | — | 0.4524 | — |
  | SYN1 | 0.00 | 0.02 | -0.28 | — | — | 1.0000 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.00 | 0.03 | -0.86 | — | — | 1.0000 | — |
  | NEFL | 0.00 | 0.03 | -0.98 | — | — | 1.0000 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.96 | 0.73 | 2.78 | 13.50 | 0.50 | 0.0265 | 0.0909 |
  | UCHL1 | 1.00 | 0.92 | 2.34 | 16.04 | 0.70 | 0.4524 | 0.0006 |
  | SOX2 | 0.21 | 0.11 | 1.24 | 12.81 | 0.53 | 0.3661 | 0.0174 |
  | NES | 0.96 | 0.73 | 2.76 | 14.77 | 0.19 | 0.0265 | 0.3159 |
  | VIM | 0.50 | 0.53 | -0.17 | 18.00 | 0.63 | 1.0000 | 0.0706 |
  | PAX6 | 0.18 | 0.05 | 2.28 | 12.18 | 0.28 | 0.0573 | 0.3794 |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 1.00 | 0.95 | 1.58 | 20.12 | 0.93 | 1.0000 | 0.0000 |
  | HOPX | 0.61 | 0.29 | 1.91 | 12.88 | -0.58 | 0.0070 | 0.1920 |
  | EOMES | 0.32 | 0.07 | 2.58 | 12.60 | 0.39 | 0.0033 | 0.8025 |
  | TBR1 | 0.04 | 0.11 | -1.13 | — | — | 0.7595 | — |
  | SATB2 | 0.00 | 0.02 | 0.02 | 17.38 | 0.70 | 1.0000 | 0.1265 |
  | BCL11B | 0.04 | 0.02 | 1.36 | — | — | 0.8597 | — |
  | CUX1 | 0.00 | 0.03 | -0.86 | — | — | 1.0000 | — |
  | CUX2 | 0.00 | 0.07 | -2.02 | — | — | 0.6105 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | 0.38 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.86 | 0.75 | 0.89 | 15.64 | 0.47 | 0.6279 | 0.1011 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.00 | 0.12 | -2.91 | — | — | 0.2296 | — |
  | SLC1A3 | 0.07 | 0.02 | 1.98 | — | — | 0.4430 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | 1.08 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.14 | 0.14 | 0.21 | — | — | 1.0000 | — |
  | AIF1 | 0.00 | 0.06 | -1.93 | — | — | 0.8044 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.18 | 0.11 | 0.97 | 13.68 | -0.08 | 0.5642 | 0.8158 |
  | COL1A2 | 0.21 | 0.18 | 0.43 | 13.91 | -0.52 | 1.0000 | 0.0909 |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.03 | -1.04 | — | — | 1.0000 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.89 | 0.51 | 2.80 | 14.05 | -0.02 | 0.0013 | 0.8948 |
  | TOP2A | 1.00 | 0.58 | 5.35 | 14.49 | 0.19 | 0.0000 | 0.3140 |
  | HMGB2 | 1.00 | 0.82 | 3.68 | 14.81 | 0.90 | 0.0517 | 0.0001 |
  | CDK1 | 0.86 | 0.70 | 1.23 | 14.33 | 0.65 | 0.3060 | 0.0001 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C3

- **n_cells**: 34
- **median_n_detected**: 1005.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | PEA15 | 0.66 | 2/3 |
  | PSAT1 | 0.61 | 3/3 |
  | HSPA5 | 0.56 | 2/3 |
  | PDIA4 | 0.53 | 3/3 |
  | FABP7 | 0.50 | 2/3 |
  | HSP90B1 | 0.48 | 2/3 |
  | P4HB | 0.45 | 2/3 |
  | LSS | 0.45 | 3/3 |
  | GLUD1 | 0.42 | 3/3 |
  | CALU | 0.42 | 3/3 |
  | MVD | 0.40 | 3/3 |
  | HMGCS1 | 0.40 | 3/3 |
  | PDIA6 | 0.40 | 3/3 |
  | TMED9 | 0.39 | 3/3 |
  | CALR | 0.39 | 2/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | IPO7 | 1.00 | 0.54 | 1.2e-07 |
  | ATP5PD | 1.00 | 0.77 | 4.2e-03 |
  | CYB5R3 | 1.00 | 0.83 | 2.6e-02 |
  | DYNLRB1 | 1.00 | 0.65 | 2.3e-05 |
  | DARS1 | 1.00 | 0.84 | 3.9e-02 |
  | TARDBP | 1.00 | 0.72 | 5.1e-04 |
  | CTSD | 1.00 | 0.81 | 1.1e-02 |
  | RPS6 | 1.00 | 0.83 | 2.6e-02 |
  | GDI1 | 1.00 | 0.83 | 2.6e-02 |
  | MARCKS | 1.00 | 0.79 | 6.9e-03 |
  | PSMB2 | 1.00 | 0.74 | 9.2e-04 |
  | ALDH7A1 | 1.00 | 0.72 | 5.1e-04 |
  | ATP5PO | 1.00 | 0.83 | 2.6e-02 |
  | LSS | 1.00 | 0.37 | 2.2e-12 |
  | FASN | 1.00 | 0.83 | 2.6e-02 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | COTL1 | 0.24 | 0.56 | 2.0e-03 |
  | RBP1 | 0.12 | 0.39 | 8.7e-03 |
  | PTMA | 0.00 | 0.25 | 1.5e-03 |
  | ELAVL2 | 0.15 | 0.39 | 2.6e-02 |
  | CORO1A | 0.03 | 0.21 | 4.7e-02 |
  | ILKAP | 0.03 | 0.20 | 4.8e-02 |
  | AK2 | 0.00 | 0.17 | 2.6e-02 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | PEA15 | +2.27 | 2.0e-15 |
  | GLUD1 | +1.85 | 2.0e-13 |
  | PSAT1 | +1.72 | 1.6e-16 |
  | CALR | +1.54 | 3.7e-13 |
  | HSPA5 | +1.52 | 1.7e-15 |
  | HMGCS1 | +1.49 | 6.9e-11 |
  | FABP7 | +1.49 | 8.8e-17 |
  | TUBB2A | +1.49 | 2.0e-06 |
  | PDIA4 | +1.47 | 9.1e-17 |
  | SEPTIN11 | +1.45 | 3.7e-13 |
  | TUBB2A | +1.44 | 2.6e-06 |
  | CALU | +1.33 | 5.7e-16 |
  | PLGRKT | +1.33 | 2.7e-03 |
  | FASN | +1.32 | 1.1e-11 |
  | P4HB | +1.31 | 2.9e-15 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | SLC9A3R1 | -3.30 | 3.8e-03 |
  | S100B | -2.85 | 9.8e-04 |
  | SAMHD1 | -2.00 | 1.2e-03 |
  | CTSZ | -1.76 | 4.2e-02 |
  | CUL3 | -1.63 | 4.7e-04 |
  | ELAVL3 | -1.40 | 8.7e-03 |
  | PRPF3 | -1.34 | 1.2e-02 |
  | SRSF5 | -1.29 | 8.0e-03 |
  | MIEN1 | -1.14 | 4.4e-02 |
  | DSC1 | -1.13 | 1.6e-02 |
  | YLPM1 | -1.04 | 2.4e-03 |
  | POU3F2 | -1.03 | 5.3e-04 |
  | ATP1B3 | -1.02 | 4.7e-03 |
  | SMARCE1 | -1.00 | 1.5e-02 |
  | NOVA2 | -0.97 | 7.0e-03 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | PEA15 | +2.67 | 6.3e-29 |
  | PSAT1 | +2.44 | 7.2e-32 |
  | EEF1A2 | +2.43 | 2.6e-03 |
  | HMGCS1 | +2.17 | 2.8e-21 |
  | RNMT | +2.17 | 1.0e-02 |
  | GLUD1 | +2.03 | 5.8e-22 |
  | HSPA5 | +2.01 | 2.3e-39 |
  | TMED9 | +1.95 | 7.5e-25 |
  | QKI | +1.93 | 6.1e-12 |
  | FABP7 | +1.93 | 2.9e-32 |
  | NDRG2 | +1.90 | 1.8e-21 |
  | PDIA4 | +1.88 | 8.0e-37 |
  | PTPRZ1 | +1.87 | 3.7e-07 |
  | PDIA6 | +1.85 | 1.5e-30 |
  | PSAP | +1.82 | 4.9e-05 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | AP2A2 | -5.11 | 1.7e-05 |
  | TF | -2.36 | 4.0e-07 |
  | VCL | -2.03 | 3.1e-04 |
  | RAB21 | -1.94 | 3.6e-02 |
  | P4HA1 | -1.75 | 2.6e-03 |
  | CNP | -1.74 | 6.5e-03 |
  | CA2 | -1.65 | 2.9e-02 |
  | STAU1 | -1.63 | 2.9e-02 |
  | KRT86 | -1.49 | 5.9e-03 |
  | RBM17 | -1.42 | 2.5e-02 |
  | SLC7A1 | -1.33 | 1.5e-02 |
  | FKBP1A | -1.31 | 1.8e-02 |
  | PGA3 | -1.27 | 7.6e-03 |
  | SERPINH1 | -1.22 | 5.0e-02 |
  | SMARCC2 | -1.20 | 5.5e-03 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_UNFOLDED_PROTEIN_RESPONSE | 0.089 | +0.041 |
  | HALLMARK_APICAL_SURFACE | 0.035 | +0.021 |
  | HALLMARK_MTORC1_SIGNALING | 0.091 | +0.020 |
  | HALLMARK_HEDGEHOG_SIGNALING | 0.104 | +0.014 |
  | HALLMARK_PI3K_AKT_MTOR_SIGNALING | 0.054 | +0.013 |
  | HALLMARK_APOPTOSIS | 0.025 | +0.010 |
  | HALLMARK_MYC_TARGETS_V2 | 0.032 | +0.010 |
  | HALLMARK_ESTROGEN_RESPONSE_LATE | 0.040 | +0.006 |
  | HALLMARK_BILE_ACID_METABOLISM | 0.012 | +0.006 |
  | HALLMARK_UV_RESPONSE_DN | 0.015 | +0.006 |
  | HALLMARK_OXIDATIVE_PHOSPHORYLATION | 0.031 | +0.005 |
  | HALLMARK_COMPLEMENT | 0.029 | +0.004 |
  | HALLMARK_ESTROGEN_RESPONSE_EARLY | 0.029 | +0.003 |
  | HALLMARK_CHOLESTEROL_HOMEOSTASIS | 0.035 | +0.003 |
  | HALLMARK_MYC_TARGETS_V1 | 0.180 | +0.003 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 2.18 | 14.52 | 0.20 | 0.6036 | 0.4804 |
  | TUBB3 | 1.00 | 0.96 | 1.48 | 16.17 | 0.54 | 0.9741 | 0.0236 |
  | DCX | 1.00 | 0.63 | 5.31 | 12.99 | -0.43 | 0.0000 | 0.0404 |
  | STMN1 | 1.00 | 1.00 | -4.51 | 15.87 | -0.13 | 1.0000 | 0.7378 |
  | STMN2 | 0.68 | 0.52 | 0.94 | 14.70 | -0.62 | 0.2798 | 0.0621 |
  | GAP43 | 0.00 | 0.08 | -2.67 | — | — | 0.3359 | — |
  | SYN1 | 0.00 | 0.02 | -0.57 | — | — | 1.0000 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.00 | 0.03 | -1.15 | — | — | 0.9741 | — |
  | NEFL | 0.03 | 0.03 | 0.42 | — | — | 1.0000 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.94 | 0.73 | 2.29 | 13.04 | 0.01 | 0.0312 | 0.6370 |
  | UCHL1 | 1.00 | 0.92 | 2.62 | 16.26 | 0.93 | 0.3358 | 0.0000 |
  | SOX2 | 0.06 | 0.11 | -0.76 | — | — | 0.7863 | — |
  | NES | 0.88 | 0.73 | 1.31 | 15.00 | 0.42 | 0.2564 | 0.1304 |
  | VIM | 0.68 | 0.52 | 0.90 | 18.39 | 1.03 | 0.2873 | 0.0030 |
  | PAX6 | 0.09 | 0.05 | 1.11 | — | — | 0.6036 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 0.97 | 0.95 | 0.20 | 20.67 | 1.49 | 1.0000 | 0.0000 |
  | HOPX | 0.94 | 0.27 | 5.12 | 14.33 | 1.05 | 0.0000 | 0.0000 |
  | EOMES | 0.00 | 0.09 | -2.72 | — | — | 0.3373 | — |
  | TBR1 | 0.00 | 0.11 | -3.07 | — | — | 0.1772 | — |
  | SATB2 | 0.00 | 0.02 | -0.27 | 17.76 | 1.06 | 1.0000 | 0.0619 |
  | BCL11B | 0.00 | 0.02 | -0.65 | — | — | 1.0000 | — |
  | CUX1 | 0.00 | 0.03 | -1.15 | — | — | 0.9741 | — |
  | CUX2 | 0.00 | 0.07 | -2.30 | — | — | 0.4574 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | 0.10 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.88 | 0.74 | 1.22 | 15.47 | 0.32 | 0.2581 | 0.1202 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.24 | 0.11 | 1.41 | 12.06 | -2.85 | 0.1881 | 0.0010 |
  | SLC1A3 | 0.24 | 0.02 | 4.31 | 12.14 | 0.24 | 0.0000 | 0.7499 |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | 0.80 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.32 | 0.13 | 1.72 | 12.05 | -0.11 | 0.0260 | 0.8395 |
  | AIF1 | 0.00 | 0.06 | -2.21 | — | — | 0.6036 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.06 | 0.11 | -0.70 | — | — | 0.9615 | — |
  | COL1A2 | 0.21 | 0.18 | 0.35 | 14.80 | 0.50 | 0.9826 | 0.3484 |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.03 | -1.32 | — | — | 0.9741 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.68 | 0.52 | 0.92 | 13.76 | -0.31 | 0.2823 | 0.3728 |
  | TOP2A | 0.74 | 0.59 | 0.89 | 14.12 | -0.23 | 0.3452 | 0.6783 |
  | HMGB2 | 0.88 | 0.82 | 0.57 | 13.43 | -0.54 | 0.8838 | 0.0165 |
  | CDK1 | 0.94 | 0.69 | 2.52 | 13.90 | 0.22 | 0.0087 | 0.3804 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C4

- **n_cells**: 26
- **median_n_detected**: 862.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | MYL6 | 0.61 | 2/3 |
  | LGALS1 | 0.60 | 3/3 |
  | MYL12A | 0.54 | 3/3 |
  | TPM4 | 0.46 | 2/3 |
  | CLIC1 | 0.43 | 3/3 |
  | LMNA | 0.39 | 3/3 |
  | FABP7 | 0.39 | 3/3 |
  | ACTN4 | 0.38 | 3/3 |
  | RAP1B | 0.37 | 3/3 |
  | FLNA | 0.36 | 3/3 |
  | PLOD1 | 0.33 | 1/3 |
  | THY1 | 0.33 | 2/3 |
  | VCL | 0.33 | 3/3 |
  | PFN1 | 0.33 | 2/3 |
  | TUBB2B | 0.33 | 2/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | MYO1B | 1.00 | 0.04 | 3.5e-31 |
  | LMNA | 1.00 | 0.18 | 9.3e-18 |
  | GSN | 1.00 | 0.37 | 2.8e-10 |
  | CLIC4 | 1.00 | 0.45 | 2.4e-08 |
  | ANXA1 | 1.00 | 0.07 | 2.3e-26 |
  | THY1 | 1.00 | 0.01 | 1.7e-38 |
  | PDLIM1 | 1.00 | 0.07 | 3.7e-27 |
  | TLN1 | 1.00 | 0.33 | 1.9e-11 |
  | S100A11 | 1.00 | 0.07 | 3.2e-26 |
  | SAE1 | 1.00 | 0.80 | 2.6e-02 |
  | ITGB1 | 1.00 | 0.11 | 7.6e-23 |
  | ISYNA1 | 1.00 | 0.72 | 3.0e-03 |
  | ALDH2 | 1.00 | 0.54 | 2.5e-06 |
  | VCL | 1.00 | 0.28 | 5.2e-13 |
  | MYOF | 1.00 | 0.04 | 1.0e-30 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | DPYSL5 | 0.04 | 0.92 | 3.0e-23 |
  | UCHL1 | 0.08 | 0.95 | 5.9e-25 |
  | MAP1B | 0.12 | 0.95 | 1.0e-22 |
  | DPYSL3 | 0.08 | 0.90 | 5.0e-19 |
  | NFIB | 0.08 | 0.85 | 7.9e-16 |
  | PSAT1 | 0.08 | 0.83 | 1.8e-14 |
  | JPT1 | 0.12 | 0.84 | 1.1e-13 |
  | DDAH1 | 0.23 | 0.93 | 1.1e-15 |
  | MYEF2 | 0.04 | 0.73 | 6.2e-12 |
  | CKB | 0.31 | 1.00 | 4.7e-26 |
  | HMGCS1 | 0.08 | 0.76 | 1.5e-11 |
  | ENO2 | 0.08 | 0.76 | 2.3e-11 |
  | ALDH6A1 | 0.12 | 0.77 | 1.3e-10 |
  | MSI1 | 0.04 | 0.70 | 1.2e-10 |
  | TMSB15A | 0.15 | 0.81 | 7.5e-11 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | LMNA | +3.99 | 2.6e-13 |
  | LGALS1 | +3.98 | 3.5e-15 |
  | ITGB6 | +3.72 | 1.6e-07 |
  | THY1 | +3.62 | 6.9e-04 |
  | VCL | +3.58 | 5.8e-14 |
  | TPM4 | +3.45 | 3.5e-10 |
  | CALD1 | +3.43 | 5.9e-07 |
  | TPM4 | +3.39 | 6.9e-07 |
  | MYL12A | +3.32 | 3.5e-15 |
  | MYL6 | +3.25 | 3.5e-15 |
  | RAP1B | +2.99 | 2.2e-12 |
  | MYH9 | +2.91 | 2.6e-14 |
  | ITGB1 | +2.70 | 3.8e-11 |
  | SLC9A3R1 | +2.65 | 6.4e-04 |
  | TMEM109 | +2.64 | 3.0e-05 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | FABP7 | -5.85 | 1.9e-08 |
  | CKB | -4.39 | 5.6e-06 |
  | TUBB2B | -3.22 | 7.3e-13 |
  | TUBB2A | -3.16 | 9.8e-11 |
  | FSCN1 | -2.84 | 3.2e-08 |
  | TUBB2A | -2.84 | 3.4e-04 |
  | FABP5 | -2.83 | 3.2e-12 |
  | PCNA | -2.61 | 2.0e-06 |
  | PSIP1 | -2.48 | 1.5e-04 |
  | SCRN1 | -1.93 | 5.9e-06 |
  | PAFAH1B3 | -1.90 | 1.8e-05 |
  | STMN2 | -1.89 | 2.2e-04 |
  | PFN2 | -1.86 | 2.1e-08 |
  | TUBB3 | -1.82 | 2.8e-11 |
  | TOP2B | -1.80 | 1.0e-03 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | PLAAT3 | +4.02 | 3.6e-02 |
  | LGALS1 | +3.91 | 1.2e-34 |
  | ITGB6 | +3.66 | 9.9e-19 |
  | LMNA | +3.64 | 4.1e-04 |
  | RAP1B | +3.60 | 1.4e-15 |
  | TPM4 | +3.52 | 7.9e-41 |
  | TPM4 | +3.52 | 7.9e-41 |
  | MYL12A | +3.49 | 6.3e-39 |
  | MYL6 | +3.38 | 2.4e-51 |
  | MYH9 | +3.34 | 6.7e-12 |
  | TLN1 | +3.09 | 2.1e-14 |
  | AK1 | +2.79 | 7.5e-19 |
  | FLNA | +2.67 | 3.5e-34 |
  | ACTN4 | +2.60 | 5.5e-37 |
  | CLIC1 | +2.60 | 6.1e-45 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | KRT71 | -6.86 | 3.3e-05 |
  | FABP7 | -4.33 | 3.6e-14 |
  | AP2A2 | -3.70 | 4.8e-02 |
  | CKB | -3.33 | 1.9e-07 |
  | FSCN1 | -2.90 | 1.1e-17 |
  | TUBB2B | -2.67 | 1.4e-20 |
  | MYLK | -2.54 | 1.7e-05 |
  | PCNA | -2.53 | 1.6e-19 |
  | FABP5 | -2.43 | 1.0e-22 |
  | STMN1 | -2.41 | 4.4e-19 |
  | AHNAK | -2.25 | 1.3e-04 |
  | PSMD9 | -2.23 | 7.0e-04 |
  | STMN2 | -2.08 | 2.7e-09 |
  | JPT1 | -2.04 | 4.8e-02 |
  | TUBB2A | -2.00 | 9.8e-03 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITI | 0.092 | +0.066 |
  | HALLMARK_APICAL_JUNCTION | 0.078 | +0.043 |
  | HALLMARK_HYPOXIA | 0.086 | +0.041 |
  | HALLMARK_HEDGEHOG_SIGNALING | 0.130 | +0.041 |
  | HALLMARK_OXIDATIVE_PHOSPHORYLATION | 0.053 | +0.028 |
  | HALLMARK_FATTY_ACID_METABOLISM | 0.069 | +0.025 |
  | HALLMARK_PI3K_AKT_MTOR_SIGNALING | 0.065 | +0.024 |
  | HALLMARK_APICAL_SURFACE | 0.038 | +0.024 |
  | HALLMARK_MITOTIC_SPINDLE | 0.056 | +0.022 |
  | HALLMARK_COMPLEMENT | 0.046 | +0.022 |
  | HALLMARK_ALLOGRAFT_REJECTION | 0.056 | +0.017 |
  | HALLMARK_GLYCOLYSIS | 0.082 | +0.017 |
  | HALLMARK_MTORC1_SIGNALING | 0.086 | +0.015 |
  | HALLMARK_TGF_BETA_SIGNALING | 0.013 | +0.012 |
  | HALLMARK_COAGULATION | 0.014 | +0.009 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 0.88 | 0.94 | -1.29 | 12.78 | -1.60 | 0.5132 | 0.0000 |
  | TUBB3 | 0.96 | 0.96 | -0.61 | 13.88 | -1.82 | 1.0000 | 0.0000 |
  | DCX | 0.04 | 0.67 | -5.11 | — | — | 0.0000 | — |
  | STMN1 | 1.00 | 1.00 | -4.90 | 14.32 | -1.68 | 1.0000 | 0.0000 |
  | STMN2 | 0.35 | 0.53 | -1.05 | 13.41 | -1.89 | 0.2624 | 0.0002 |
  | GAP43 | 0.00 | 0.08 | -2.28 | — | — | 0.6067 | — |
  | SYN1 | 0.00 | 0.02 | -0.17 | — | — | 1.0000 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.00 | 0.03 | -0.75 | — | — | 1.0000 | — |
  | NEFL | 0.00 | 0.03 | -0.87 | — | — | 1.0000 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.08 | 0.76 | -4.94 | — | — | 0.0000 | — |
  | UCHL1 | 0.08 | 0.95 | -7.51 | — | — | 0.0000 | — |
  | SOX2 | 0.00 | 0.12 | -2.81 | — | — | 0.3338 | — |
  | NES | 0.50 | 0.75 | -1.56 | 14.07 | -0.53 | 0.0517 | 0.0022 |
  | VIM | 0.62 | 0.53 | 0.50 | 17.61 | 0.24 | 0.8094 | 0.4542 |
  | PAX6 | 0.08 | 0.05 | 0.96 | — | — | 0.7786 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 0.62 | 0.96 | -4.04 | 13.39 | -5.85 | 0.0000 | 0.0000 |
  | HOPX | 0.04 | 0.31 | -2.92 | — | — | 0.0105 | — |
  | EOMES | 0.00 | 0.09 | -2.32 | — | — | 0.4475 | — |
  | TBR1 | 0.00 | 0.11 | -2.68 | — | — | 0.3327 | — |
  | SATB2 | 0.00 | 0.02 | 0.13 | 16.87 | 0.17 | 1.0000 | 0.8879 |
  | BCL11B | 0.00 | 0.02 | -0.26 | — | — | 1.0000 | — |
  | CUX1 | 0.00 | 0.03 | -0.75 | — | — | 1.0000 | — |
  | CUX2 | 0.00 | 0.07 | -1.91 | — | — | 0.7804 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | 0.49 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.62 | 0.75 | -0.96 | 15.45 | 0.27 | 0.3519 | 0.3964 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.00 | 0.12 | -2.81 | — | — | 0.3338 | — |
  | SLC1A3 | 0.00 | 0.03 | -0.49 | — | — | 1.0000 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | 1.19 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.35 | 0.13 | 1.86 | 11.20 | -0.98 | 0.0271 | 0.0057 |
  | AIF1 | 0.00 | 0.06 | -1.82 | — | — | 0.7786 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.12 | 0.11 | 0.28 | — | — | 1.0000 | — |
  | COL1A2 | 0.08 | 0.18 | -1.11 | — | — | 0.6696 | — |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.92 | 0.00 | 11.11 | — | — | 0.0000 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.04 | 0.54 | -4.32 | — | — | 0.0000 | — |
  | TOP2A | 0.08 | 0.61 | -3.96 | — | — | 0.0000 | — |
  | HMGB2 | 0.62 | 0.83 | -1.63 | 13.00 | -0.96 | 0.0690 | 0.0002 |
  | CDK1 | 0.38 | 0.72 | -1.98 | 12.72 | -1.03 | 0.0046 | 0.0015 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C5

- **n_cells**: 29
- **median_n_detected**: 1009.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | DBI | 0.67 | 2/3 |
  | ANXA5 | 0.49 | 2/3 |
  | IDH1 | 0.47 | 3/3 |
  | HIP1 | 0.33 | 3/3 |
  | COL4A1 | 0.33 | 3/3 |
  | DCLK2 | 0.28 | 3/3 |
  | NME1 | 0.27 | 2/3 |
  | PDIA6 | 0.27 | 2/3 |
  | S100B | 0.27 | 3/3 |
  | ENO1 | 0.26 | 2/3 |
  | EEF2 | 0.26 | 2/3 |
  | MTHFD1 | 0.25 | 3/3 |
  | SRI | 0.25 | 2/3 |
  | PDIA4 | 0.25 | 3/3 |
  | P4HB | 0.25 | 2/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | IPO7 | 1.00 | 0.54 | 3.4e-06 |
  | MCM7 | 1.00 | 0.75 | 6.4e-03 |
  | MCM5 | 1.00 | 0.75 | 6.4e-03 |
  | MCM4 | 1.00 | 0.74 | 4.0e-03 |
  | SHMT2 | 1.00 | 0.76 | 6.6e-03 |
  | RANBP1 | 1.00 | 0.80 | 2.4e-02 |
  | PAFAH1B1 | 1.00 | 0.67 | 4.3e-04 |
  | NES | 1.00 | 0.73 | 2.6e-03 |
  | ATIC | 1.00 | 0.76 | 6.5e-03 |
  | MCM3 | 1.00 | 0.75 | 6.4e-03 |
  | SRM | 1.00 | 0.73 | 4.0e-03 |
  | GFAP | 1.00 | 0.74 | 4.0e-03 |
  | MTHFD1 | 1.00 | 0.81 | 3.9e-02 |
  | YBX1 | 1.00 | 0.79 | 1.6e-02 |
  | MCM2 | 1.00 | 0.77 | 1.0e-02 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | PTBP2 | 0.31 | 0.84 | 7.6e-08 |
  | UBE2N | 0.10 | 0.50 | 3.3e-04 |
  | EEF1A1 | 0.10 | 0.50 | 3.3e-04 |
  | UBE2V1 | 0.10 | 0.50 | 3.3e-04 |
  | RBMX | 0.10 | 0.50 | 3.3e-04 |
  | KHDRBS1 | 0.10 | 0.50 | 3.5e-04 |
  | RAB5A | 0.10 | 0.50 | 3.5e-04 |
  | CALM1 | 0.10 | 0.49 | 6.7e-04 |
  | KRT13 | 0.10 | 0.49 | 6.7e-04 |
  | MAGOH | 0.10 | 0.48 | 7.4e-04 |
  | NEDD4L | 0.24 | 0.62 | 1.3e-03 |
  | ELAVL3 | 0.10 | 0.47 | 1.4e-03 |
  | SRSF4 | 0.10 | 0.47 | 1.4e-03 |
  | H3-7 | 0.10 | 0.46 | 1.5e-03 |
  | TUBB2A | 0.10 | 0.45 | 1.6e-03 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | DBI | +1.78 | 3.9e-15 |
  | ANXA5 | +1.44 | 5.3e-13 |
  | LGALS1 | +1.42 | 1.8e-05 |
  | IDH1 | +1.38 | 2.6e-12 |
  | COTL1 | +1.34 | 7.6e-03 |
  | VIM | +1.34 | 3.0e-06 |
  | NME1 | +1.26 | 1.6e-11 |
  | ENO1 | +1.12 | 1.6e-11 |
  | ACTN4 | +1.11 | 3.3e-08 |
  | SNRNP40 | +1.08 | 1.4e-03 |
  | LGALS3 | +1.08 | 1.4e-06 |
  | EEF2 | +1.08 | 3.5e-12 |
  | NES | +1.07 | 7.4e-08 |
  | FSCN1 | +1.03 | 5.5e-06 |
  | DCLK2 | +1.03 | 2.5e-05 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | ALB | -2.39 | 1.8e-02 |
  | HMGA1 | -1.81 | 1.7e-03 |
  | PDXP | -1.77 | 2.8e-02 |
  | RAB35 | -1.60 | 3.1e-02 |
  | CEP170 | -1.51 | 6.3e-04 |
  | S100B | -1.39 | 1.0e-02 |
  | CASP14 | -1.23 | 3.0e-02 |
  | DAZAP1 | -1.19 | 8.2e-03 |
  | DYNLT1 | -1.19 | 2.5e-02 |
  | H1-0 | -1.18 | 8.3e-03 |
  | DKC1 | -1.16 | 4.7e-02 |
  | EIF4G2 | -1.15 | 7.5e-03 |
  | BUB3 | -1.14 | 8.2e-05 |
  | FKBP1A | -1.12 | 7.1e-03 |
  | ERH | -1.10 | 1.3e-04 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | DBI | +2.04 | 8.6e-31 |
  | ANXA5 | +1.91 | 4.1e-26 |
  | EML5 | +1.91 | 6.5e-03 |
  | S100B | +1.73 | 4.4e-02 |
  | IDH1 | +1.70 | 5.4e-29 |
  | PITPNB | +1.65 | 1.1e-02 |
  | QKI | +1.63 | 2.9e-09 |
  | SCARB2 | +1.56 | 6.3e-04 |
  | FHL1 | +1.54 | 3.0e-05 |
  | PDIA6 | +1.48 | 2.1e-22 |
  | LGALS3 | +1.47 | 1.0e-11 |
  | TUBB2A | +1.45 | 2.5e-03 |
  | TUBB2A | +1.45 | 2.5e-03 |
  | NME2P1 | +1.45 | 6.8e-03 |
  | ZNF148 | +1.45 | 4.4e-02 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | AP2A2 | -5.17 | 3.8e-03 |
  | NDRG1 | -2.08 | 4.4e-02 |
  | AHNAK | -1.95 | 2.1e-02 |
  | RAB21 | -1.77 | 4.3e-04 |
  | VCL | -1.73 | 2.4e-02 |
  | FKBP9 | -1.62 | 1.0e-02 |
  | CEP170 | -1.36 | 6.2e-05 |
  | TMPO | -1.19 | 9.5e-04 |
  | AIFM1 | -1.13 | 1.6e-02 |
  | FKBP1A | -1.10 | 1.9e-02 |
  | ZNF263 | -1.05 | 3.2e-02 |
  | CSRP2 | -1.03 | 5.6e-03 |
  | GRM3 | -0.98 | 3.2e-02 |
  | DKC1 | -0.91 | 3.2e-02 |
  | XRN2 | -0.85 | 2.1e-02 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_MTORC1_SIGNALING | 0.092 | +0.021 |
  | HALLMARK_GLYCOLYSIS | 0.085 | +0.021 |
  | HALLMARK_PI3K_AKT_MTOR_SIGNALING | 0.059 | +0.018 |
  | HALLMARK_HYPOXIA | 0.062 | +0.017 |
  | HALLMARK_UNFOLDED_PROTEIN_RESPONSE | 0.065 | +0.016 |
  | HALLMARK_CHOLESTEROL_HOMEOSTASIS | 0.046 | +0.014 |
  | HALLMARK_OXIDATIVE_PHOSPHORYLATION | 0.038 | +0.012 |
  | HALLMARK_ALLOGRAFT_REJECTION | 0.051 | +0.012 |
  | HALLMARK_COMPLEMENT | 0.036 | +0.011 |
  | HALLMARK_MYC_TARGETS_V1 | 0.187 | +0.010 |
  | HALLMARK_FATTY_ACID_METABOLISM | 0.053 | +0.010 |
  | HALLMARK_APICAL_SURFACE | 0.023 | +0.008 |
  | HALLMARK_INTERFERON_ALPHA_RESPONSE | 0.016 | +0.007 |
  | HALLMARK_INTERFERON_GAMMA_RESPONSE | 0.013 | +0.007 |
  | HALLMARK_ANDROGEN_RESPONSE | 0.024 | +0.006 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 1.95 | 15.03 | 0.73 | 0.8435 | 0.0011 |
  | TUBB3 | 1.00 | 0.96 | 1.24 | 15.59 | -0.09 | 1.0000 | 0.5292 |
  | DCX | 0.83 | 0.64 | 1.31 | 13.19 | -0.21 | 0.2131 | 0.5721 |
  | STMN1 | 1.00 | 1.00 | -4.74 | 16.38 | 0.45 | 1.0000 | 0.0382 |
  | STMN2 | 0.86 | 0.51 | 2.44 | 15.19 | -0.11 | 0.0029 | 0.7947 |
  | GAP43 | 0.34 | 0.07 | 2.81 | 12.77 | 0.05 | 0.0007 | 0.7128 |
  | SYN1 | 0.00 | 0.02 | -0.33 | — | — | 1.0000 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.00 | 0.03 | -0.92 | — | — | 1.0000 | — |
  | NEFL | 0.03 | 0.03 | 0.66 | — | — | 1.0000 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.79 | 0.73 | 0.39 | 13.02 | -0.01 | 1.0000 | 0.6664 |
  | UCHL1 | 1.00 | 0.92 | 2.39 | 14.62 | -0.81 | 0.5094 | 0.0002 |
  | SOX2 | 0.03 | 0.12 | -1.31 | — | — | 0.6545 | — |
  | NES | 1.00 | 0.73 | 4.46 | 15.62 | 1.07 | 0.0026 | 0.0000 |
  | VIM | 0.90 | 0.52 | 2.83 | 18.63 | 1.34 | 0.0007 | 0.0000 |
  | PAX6 | 0.00 | 0.05 | -1.71 | — | — | 0.8413 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 1.00 | 0.95 | 1.63 | 19.08 | -0.14 | 0.8413 | 0.8884 |
  | HOPX | 0.79 | 0.28 | 3.21 | 13.67 | 0.34 | 0.0000 | 0.2229 |
  | EOMES | 0.03 | 0.08 | -0.83 | — | — | 0.9643 | — |
  | TBR1 | 0.00 | 0.11 | -2.84 | — | — | 0.2683 | — |
  | SATB2 | 0.00 | 0.02 | -0.03 | 17.13 | 0.43 | 1.0000 | 0.4367 |
  | BCL11B | 0.00 | 0.02 | -0.42 | — | — | 1.0000 | — |
  | CUX1 | 0.00 | 0.03 | -0.92 | — | — | 1.0000 | — |
  | CUX2 | 0.00 | 0.07 | -2.07 | — | — | 0.6545 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | 0.33 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 1.00 | 0.74 | 4.37 | 16.06 | 0.94 | 0.0040 | 0.0000 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.79 | 0.09 | 5.23 | 13.87 | -1.39 | 0.0000 | 0.0105 |
  | SLC1A3 | 0.10 | 0.02 | 2.55 | — | — | 0.1528 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | 1.03 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.52 | 0.12 | 2.93 | 12.02 | -0.10 | 0.0000 | 0.6624 |
  | AIF1 | 0.00 | 0.06 | -1.98 | — | — | 0.8451 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.03 | 0.11 | -1.26 | — | — | 0.8142 | — |
  | COL1A2 | 0.07 | 0.18 | -1.29 | — | — | 0.4857 | — |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.03 | -1.09 | — | — | 1.0000 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.76 | 0.52 | 1.49 | 14.17 | 0.11 | 0.0782 | 0.9735 |
  | TOP2A | 0.79 | 0.59 | 1.33 | 14.20 | -0.11 | 0.1654 | 0.7007 |
  | HMGB2 | 0.93 | 0.82 | 1.29 | 13.74 | -0.21 | 0.4857 | 0.3392 |
  | CDK1 | 0.97 | 0.69 | 3.06 | 14.29 | 0.62 | 0.0074 | 0.0006 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C6

- **n_cells**: 26
- **median_n_detected**: 862.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | VIM | 0.78 | 3/3 |
  | RAB2A | 0.61 | 3/3 |
  | UBE2N | 0.61 | 3/3 |
  | RPS27A | 0.50 | 3/3 |
  | STMN2 | 0.50 | 3/3 |
  | UBE2V2 | 0.46 | 3/3 |
  | MAGOHB | 0.42 | 3/3 |
  | RPL22 | 0.37 | 3/3 |
  | DSTN | 0.34 | 3/3 |
  | DYNLL1 | 0.33 | 3/3 |
  | KHDRBS1 | 0.32 | 3/3 |
  | SYNE2 | 0.31 | 3/3 |
  | NEDD4L | 0.31 | 3/3 |
  | RPL36AL | 0.31 | 3/3 |
  | RBMX | 0.30 | 3/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | TUBB4A | 1.00 | 0.63 | 1.2e-03 |
  | ISYNA1 | 1.00 | 0.72 | 1.7e-02 |
  | BANF1 | 1.00 | 0.74 | 2.5e-02 |
  | TMSB15A | 1.00 | 0.78 | 4.7e-02 |
  | MAP4 | 1.00 | 0.75 | 2.5e-02 |
  | CDC42 | 1.00 | 0.77 | 4.7e-02 |
  | DYNLL1 | 1.00 | 0.60 | 6.3e-04 |
  | ALDH7A1 | 1.00 | 0.72 | 1.7e-02 |
  | FEN1 | 1.00 | 0.66 | 3.4e-03 |
  | ALDH6A1 | 1.00 | 0.74 | 2.5e-02 |
  | CACYBP | 1.00 | 0.77 | 4.7e-02 |
  | RTRAF | 1.00 | 0.77 | 3.3e-02 |
  | MAT2B | 1.00 | 0.73 | 1.7e-02 |
  | CPSF7 | 1.00 | 0.74 | 2.5e-02 |
  | SYNE2 | 1.00 | 0.59 | 5.2e-04 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | RPS27A | 0.19 | 0.61 | 2.0e-03 |
  | TPM4 | 0.19 | 0.61 | 3.3e-03 |
  | RAB2A | 0.19 | 0.57 | 9.5e-03 |
  | PPP2CA | 0.19 | 0.55 | 1.7e-02 |
  | GATAD2A | 0.19 | 0.54 | 1.7e-02 |
  | VIM | 0.19 | 0.54 | 1.7e-02 |
  | STMN2 | 0.19 | 0.53 | 1.9e-02 |
  | UBE2V2 | 0.19 | 0.52 | 2.7e-02 |
  | UBE2N | 0.19 | 0.52 | 2.7e-02 |
  | EEF1A1 | 0.19 | 0.52 | 2.7e-02 |
  | CALM1 | 0.19 | 0.52 | 2.7e-02 |
  | RBMX | 0.19 | 0.52 | 2.7e-02 |
  | KHDRBS1 | 0.19 | 0.52 | 2.8e-02 |
  | RAB5C | 0.19 | 0.52 | 2.8e-02 |
  | HSPA1A | 0.19 | 0.51 | 2.8e-02 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | STMN2 | +2.95 | 4.2e-04 |
  | ALB | +2.84 | 4.7e-02 |
  | RBMX | +2.28 | 2.4e-04 |
  | KHDRBS1 | +2.23 | 9.4e-05 |
  | TUBB4B | +2.12 | 2.8e-02 |
  | RAB2A | +2.00 | 7.4e-05 |
  | TAGLN3 | +1.99 | 8.3e-03 |
  | MAGOHB | +1.94 | 5.0e-04 |
  | UBE2N | +1.89 | 1.4e-04 |
  | RBP1 | +1.87 | 5.8e-03 |
  | SEPTIN3 | +1.86 | 4.4e-02 |
  | DAZAP1 | +1.85 | 1.1e-02 |
  | RPS27A | +1.76 | 1.4e-04 |
  | RBM39 | +1.76 | 4.9e-03 |
  | H3-4 | +1.74 | 6.5e-03 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | VIM | -4.65 | 3.8e-03 |
  | HNRNPAB | -1.47 | 4.1e-02 |
  | ARHGDIA | -1.45 | 2.2e-02 |
  | IST1 | -1.17 | 3.7e-03 |
  | ITGB1 | -1.06 | 4.4e-02 |
  | RRM2 | -0.93 | 4.5e-02 |
  | ERP29 | -0.68 | 3.8e-02 |
  | NUP153 | -0.58 | 7.3e-03 |
  | MAD2L1 | -0.41 | 4.9e-02 |
  | SHMT2 | -0.38 | 1.2e-02 |
  | MYBPC1 | -0.24 | 5.0e-02 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | NUCKS1 | +3.16 | 1.8e-02 |
  | LSM4 | +2.87 | 2.0e-02 |
  | H1-2 | +2.83 | 5.8e-04 |
  | BRD3 | +2.76 | 8.8e-03 |
  | RAB2A | +2.64 | 4.2e-22 |
  | TAGLN3 | +2.60 | 2.9e-09 |
  | STMN2 | +2.54 | 3.7e-12 |
  | OSBPL9 | +2.54 | 1.3e-04 |
  | MAPT | +2.50 | 2.5e-03 |
  | NFIX | +2.45 | 4.8e-03 |
  | WTAP | +2.45 | 4.1e-02 |
  | SEPTIN5 | +2.39 | 4.8e-04 |
  | GNA13 | +2.39 | 4.1e-05 |
  | PDHX | +2.31 | 1.4e-02 |
  | EEF1A2 | +2.31 | 5.9e-03 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | TERF2IP | -5.42 | 1.4e-04 |
  | VIM | -4.26 | 2.0e-15 |
  | VIM | -4.26 | 2.0e-15 |
  | AHNAK | -3.20 | 7.6e-04 |
  | CNP | -2.72 | 9.4e-03 |
  | GTF3C5 | -2.13 | 1.4e-03 |
  | RAB21 | -2.12 | 7.6e-05 |
  | VCL | -1.94 | 1.4e-03 |
  | SERPINH1 | -1.73 | 1.0e-02 |
  | CDSN | -1.40 | 3.3e-02 |
  | C11orf54 | -1.36 | 2.4e-03 |
  | FLRT3 | -1.25 | 1.6e-04 |
  | SLC7A1 | -1.08 | 7.0e-03 |
  | SMARCC2 | -1.06 | 1.8e-02 |
  | SMARCC2 | -1.06 | 1.8e-02 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_E2F_TARGETS | 0.092 | +0.014 |
  | HALLMARK_G2M_CHECKPOINT | 0.100 | +0.012 |
  | HALLMARK_MYC_TARGETS_V1 | 0.189 | +0.012 |
  | HALLMARK_UV_RESPONSE_DN | 0.018 | +0.009 |
  | HALLMARK_ESTROGEN_RESPONSE_EARLY | 0.033 | +0.007 |
  | HALLMARK_UNFOLDED_PROTEIN_RESPONSE | 0.056 | +0.006 |
  | HALLMARK_ESTROGEN_RESPONSE_LATE | 0.038 | +0.005 |
  | HALLMARK_CHOLESTEROL_HOMEOSTASIS | 0.037 | +0.005 |
  | HALLMARK_P53_PATHWAY | 0.028 | +0.005 |
  | HALLMARK_HEDGEHOG_SIGNALING | 0.095 | +0.004 |
  | HALLMARK_PANCREAS_BETA_CELLS | 0.010 | +0.004 |
  | HALLMARK_SPERMATOGENESIS | 0.010 | +0.003 |
  | HALLMARK_APICAL_SURFACE | 0.017 | +0.002 |
  | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | 0.045 | +0.002 |
  | HALLMARK_APOPTOSIS | 0.016 | +0.001 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 1.79 | 14.05 | -0.33 | 0.7453 | 0.6143 |
  | TUBB3 | 1.00 | 0.96 | 1.08 | 15.87 | 0.22 | 0.9506 | 0.3250 |
  | DCX | 0.69 | 0.65 | 0.24 | 12.55 | -0.84 | 1.0000 | 0.2930 |
  | STMN1 | 1.00 | 1.00 | -4.90 | 15.90 | -0.10 | 1.0000 | 0.5440 |
  | STMN2 | 0.19 | 0.53 | -2.16 | 18.21 | 2.95 | 0.0190 | 0.0004 |
  | GAP43 | 0.19 | 0.08 | 1.61 | 12.80 | 0.10 | 0.2380 | 0.6587 |
  | SYN1 | 0.08 | 0.02 | 2.46 | — | — | 0.3237 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.04 | 0.03 | 0.95 | — | — | 0.9008 | — |
  | NEFL | 0.04 | 0.03 | 0.83 | — | — | 0.9277 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.77 | 0.74 | 0.19 | 12.73 | -0.30 | 1.0000 | 0.9769 |
  | UCHL1 | 1.00 | 0.92 | 2.23 | 15.49 | 0.12 | 0.5679 | 0.2930 |
  | SOX2 | 0.23 | 0.11 | 1.37 | 11.82 | -0.55 | 0.2683 | 0.2609 |
  | NES | 0.81 | 0.74 | 0.49 | 15.20 | 0.63 | 0.8601 | 0.0160 |
  | VIM | 0.19 | 0.54 | -2.20 | 12.78 | -4.65 | 0.0173 | 0.0038 |
  | PAX6 | 0.04 | 0.05 | 0.13 | — | — | 1.0000 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 1.00 | 0.95 | 1.47 | 19.65 | 0.46 | 0.9506 | 0.0058 |
  | HOPX | 0.50 | 0.29 | 1.27 | 13.70 | 0.28 | 0.1747 | 0.9140 |
  | EOMES | 0.00 | 0.09 | -2.32 | — | — | 0.4477 | — |
  | TBR1 | 0.19 | 0.10 | 1.18 | 14.00 | 0.40 | 0.4672 | 0.4981 |
  | SATB2 | 0.00 | 0.02 | 0.13 | 16.67 | -0.05 | 1.0000 | 0.6143 |
  | BCL11B | 0.12 | 0.02 | 3.00 | — | — | 0.1233 | — |
  | CUX1 | 0.00 | 0.03 | -0.75 | — | — | 1.0000 | — |
  | CUX2 | 0.12 | 0.06 | 1.16 | — | — | 0.5376 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.15 | 0.01 | 4.60 | — | — | 0.0087 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.96 | 0.74 | 2.56 | 14.44 | -0.77 | 0.0998 | 0.1864 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.00 | 0.12 | -2.81 | — | — | 0.3528 | — |
  | SLC1A3 | 0.00 | 0.03 | -0.49 | — | — | 1.0000 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | 1.19 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.12 | 0.14 | -0.10 | — | — | 1.0000 | — |
  | AIF1 | 0.00 | 0.06 | -1.82 | — | — | 0.7453 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.04 | 0.11 | -1.09 | — | — | 0.6971 | — |
  | COL1A2 | 0.15 | 0.18 | -0.12 | — | — | 1.0000 | — |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.03 | -0.93 | — | — | 1.0000 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.81 | 0.52 | 1.88 | 14.22 | 0.15 | 0.0633 | 0.3913 |
  | TOP2A | 0.81 | 0.59 | 1.44 | 14.70 | 0.45 | 0.1695 | 0.0380 |
  | HMGB2 | 1.00 | 0.82 | 3.57 | 13.71 | -0.24 | 0.0925 | 0.7883 |
  | CDK1 | 0.96 | 0.70 | 2.89 | 14.04 | 0.36 | 0.0368 | 0.1052 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C7

- **n_cells**: 84
- **median_n_detected**: 696.5
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | TUBB2A | 0.37 | 2/3 |
  | RAB5C | 0.36 | 2/3 |
  | RRM1 | 0.33 | 2/3 |
  | FBL | 0.33 | 1/3 |
  | UBE2V2 | 0.29 | 2/3 |
  | ALB | 0.29 | 1/3 |
  | PSMA6 | 0.28 | 2/3 |
  | DDAH1 | 0.27 | 2/3 |
  | GSTM2 | 0.27 | 2/3 |
  | UBE2N | 0.27 | 2/3 |
  | MCM6 | 0.26 | 2/3 |
  | PSAT1 | 0.26 | 2/3 |
  | MCM7 | 0.26 | 2/3 |
  | MCM4 | 0.26 | 2/3 |
  | RPS2 | 0.26 | 2/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | DDAH1 | 1.00 | 0.90 | 1.6e-03 |
  | ACAT2 | 1.00 | 0.92 | 7.5e-03 |
  | HSPB1 | 1.00 | 0.81 | 2.1e-06 |
  | DPYSL5 | 1.00 | 0.88 | 4.6e-04 |
  | FDPS | 1.00 | 0.90 | 3.4e-03 |
  | UCHL1 | 1.00 | 0.91 | 5.1e-03 |
  | MAP2 | 1.00 | 0.93 | 3.4e-02 |
  | PHGDH | 1.00 | 0.90 | 2.4e-03 |
  | PGRMC1 | 1.00 | 0.94 | 4.6e-02 |
  | RAB6A | 1.00 | 0.88 | 4.6e-04 |
  | MAP1B | 1.00 | 0.91 | 5.1e-03 |
  | IDH2 | 1.00 | 0.90 | 2.4e-03 |
  | MCM2 | 1.00 | 0.75 | 1.1e-08 |
  | SRP14 | 1.00 | 0.89 | 1.1e-03 |
  | SLC25A6 | 1.00 | 0.90 | 2.4e-03 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | CKAP4 | 0.10 | 0.56 | 2.8e-14 |
  | UBE2V1 | 0.10 | 0.53 | 1.5e-13 |
  | UBE2N | 0.10 | 0.53 | 1.5e-13 |
  | RBMX | 0.10 | 0.53 | 1.5e-13 |
  | EEF1A1 | 0.10 | 0.53 | 1.5e-13 |
  | KHDRBS1 | 0.10 | 0.53 | 1.8e-13 |
  | RAB5A | 0.10 | 0.53 | 1.9e-13 |
  | ARL3 | 0.25 | 0.68 | 2.5e-12 |
  | CALM1 | 0.10 | 0.52 | 4.3e-13 |
  | KRT13 | 0.10 | 0.52 | 4.5e-13 |
  | DYNLL1 | 0.24 | 0.65 | 2.8e-11 |
  | MAGOH | 0.10 | 0.51 | 2.3e-12 |
  | TP53BP1 | 0.24 | 0.65 | 5.4e-11 |
  | IGF2BP2 | 0.08 | 0.49 | 2.0e-12 |
  | ERP44 | 0.23 | 0.63 | 8.0e-11 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | FBL | +4.11 | 4.3e-07 |
  | FKBP9 | +2.32 | 2.4e-04 |
  | VCL | +2.06 | 1.4e-02 |
  | DNHD1 | +2.03 | 1.9e-02 |
  | CNP | +1.80 | 2.7e-02 |
  | SERPINH1 | +1.80 | 3.2e-02 |
  | CA2 | +1.63 | 6.3e-03 |
  | HNRNPAB | +1.27 | 2.7e-07 |
  | MYDGF | +1.22 | 3.0e-02 |
  | C11orf54 | +1.07 | 1.5e-03 |
  | H1-1 | +1.03 | 2.1e-04 |
  | LMNB2 | +0.96 | 6.7e-05 |
  | HARS1 | +0.94 | 4.0e-02 |
  | INSR | +0.91 | 7.5e-03 |
  | SDCBP2 | +0.82 | 3.0e-02 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | S100B | -3.05 | 2.2e-02 |
  | ALB | -2.95 | 2.6e-08 |
  | DSC1 | -2.55 | 9.2e-06 |
  | SEPTIN5 | -1.89 | 2.0e-02 |
  | SSBP2 | -1.89 | 1.4e-02 |
  | STXBP1 | -1.88 | 4.1e-04 |
  | PER1 | -1.72 | 1.8e-02 |
  | PHF6 | -1.65 | 4.1e-04 |
  | SRSF4 | -1.62 | 4.1e-03 |
  | PSMA6 | -1.51 | 1.0e-14 |
  | MAPRE2 | -1.47 | 1.9e-02 |
  | ATP5PD | -1.46 | 1.7e-13 |
  | TLN1 | -1.40 | 3.6e-04 |
  | TAGLN3 | -1.39 | 3.1e-05 |
  | LARP1 | -1.35 | 3.0e-03 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_E2F_TARGETS | 0.097 | +0.021 |
  | HALLMARK_G2M_CHECKPOINT | 0.102 | +0.015 |
  | HALLMARK_MYC_TARGETS_V1 | 0.190 | +0.013 |
  | HALLMARK_DNA_REPAIR | 0.013 | +0.006 |
  | HALLMARK_ESTROGEN_RESPONSE_LATE | 0.037 | +0.004 |
  | HALLMARK_P53_PATHWAY | 0.026 | +0.002 |
  | HALLMARK_MYC_TARGETS_V2 | 0.024 | +0.002 |
  | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | 0.045 | +0.002 |
  | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITI | 0.029 | +0.002 |
  | HALLMARK_ALLOGRAFT_REJECTION | 0.041 | +0.001 |
  | HALLMARK_HEDGEHOG_SIGNALING | 0.092 | +0.001 |
  | HALLMARK_MITOTIC_SPINDLE | 0.036 | +0.001 |
  | HALLMARK_GLYCOLYSIS | 0.066 | +0.001 |
  | HALLMARK_TGF_BETA_SIGNALING | 0.002 | +0.000 |
  | HALLMARK_INTERFERON_ALPHA_RESPONSE | 0.009 | +0.000 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.93 | 3.58 | 14.59 | 0.29 | 0.0335 | 0.0358 |
  | TUBB3 | 1.00 | 0.96 | 2.87 | 15.58 | -0.10 | 0.1371 | 0.2175 |
  | DCX | 0.71 | 0.64 | 0.46 | 13.59 | 0.25 | 0.3738 | 0.7612 |
  | STMN1 | 1.00 | 1.00 | -3.12 | 16.59 | 0.73 | 1.0000 | 0.0000 |
  | STMN2 | 0.90 | 0.48 | 3.29 | 15.43 | 0.22 | 0.0000 | 0.2135 |
  | GAP43 | 0.02 | 0.09 | -1.66 | — | — | 0.1197 | — |
  | SYN1 | 0.00 | 0.02 | -1.95 | — | — | 0.5510 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.02 | 0.03 | -0.06 | — | — | 1.0000 | — |
  | NEFL | 0.02 | 0.03 | -0.18 | — | — | 1.0000 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.90 | 0.72 | 1.83 | 13.28 | 0.29 | 0.0008 | 0.6543 |
  | UCHL1 | 1.00 | 0.91 | 4.02 | 14.91 | -0.62 | 0.0051 | 0.0000 |
  | SOX2 | 0.00 | 0.13 | -4.60 | — | — | 0.0003 | — |
  | NES | 0.93 | 0.72 | 2.26 | 14.24 | -0.42 | 0.0001 | 0.0000 |
  | VIM | 0.90 | 0.49 | 3.25 | 17.25 | -0.19 | 0.0000 | 0.8235 |
  | PAX6 | 0.00 | 0.06 | -3.34 | — | — | 0.0465 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 1.00 | 0.95 | 3.26 | 18.78 | -0.52 | 0.0668 | 0.0000 |
  | HOPX | 0.11 | 0.32 | -1.91 | 12.66 | -0.79 | 0.0001 | 0.1060 |
  | EOMES | 0.26 | 0.06 | 2.42 | 12.61 | 0.48 | 0.0000 | 0.0335 |
  | TBR1 | 0.00 | 0.12 | -4.47 | — | — | 0.0007 | — |
  | SATB2 | 0.00 | 0.02 | -1.66 | 16.84 | 0.14 | 0.5445 | 0.5686 |
  | BCL11B | 0.01 | 0.02 | -0.35 | — | — | 1.0000 | — |
  | CUX1 | 0.00 | 0.03 | -2.54 | — | — | 0.2765 | — |
  | CUX2 | 0.00 | 0.07 | -3.70 | — | — | 0.0236 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | -1.29 | — | — | 0.7669 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.89 | 0.73 | 1.53 | 15.38 | 0.26 | 0.0042 | 0.0143 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.06 | 0.12 | -0.96 | 11.74 | -3.05 | 0.2597 | 0.0225 |
  | SLC1A3 | 0.01 | 0.03 | -0.60 | — | — | 0.8558 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | -0.59 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.07 | 0.14 | -1.03 | 13.88 | 1.80 | 0.1405 | 0.0273 |
  | AIF1 | 0.01 | 0.07 | -1.97 | — | — | 0.1154 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.17 | 0.10 | 0.85 | 13.67 | -0.22 | 0.1891 | 0.5483 |
  | COL1A2 | 0.02 | 0.19 | -3.00 | — | — | 0.0001 | — |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.04 | -2.71 | — | — | 0.1949 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.92 | 0.48 | 3.48 | 14.45 | 0.42 | 0.0000 | 0.0355 |
  | TOP2A | 0.90 | 0.56 | 2.81 | 14.63 | 0.40 | 0.0000 | 0.0295 |
  | HMGB2 | 0.98 | 0.81 | 3.00 | 14.48 | 0.62 | 0.0001 | 0.0000 |
  | CDK1 | 0.77 | 0.70 | 0.55 | 13.98 | 0.32 | 0.2818 | 0.0005 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C8

- **n_cells**: 25
- **median_n_detected**: 766.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | PEA15 | 0.61 | 2/3 |
  | ELAVL3 | 0.42 | 2/3 |
  | DCX | 0.38 | 3/3 |
  | PSAT1 | 0.37 | 3/3 |
  | LGALS3 | 0.37 | 3/3 |
  | TOP2A | 0.33 | 3/3 |
  | MSI1 | 0.33 | 3/3 |
  | RAB35 | 0.31 | 1/3 |
  | MKI67 | 0.28 | 3/3 |
  | ATOX1 | 0.27 | 2/3 |
  | JPT1 | 0.27 | 3/3 |
  | ARL3 | 0.26 | 2/3 |
  | CEP170 | 0.26 | 1/3 |
  | RAC2 | 0.26 | 1/3 |
  | COTL1 | 0.25 | 2/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | MSI1 | 1.00 | 0.66 | 2.1e-02 |
  | LGALS3 | 1.00 | 0.60 | 2.3e-03 |
  | TOP2A | 1.00 | 0.58 | 2.2e-03 |
  | QKI | 1.00 | 0.68 | 3.0e-02 |
  | MKI67 | 0.96 | 0.51 | 2.1e-03 |
  | PPM1G | 0.76 | 0.39 | 4.0e-02 |
  | SMARCC2 | 0.76 | 0.31 | 3.1e-03 |
  | IQGAP2 | 0.68 | 0.29 | 2.1e-02 |
  | POLR2A | 0.60 | 0.23 | 2.3e-02 |
  | KIF22 | 0.60 | 0.25 | 3.3e-02 |
  | GNA13 | 0.60 | 0.24 | 3.0e-02 |
  | BANF1 | 0.32 | 0.76 | 3.1e-03 |
  | DPYSL4 | 0.24 | 0.65 | 1.3e-02 |
  | NDUFA4 | 0.20 | 0.64 | 6.2e-03 |
  | KIF5B | 0.16 | 0.52 | 4.0e-02 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | UQCRC1 | 0.16 | 0.67 | 1.0e-03 |
  | COTL1 | 0.08 | 0.57 | 1.5e-03 |
  | BANF1 | 0.32 | 0.76 | 3.1e-03 |
  | NDUFA4 | 0.20 | 0.64 | 6.2e-03 |
  | DPYSL4 | 0.24 | 0.65 | 1.3e-02 |
  | IRGQ | 0.12 | 0.53 | 1.2e-02 |
  | RAB35 | 0.00 | 0.40 | 3.1e-03 |
  | KIF5B | 0.16 | 0.52 | 4.0e-02 |
  | CEP170 | 0.00 | 0.36 | 6.2e-03 |
  | RAC2 | 0.00 | 0.36 | 6.2e-03 |
  | COX6B1 | 0.12 | 0.48 | 4.0e-02 |
  | OSBPL3 | 0.08 | 0.42 | 4.0e-02 |
  | DNM1L | 0.08 | 0.42 | 4.0e-02 |
  | GPX1 | 0.04 | 0.38 | 3.0e-02 |
  | KLC1 | 0.00 | 0.31 | 3.0e-02 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | PEA15 | +1.26 | 6.8e-04 |
  | SMC1A | +0.68 | 9.6e-03 |
  | FABP5 | +0.63 | 9.7e-05 |
  | PSAT1 | +0.55 | 4.8e-02 |
  | PHGDH | +0.50 | 1.9e-02 |
  | RPL9P9 | +0.50 | 1.8e-03 |
  | SYNE1 | +0.49 | 4.5e-03 |
  | RPL10A | +0.48 | 1.4e-02 |
  | FAU | +0.46 | 4.5e-03 |
  | RPS25 | +0.45 | 2.2e-03 |
  | RPL23 | +0.45 | 1.4e-02 |
  | RPSA | +0.43 | 5.6e-03 |
  | RPS13 | +0.43 | 1.7e-02 |
  | RPL4 | +0.43 | 3.9e-03 |
  | SRI | +0.42 | 3.2e-02 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | DCX | -1.74 | 5.5e-03 |
  | ELAVL3 | -1.66 | 1.2e-03 |
  | ARPC3 | -1.47 | 5.6e-03 |
  | TSN | -1.45 | 1.5e-02 |
  | MARCKSL1 | -1.37 | 4.2e-02 |
  | ATOX1 | -1.36 | 1.4e-03 |
  | CORO1B | -1.33 | 3.0e-02 |
  | LBR | -1.29 | 5.5e-03 |
  | OPA1 | -1.28 | 4.2e-02 |
  | OSBPL9 | -1.17 | 3.5e-02 |
  | ARPC2 | -1.16 | 1.6e-02 |
  | PSMA6 | -1.15 | 4.4e-02 |
  | WBP11 | -1.13 | 3.0e-02 |
  | NCAM1 | -1.13 | 2.8e-02 |
  | ACTR1A | -1.11 | 5.2e-03 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | PSMC4 | +1.74 | 4.5e-02 |
  | ARCN1 | +1.69 | 6.7e-03 |
  | PEA15 | +1.51 | 1.8e-15 |
  | OGT | +1.43 | 1.2e-02 |
  | PSAT1 | +1.34 | 1.1e-13 |
  | FHL1 | +1.15 | 1.0e-02 |
  | TMED7 | +1.11 | 1.8e-02 |
  | TUBB2A | +1.10 | 3.8e-03 |
  | TUBB2A | +1.10 | 3.8e-03 |
  | QKI | +1.02 | 9.4e-04 |
  | H1-0 | +0.98 | 3.1e-04 |
  | ALDH6A1 | +0.95 | 7.7e-07 |
  | RANBP3 | +0.93 | 5.8e-03 |
  | RPS28 | +0.92 | 3.1e-04 |
  | NDRG2 | +0.92 | 9.6e-07 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | TERF2IP | -7.40 | 2.0e-02 |
  | FKBP9 | -1.93 | 9.0e-04 |
  | CD207 | -1.58 | 3.3e-02 |
  | PRC1 | -1.57 | 3.9e-03 |
  | SMAD4 | -1.53 | 2.5e-02 |
  | DCX | -1.36 | 1.3e-04 |
  | SMARCC2 | -1.26 | 1.5e-02 |
  | SMARCC2 | -1.26 | 1.5e-02 |
  | TMSB10 | -1.20 | 6.1e-05 |
  | ELAVL3 | -1.13 | 1.4e-05 |
  | DCTN2 | -1.10 | 2.6e-02 |
  | CSRP2 | -1.08 | 7.2e-03 |
  | RBP1 | -1.07 | 8.3e-03 |
  | STMN1 | -1.05 | 2.0e-12 |
  | EOMES | -1.05 | 1.8e-02 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_MYC_TARGETS_V1 | 0.209 | +0.032 |
  | HALLMARK_E2F_TARGETS | 0.091 | +0.013 |
  | HALLMARK_G2M_CHECKPOINT | 0.100 | +0.012 |
  | HALLMARK_UNFOLDED_PROTEIN_RESPONSE | 0.061 | +0.011 |
  | HALLMARK_MYC_TARGETS_V2 | 0.032 | +0.010 |
  | HALLMARK_APICAL_SURFACE | 0.023 | +0.008 |
  | HALLMARK_P53_PATHWAY | 0.030 | +0.007 |
  | HALLMARK_APOPTOSIS | 0.021 | +0.007 |
  | HALLMARK_DNA_REPAIR | 0.013 | +0.005 |
  | HALLMARK_ANDROGEN_RESPONSE | 0.022 | +0.005 |
  | HALLMARK_ESTROGEN_RESPONSE_LATE | 0.036 | +0.002 |
  | HALLMARK_NOTCH_SIGNALING | 0.002 | +0.002 |
  | HALLMARK_ALLOGRAFT_REJECTION | 0.041 | +0.002 |
  | HALLMARK_WNT_BETA_CATENIN_SIGNALING | 0.002 | +0.001 |
  | HALLMARK_MYOGENESIS | 0.019 | +0.001 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 1.73 | 13.91 | -0.47 | 1.0000 | 0.0067 |
  | TUBB3 | 1.00 | 0.96 | 1.02 | 14.88 | -0.82 | 1.0000 | 0.0003 |
  | DCX | 0.36 | 0.66 | -1.74 | 11.66 | -1.74 | 0.1549 | 0.0055 |
  | STMN1 | 1.00 | 1.00 | -4.96 | 15.61 | -0.40 | 1.0000 | 0.0363 |
  | STMN2 | 0.76 | 0.52 | 1.50 | 14.53 | -0.81 | 0.3419 | 0.0045 |
  | GAP43 | 0.00 | 0.08 | -2.22 | — | — | 0.9060 | — |
  | SYN1 | 0.00 | 0.02 | -0.11 | — | — | 1.0000 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.00 | 0.03 | -0.70 | — | — | 1.0000 | — |
  | NEFL | 0.00 | 0.03 | -0.81 | — | — | 1.0000 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.80 | 0.73 | 0.43 | 12.53 | -0.54 | 1.0000 | 0.0127 |
  | UCHL1 | 1.00 | 0.92 | 2.17 | 15.26 | -0.12 | 0.9060 | 0.2320 |
  | SOX2 | 0.04 | 0.11 | -1.09 | — | — | 1.0000 | — |
  | NES | 0.84 | 0.74 | 0.78 | 14.47 | -0.13 | 1.0000 | 0.2226 |
  | VIM | 0.80 | 0.52 | 1.78 | 17.03 | -0.40 | 0.2054 | 0.6724 |
  | PAX6 | 0.08 | 0.05 | 1.02 | — | — | 1.0000 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 1.00 | 0.95 | 1.41 | 19.48 | 0.28 | 1.0000 | 0.0590 |
  | HOPX | 0.48 | 0.29 | 1.15 | 13.26 | -0.19 | 0.5870 | 0.6765 |
  | EOMES | 0.08 | 0.08 | 0.22 | — | — | 1.0000 | — |
  | TBR1 | 0.00 | 0.11 | -2.62 | — | — | 0.6638 | — |
  | SATB2 | 0.00 | 0.02 | 0.18 | 16.99 | 0.28 | 1.0000 | 0.6888 |
  | BCL11B | 0.00 | 0.02 | -0.20 | — | — | 1.0000 | — |
  | CUX1 | 0.00 | 0.03 | -0.70 | — | — | 1.0000 | — |
  | CUX2 | 0.00 | 0.07 | -1.85 | — | — | 1.0000 | — |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.00 | 0.01 | 0.55 | — | — | 1.0000 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.92 | 0.74 | 1.69 | 15.52 | 0.35 | 0.5580 | 0.2702 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.00 | 0.12 | -2.75 | — | — | 0.6638 | — |
  | SLC1A3 | 0.00 | 0.03 | -0.43 | — | — | 1.0000 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | 1.25 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.12 | 0.14 | -0.04 | — | — | 1.0000 | — |
  | AIF1 | 0.00 | 0.06 | -1.76 | — | — | 1.0000 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.08 | 0.11 | -0.22 | — | — | 1.0000 | — |
  | COL1A2 | 0.12 | 0.18 | -0.49 | — | — | 1.0000 | — |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.03 | -0.87 | — | — | 1.0000 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.96 | 0.51 | 3.96 | 13.78 | -0.30 | 0.0021 | 0.1909 |
  | TOP2A | 1.00 | 0.58 | 5.18 | 14.37 | 0.06 | 0.0022 | 0.8032 |
  | HMGB2 | 1.00 | 0.82 | 3.51 | 14.36 | 0.43 | 0.2598 | 0.5174 |
  | CDK1 | 0.84 | 0.70 | 1.03 | 13.94 | 0.24 | 0.8410 | 0.2551 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*

## Cluster C9

- **n_cells**: 19
- **median_n_detected**: 859.0
- **consensus markers**:

  | Gene | score | n_modalities |
  |------|-------|--------------|
  | MAGOHB | 0.53 | 3/3 |
  | UBE2N | 0.53 | 3/3 |
  | STMN2 | 0.52 | 3/3 |
  | RAB2A | 0.49 | 3/3 |
  | UBE2V2 | 0.44 | 3/3 |
  | VIM | 0.41 | 3/3 |
  | PHIP | 0.33 | 3/3 |
  | RPS27A | 0.33 | 3/3 |
  | TAGLN3 | 0.30 | 3/3 |
  | RBMX | 0.27 | 3/3 |
  | ACOT7 | 0.27 | 2/3 |
  | SEPTIN3 | 0.26 | 2/3 |
  | KHDRBS1 | 0.25 | 3/3 |
  | H3-4 | 0.23 | 3/3 |
  | TERF2IP | 0.22 | 3/3 |

- **detection markers**:

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | UGDH | 1.00 | 0.65 | 1.3e-02 |
  | RCN2 | 1.00 | 0.74 | 4.7e-02 |
  | ALDH6A1 | 1.00 | 0.75 | 4.7e-02 |
  | CPNE1 | 1.00 | 0.66 | 1.7e-02 |
  | GGH | 1.00 | 0.73 | 4.7e-02 |
  | TMED10 | 1.00 | 0.74 | 4.7e-02 |
  | ISYNA1 | 1.00 | 0.72 | 3.3e-02 |
  | TMED9 | 1.00 | 0.75 | 4.8e-02 |
  | HMGCS1 | 1.00 | 0.73 | 4.7e-02 |
  | EIF1 | 1.00 | 0.74 | 4.7e-02 |
  | MSH2 | 1.00 | 0.52 | 2.1e-03 |
  | GNAI2 | 1.00 | 0.70 | 2.4e-02 |
  | PPP1CA | 1.00 | 0.70 | 2.3e-02 |
  | MAT2B | 1.00 | 0.73 | 4.7e-02 |
  | ECI1 | 1.00 | 0.63 | 1.3e-02 |

- **depleted detection markers** (under-detected vs other clusters):

  | Gene | det_rate_in | det_rate_out | qvalue |
  |------|-------------|--------------|--------|
  | RPS21 | 0.68 | 0.92 | 3.6e-02 |

- **intensity markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | SNW1 | +3.14 | 2.1e-02 |
  | STMN2 | +2.71 | 5.9e-05 |
  | TMSB10 | +2.49 | 3.3e-02 |
  | ELAVL3 | +2.41 | 8.4e-03 |
  | KHDRBS1 | +2.08 | 9.6e-05 |
  | RBMX | +2.07 | 4.4e-05 |
  | MAGOHB | +2.01 | 1.8e-05 |
  | NOVA2 | +2.00 | 1.6e-03 |
  | NCALD | +1.96 | 1.9e-02 |
  | ACOT7 | +1.91 | 8.4e-04 |
  | CENPV | +1.88 | 1.5e-02 |
  | TAGLN3 | +1.87 | 9.4e-04 |
  | RAB2A | +1.81 | 9.6e-05 |
  | OTUB1 | +1.72 | 6.8e-03 |
  | HSPBP1 | +1.65 | 7.4e-03 |

- **depleted intensity markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | VIM | -4.01 | 7.2e-03 |
  | HMGA2 | -1.65 | 3.6e-03 |
  | CTNNB1 | -1.05 | 4.1e-02 |
  | FLG2 | -1.03 | 4.8e-02 |
  | IGF2BP1 | -0.92 | 4.5e-02 |
  | ANXA6 | -0.69 | 8.9e-03 |
  | SMC4 | -0.52 | 2.4e-02 |
  | RRM1 | -0.43 | 4.3e-02 |
  | SHMT2 | -0.37 | 3.0e-02 |

- **scplainer markers** (up-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | NUCKS1 | +2.95 | 4.1e-02 |
  | NFIX | +2.90 | 6.1e-07 |
  | TAGLN3 | +2.55 | 1.7e-10 |
  | BRD3 | +2.52 | 7.4e-03 |
  | LSM4 | +2.48 | 8.4e-03 |
  | PSMC4 | +2.43 | 1.0e-04 |
  | H1-2 | +2.43 | 8.3e-04 |
  | RAB2A | +2.40 | 9.9e-23 |
  | STMN2 | +2.39 | 1.3e-15 |
  | NOVA2 | +2.35 | 4.5e-04 |
  | NOVA2 | +2.35 | 4.5e-04 |
  | WTAP | +2.28 | 3.4e-02 |
  | GNA13 | +2.26 | 3.8e-08 |
  | SEPTIN3 | +2.25 | 7.4e-12 |
  | PDHX | +2.22 | 1.7e-04 |

- **depleted scplainer markers** (down-regulated):

  | Gene | log2FC | qvalue |
  |------|--------|--------|
  | TERF2IP | -5.72 | 6.1e-06 |
  | VIM | -3.31 | 7.4e-09 |
  | VIM | -3.31 | 7.4e-09 |
  | VCL | -2.57 | 1.1e-02 |
  | RAB21 | -1.87 | 8.3e-04 |
  | C11orf54 | -1.59 | 9.1e-04 |
  | FKBP9 | -1.50 | 8.1e-04 |
  | RBM42 | -1.32 | 1.3e-02 |
  | RNASEH2C | -1.29 | 1.1e-02 |
  | GTF3C5 | -1.20 | 2.5e-03 |
  | SOGA1 | -1.20 | 5.0e-02 |
  | SMARCC2 | -1.16 | 4.9e-02 |
  | SMARCC2 | -1.16 | 4.9e-02 |
  | PCNA | -1.08 | 2.8e-04 |
  | KRT25 | -1.04 | 2.2e-03 |

- **AUCell pathways**:

  | Pathway | AUC_in | delta_vs_rest |
  |---------|--------|---------------|
  | HALLMARK_E2F_TARGETS | 0.096 | +0.018 |
  | HALLMARK_HEDGEHOG_SIGNALING | 0.107 | +0.016 |
  | HALLMARK_G2M_CHECKPOINT | 0.104 | +0.016 |
  | HALLMARK_MYC_TARGETS_V1 | 0.193 | +0.015 |
  | HALLMARK_ESTROGEN_RESPONSE_EARLY | 0.034 | +0.009 |
  | HALLMARK_PANCREAS_BETA_CELLS | 0.014 | +0.007 |
  | HALLMARK_ESTROGEN_RESPONSE_LATE | 0.037 | +0.003 |
  | HALLMARK_UV_RESPONSE_DN | 0.013 | +0.003 |
  | HALLMARK_ANDROGEN_RESPONSE | 0.020 | +0.003 |
  | HALLMARK_CHOLESTEROL_HOMEOSTASIS | 0.034 | +0.002 |
  | HALLMARK_SPERMATOGENESIS | 0.010 | +0.002 |
  | HALLMARK_WNT_BETA_CATENIN_SIGNALING | 0.003 | +0.002 |
  | HALLMARK_ADIPOGENESIS | 0.010 | +0.002 |
  | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | 0.045 | +0.001 |
  | HALLMARK_HEME_METABOLISM | 0.015 | +0.000 |

- **researcher markers of interest** (per-marker stats):

  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |
  |------|-------------|--------------|---------|-------------------|------------|----------|----------|
  | MAP2 | 1.00 | 0.94 | 1.33 | 14.25 | -0.10 | 0.9678 | 0.6544 |
  | TUBB3 | 1.00 | 0.96 | 0.63 | 15.69 | 0.02 | 1.0000 | 0.2125 |
  | DCX | 0.68 | 0.65 | 0.17 | 15.85 | 2.47 | 1.0000 | 0.3625 |
  | STMN1 | 1.00 | 1.00 | -5.36 | 15.85 | -0.14 | 1.0000 | 0.2027 |
  | STMN2 | 0.37 | 0.53 | -0.89 | 17.97 | 2.71 | 0.5237 | 0.0001 |
  | GAP43 | 0.21 | 0.08 | 1.78 | — | — | 0.1959 | — |
  | SYN1 | 0.11 | 0.02 | 2.95 | — | — | 0.1774 | — |
  | SNAP25 | — | — | — | — | — | — | — |
  | NEFM | 0.16 | 0.03 | 2.94 | — | — | 0.0813 | — |
  | NEFL | 0.11 | 0.03 | 2.18 | — | — | 0.3233 | — |
  | NEFH | — | — | — | — | — | — | — |
  | ENO2 | 0.74 | 0.74 | -0.08 | 13.81 | 0.78 | 1.0000 | 0.5060 |
  | UCHL1 | 1.00 | 0.92 | 1.77 | 15.33 | -0.05 | 0.7167 | 0.1898 |
  | SOX2 | 0.32 | 0.11 | 1.99 | 12.90 | 0.62 | 0.0749 | 0.0487 |
  | NES | 0.68 | 0.74 | -0.45 | 15.13 | 0.55 | 0.9678 | 0.2051 |
  | VIM | 0.37 | 0.53 | -0.92 | 13.42 | -4.01 | 0.4065 | 0.0072 |
  | PAX6 | 0.21 | 0.05 | 2.56 | — | — | 0.0703 | — |
  | HES1 | — | — | — | — | — | — | — |
  | FABP7 | 1.00 | 0.95 | 1.02 | 19.64 | 0.44 | 1.0000 | 0.0195 |
  | HOPX | 0.42 | 0.30 | 0.81 | 13.25 | -0.20 | 0.6195 | 0.6911 |
  | EOMES | 0.00 | 0.09 | -1.87 | — | — | 0.7167 | — |
  | TBR1 | 0.37 | 0.10 | 2.47 | 13.83 | 0.34 | 0.0245 | 0.1538 |
  | SATB2 | 0.00 | 0.02 | 0.58 | 16.34 | -0.38 | 1.0000 | 0.6684 |
  | BCL11B | 0.11 | 0.02 | 2.86 | — | — | 0.1918 | — |
  | CUX1 | 0.00 | 0.03 | -0.30 | — | — | 1.0000 | — |
  | CUX2 | 0.32 | 0.06 | 2.97 | 13.03 | -0.05 | 0.0161 | 0.9631 |
  | RELN | — | — | — | — | — | — | — |
  | DLX1 | — | — | — | — | — | — | — |
  | DLX2 | 0.11 | 0.01 | 3.74 | — | — | 0.0946 | — |
  | GAD1 | — | — | — | — | — | — | — |
  | GAD2 | — | — | — | — | — | — | — |
  | GFAP | 0.95 | 0.74 | 2.08 | 15.33 | 0.14 | 0.1918 | 0.8653 |
  | AQP4 | — | — | — | — | — | — | — |
  | S100B | 0.00 | 0.12 | -2.35 | — | — | 0.3820 | — |
  | SLC1A3 | 0.00 | 0.03 | -0.03 | — | — | 1.0000 | — |
  | SOX10 | — | — | — | — | — | — | — |
  | OLIG1 | — | — | — | — | — | — | — |
  | OLIG2 | — | — | — | — | — | — | — |
  | PLP1 | — | — | — | — | — | — | — |
  | MBP | 0.00 | 0.01 | 1.65 | — | — | 1.0000 | — |
  | MOG | — | — | — | — | — | — | — |
  | CNP | 0.00 | 0.14 | -2.67 | — | — | 0.2732 | — |
  | AIF1 | 0.00 | 0.06 | -1.36 | — | — | 0.9678 | — |
  | CSF1R | — | — | — | — | — | — | — |
  | CX3CR1 | — | — | — | — | — | — | — |
  | TREM2 | — | — | — | — | — | — | — |
  | CD68 | — | — | — | — | — | — | — |
  | PECAM1 | — | — | — | — | — | — | — |
  | VWF | — | — | — | — | — | — | — |
  | CDH5 | — | — | — | — | — | — | — |
  | CLDN5 | — | — | — | — | — | — | — |
  | COL1A1 | 0.11 | 0.11 | 0.22 | — | — | 1.0000 | — |
  | COL1A2 | 0.21 | 0.18 | 0.44 | — | — | 1.0000 | — |
  | LUM | — | — | — | — | — | — | — |
  | DCN | — | — | — | — | — | — | — |
  | PDGFRB | 0.00 | 0.03 | -0.47 | — | — | 1.0000 | — |
  | ACTA2 | — | — | — | — | — | — | — |
  | MKI67 | 0.58 | 0.52 | 0.30 | 14.20 | 0.14 | 0.9805 | 0.7387 |
  | TOP2A | 0.58 | 0.60 | -0.13 | 14.70 | 0.43 | 1.0000 | 0.1397 |
  | HMGB2 | 1.00 | 0.82 | 3.11 | 13.75 | -0.19 | 0.1340 | 0.3746 |
  | CDK1 | 0.84 | 0.70 | 1.01 | 13.74 | 0.03 | 0.4818 | 0.7944 |

  *det_rate: fraction of cells with protein detected (0–1); log2_OR: log2 odds ratio for detection specificity; median_log2int_in: median log2(intensity+1) in detected cells of this cluster; log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); — = not tested (below min_cells threshold or gene not matched)*
