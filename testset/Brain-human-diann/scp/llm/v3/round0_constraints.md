# Round 0 constraints — brain

*Generated: 2026-03-07T20:37:10.225989*
*Prompt hash: 0b1a5b32289b4955*

## 1. Developmental / tissue vocabulary constraints

**Do NOT use mature adult brain labels** that imply terminal differentiation or postnatal states, because GW13–19 cortex/telencephalon is dominated by progenitors and immature neurons/glia. Avoid: *mature astrocyte* (homeostatic ALDH1L1/AQP4-rich adult-like), *mature oligodendrocyte / myelinating OL* (MBP/PLP1-high myelin program), *fully mature microglia* (adult homeostatic framing), *mature endothelial/pericyte* (adult BBB specialization language), and adult neuronal subtype labels that presume fully established circuitry (e.g., “IT L2/3”, “PT”, “Corticothalamic” as mature classes).

**Prefer developmental terminology**, e.g.:
- **Neural progenitors**: neuroepithelial-like / **apical radial glia (aRG)**, **basal/intermediate progenitors (IP)**, **outer radial glia (oRG)** (especially relevant in human GZ).
- **Neuronal lineage**: **newborn neurons / neuroblasts**, **immature excitatory neurons** (CP-enriched), **immature interneurons** (notably from ganglionic eminences at GW13), **migrating interneurons**.
- **Glial lineage**: **glial progenitors / gliogenic radial glia**, **OPC-like** (precursor) rather than mature OL; **astroglial precursor** rather than mature astrocyte.
- **Non-neural**: **endothelial**, **pericyte**, **vascular smooth muscle–like**, **meningeal fibroblast–like**, **choroid plexus/ependymal-like** (developmental epithelia), **microglia** (developmental/amoeboid framing is acceptable).

**Expected compartments by dissection**: GZ should enrich progenitors (RG/IP/oRG) and cycling states; CP should enrich post-mitotic immature excitatory neurons and early glia. GW13 whole telencephalon including ganglionic eminences increases expectation of interneuron lineage and migrating interneurons.

---

## 2. Expected ambient signals

Fresh prenatal brain dissociation and single-cell proteomics commonly yield **ambient carryover** from abundant cytosolic and blood/ECM proteins that can appear across many clusters and be **non-discriminative**:
- **Ribosomal proteins** (RPL/RPS), **core histones** (H2A/H2B/H3/H4), **actin/tubulin** (ACTB, TUBB), **GAPDH/ENO1/PKM**, **HSPs** (HSPA8/HSP90), **ubiquitin/proteasome components**.
- **Mitochondrial proteins** (ATP synthase, VDACs) from stressed/lysed cells.
- **Plasma/blood contaminants** depending on perfusion: **albumin (ALB)**, **hemoglobins (HBA/HBB)**, **transferrin (TF)**, **apolipoproteins**—often broadly detected and not cell-type specific.
- **ECM/meningeal proteins** (collagens, fibronectin) can be ambient if meninges/vasculature are present in dissections.

Treat these as background unless strongly enriched in a biologically coherent population.

---

## 3. Non-self protein acquisition

Populations expected here that can **acquire non-self proteins**:
- **Microglia / macrophage-like cells**: phagocytosis of apoptotic progenitors/neurons and synaptic material; would show mixed neuronal/progenitor proteins alongside myeloid programs.
- **Endothelial/perivascular cells**: uptake/transcytosis of plasma proteins (albumin, transferrin, apolipoproteins), especially if blood contamination exists.
- **Radial glia/progenitors** can contain engulfed debris during development (less than microglia but possible), producing low-level “foreign” neuronal markers.

Mechanisms: **phagocytosis** (dominant), **trogocytosis** (immune-like nibbling; possible for microglia), and less likely **NETosis** (neutrophils are not expected in healthy prenatal brain; if present, consider contamination).

In data, this appears as **co-detection of lineage-incongruent proteins** (e.g., neuronal structural proteins within microglia-like profiles) without corresponding full transcriptional/proteomic programs of that lineage.

---

## 4. Expected artefacts

Common artefacts in fresh prenatal brain SCP:
- **Dissociation stress / ischemia**: elevated heat-shock/chaperones, immediate-early/stress proteins, proteolysis signatures; can form “stress clusters” or blur boundaries.
- **Cell cycle effects**: strong separation of cycling vs non-cycling progenitors (S/G2M) that can dominate clustering; interpret as states within progenitor lineages.
- **Batch/dissection effects**: GW13 vs GW15/19 and GZ vs CP can drive global shifts (developmental time, region, cell size/protein content), potentially clustering by sample rather than biology.
- **Doublets/multiplets**: especially in dense fetal tissue; manifest as hybrid profiles (e.g., progenitor + neuron, neuron + endothelial) with unusually broad marker sets.
- **Ambient/blood contamination gradients**: clusters may separate by degree of ALB/HBB/ECM carryover rather than true identity.
- **Low-protein/low-coverage cells**: can cluster by missingness, appearing as ambiguous “low signal” groups.