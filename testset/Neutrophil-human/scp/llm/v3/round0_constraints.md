# Round 0 constraints (auto-generated)

*Generated: 2026-03-10T09:26:24.473890+00:00*
*Model: gpt-5.2*

### 1. Developmental / tissue vocabulary constraints
- **Do NOT label as non-neutrophil mature cell types**: neuron, astrocyte, oligodendrocyte, microglia/macrophage, endothelial, pericyte, fibroblast, epithelial/keratinocyte, T/NK/B cell, dendritic cell, monocyte, platelet, RBC. The FACS gate (CD45+CD66b+CD49d-) makes these biologically implausible as true identities; such proteins must be interpreted as **cargo/contact/ambient**.
- **Preferred neutrophil-state terminology (use “state” not “type”)**:  
  - *Tumour-associated neutrophil (TAN)* as the umbrella label.  
  - *Immature / banded-like / granulopoiesis-skewed* (reflecting developmental granule program carryover rather than lineage mixing).  
  - *Activated/primed/armed*, *degranulating*, *phagocytic*, *immunosuppressive-like*, *NET-forming (vital)*, *NETosis/lytic remnant*, *IFN-stimulated*, *hypoxic/metabolic-stressed*, *antigen-presentation–associated (HLA-rich)* (still neutrophils).  
- **Expected populations**: multiple TAN functional states plus **patient/batch-dominated** separations. “Contaminant cell types” should only be invoked as **acquired material** (Section 3), not as labels.

### 2. Expected ambient signals
Enzymatic digestion + sorting of highly granulated neutrophils causes widespread extracellular release and re-uptake/adsorption of abundant proteins. Likely **ambient/non-discriminative across essentially all cells** include:
- **Neutrophil granule proteins**: MPO, ELANE, PRTN3, CTSG, AZU1, LTF, CAMP, BPI, DEFA1/3, S100A8, S100A9, S100A12, LGALS3, MMP9 (and related granule/secretory proteins). Presence/detection is expected everywhere; **only relative intensity patterns** (enrichment vs depletion) should be used for state calls.
- **Common plasma/ECM carryover** from tumor vasculature and dissociation: ALB, APOA1/APOA2, APOE, TF, FGA/FGB/FGG, HP, HPR, SERPINA1, SERPINC1, C3/C4, FN1, various collagens. These are typically **non-cell-type-discriminative** in this gated dataset.
- **Housekeeping stress/background**: actins/tubulins, histones at low levels (but see NETosis below), ribosomal proteins. Treat as technical/contextual unless strongly patterned.

### 3. Non-self protein acquisition
Neutrophils can display proteins not synthesized by them due to:
- **Phagocytosis**: uptake of tumor/brain debris, RBCs, platelets, and microbial-like material (if any). Data signature: co-detection of neutrophil granule proteins with **tissue-specific cargos** (e.g., hemoglobins; neural proteins; keratins) plus lysosomal/proteolytic machinery (CTSD/CTSS, LAMPs) and sometimes reduced “cytosolic purity”.
- **Trogocytosis / membrane nibbling**: acquisition of surface/membrane proteins from tumor, endothelial, or immune cells. Data signature: patchy presence of **membrane-associated markers** atypical for neutrophils (e.g., EPCAM/keratins as contact proxies; endothelial adhesion proteins) without a full non-neutrophil proteome.
- **NET formation (vital) and NETosis (lytic remnants)**: extracellular DNA/histone complexes and granule proteins. Data signature: elevated histones (H2A/H2B/H3/H4), chromatin proteins, sometimes PAD4-associated patterns, plus granule proteins; “remnant” profiles may show loss of typical cytosolic proteins with persistence of chromatin/granule components.
- **Platelet cloaking / aggregates**: neutrophil–platelet interactions can add platelet proteins (PF4, PPBP, GP1BA, ITGA2B). Interpret as interaction state or doublet/aggregate depending on breadth and intensity.

### 4. Expected artefacts
- **Dissociation/handling activation**: enzymatic digestion and sorting can induce rapid degranulation/priming. Manifests as broad shifts in granule protein intensities, stress proteins (HSPs), oxidative enzymes, and adhesion-related proteins—potentially forming clusters reflecting “processing-induced activation” rather than in vivo TAN biology.
- **Ambient granule contamination**: pervasive granule proteins flatten biological separation if using detection rate; clustering may instead reflect subtle intensity differences or missingness. Avoid over-interpreting ubiquitous MPO/ELANE/S100A8/A9 presence.
- **Batch/patient effects (6 patients = 6 batches)**: strong between-patient proteomic baselines (therapy, steroid use, ischemia time, tumor necrosis, blood content) can drive clusters. Expect clusters enriched for one patient’s global intensity scaling, plasma contamination, or ischemia/stress.
- **Doublets/aggregates**: despite gating, neutrophil–platelet or neutrophil–tumor debris aggregates can appear. Signature: simultaneous strong platelet/plasma/tissue proteins plus intact neutrophil program, unusually high total protein intensity, and broadened proteome complexity.
- **Cell damage/lysis and “NET remnants” vs low-quality cells**: ruptured cells can look like NETosis (histone-rich, cytosol-poor). Distinguish by consistency of NET-associated proteins vs random loss and extreme missingness.