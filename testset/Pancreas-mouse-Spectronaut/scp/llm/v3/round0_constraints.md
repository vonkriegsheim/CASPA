# Round 0 constraints (auto-generated)

*Generated: 2026-03-09T10:31:41.481891+00:00*
*Model: gpt-5.2*

## 1. Developmental / tissue vocabulary constraints
- **Do NOT use developmental/embryonic progenitor labels** (e.g., “endocrine progenitor,” “ductal progenitor,” “acinar progenitor,” “multipotent pancreatic progenitor,” “Ngn3+ endocrine precursor”). This is an **adult injury/regeneration** setting (2 days post-caerulein), so observed “immature-like” states should be framed as **injury-associated / dedifferentiated / regenerating** rather than embryonic.
- Avoid **fully mature islet subtype overcommitment** unless strongly supported (e.g., “gamma/PP cell,” “epsilon cell”) because endocrine cells may be under-sampled in pancreas dissociations and SCP coverage can be sparse.
- Prefer terminology such as:
  - **“injured/stressed acinar,” “acinar-to-ductal metaplasia (ADM)-like,” “duct-like/reactive ductal,” “regenerating acinar,” “inflammatory epithelial program,” “cycling epithelial”**.
  - For mesenchyme: **“activated fibroblast / myofibroblast (stellate cell–like),” “ECM-remodeling fibroblast”** rather than developmental mesenchymal progenitors.
- **Expected cell types (biologically plausible at day 2 caerulein):**
  - **Acinar** (including stressed/injured states, partial dedifferentiation/ADM-like programs)
  - **Ductal** (reactive/expanded)
  - **Endocrine** (beta/alpha/delta potentially present but not necessarily abundant)
  - **Immune infiltrates**: **neutrophils, monocytes/macrophages, dendritic cells, T cells, B cells, NK cells**
  - **Stromal**: **pancreatic stellate cells/activated fibroblasts, endothelial cells, pericytes/smooth muscle**
  - Potential **acinar debris/low-content** events due to injury.

## 2. Expected ambient signals
Likely **non-discriminative** proteins appearing across many/all clusters due to tissue lysis, high-abundance pancreas proteins, and blood/serum carryover:
- **Acinar digestive enzymes** as pervasive ambient: **PRSS1/2 (trypsinogens), CPA1/CPA2, CPB1, CTRB1/2, CEL, AMY2, REG family (e.g., REG3)**. In pancreatitis models these can be abundant extracellularly and stick to cells.
- **Cytoskeletal/housekeeping**: **ACTB, TUBB, GAPDH, ENO1**, ribosomal proteins—often reflect general content/contamination more than identity.
- **Mitochondrial proteins** can look globally elevated with stress/low-quality events (not cell-type specific).
- **Blood-derived ambient** (depending on perfusion): **ALB, hemoglobin subunits (HBA/HBB), transferrin**, complement proteins—generally not informative for pancreas cell identity.
Interpretation constraint: treat “acinar enzyme positivity” in immune/stromal clusters cautiously as likely ambient or uptake (see below), not true expression.

## 3. Non-self protein acquisition
Populations expected to **acquire proteins they did not synthesize**:
- **Macrophages/monocytes/DCs**: phagocytosis/efferocytosis of dying acinar cells → detectable **acinar enzymes + lysosomal proteins (LAMP1/2, CTSD/CTSB)** together. Would appear as immune markers coexisting with strong acinar cargo.
- **Neutrophils**: uptake of debris; plus possible **degranulation/NET-associated proteins** (e.g., MPO, ELANE, S100A8/A9) that can also become ambient and adhere to other cells.
- **Stellate cells/fibroblasts** may internalize extracellular proteins during injury remodeling, but typically less dramatic than professional phagocytes.
Data signature: **mixed lineage protein sets** (immune identity markers plus acinar digestive enzymes) without coherent epithelial programs; often coupled with high endo/lysosomal pathways.

## 4. Expected artefacts
Common SCP + injured pancreas artefacts and how they manifest:
- **Dissociation stress / injury response inflation**: broad upregulation of heat shock/chaperones (**HSP90/HSPA**), immediate-early response proteins, oxidative stress signatures across multiple clusters, potentially forming “stress clusters” not true cell types.
- **Ambient acinar enzyme dominance**: can blur boundaries and pull unrelated cells together by shared high-abundance enzyme signal; may create misleading “acinar-like” clusters among immune/stromal cells.
- **Doublets/multiplets**: especially **acinar–immune** or **acinar–ductal** mixtures in inflamed tissue; appear as simultaneous strong epithelial and immune/stromal marker panels with elevated total protein intensity/features.
- **Batch / run effects**: SCP is sensitive to plate position, carrier effects, and instrument drift; clusters may segregate by run/plate rather than biology, often with systematic shifts in global intensity, missingness, or mitochondrial/ribosomal proportions.
- **Low-content/debris events**: injured pancreas yields fragile acinar cells; can produce clusters with low protein diversity dominated by a few abundant proteins (actin/tubulin/enzymes), representing dying cells or debris rather than a distinct population.