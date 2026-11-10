"""Export Tabula Muris Senis droplet pancreas h5ad (raw UMI counts + official annotations) to MTX for Seurat."""
import sys, os
sys.stdout.reconfigure(encoding="utf-8")
import scanpy as sc, scipy.io, scipy.sparse as sp, pandas as pd, numpy as np
src = sys.argv[1]; out = sys.argv[2]; os.makedirs(out, exist_ok=True)
a = sc.read_h5ad(src)
print("adata", a.shape, "raw:", None if a.raw is None else a.raw.shape)
print("obs columns:", list(a.obs.columns))
X = a.raw.X if a.raw is not None else a.X
genes = (a.raw.var_names if a.raw is not None else a.var_names)
X = sp.csr_matrix(X)
print("value range", X.min(), X.max(), "integer-like:", bool(np.allclose(X.data[:10000], np.round(X.data[:10000]))))
scipy.io.mmwrite(os.path.join(out, "counts.mtx"), X.T.tocoo())  # genes x cells
pd.Series(genes).to_csv(os.path.join(out, "genes.tsv"), index=False, header=False)
keep = [c for c in ["cell_ontology_class", "free_annotation", "age", "sex", "mouse.id", "tissue", "method", "n_genes", "n_counts"] if c in a.obs.columns]
obs = a.obs[keep].copy(); obs.index.name = "cell"
obs.to_csv(os.path.join(out, "meta.tsv"), sep="\t")
print(obs["cell_ontology_class"].value_counts().to_string())
print("ages:", obs["age"].value_counts().to_dict() if "age" in obs else None)
