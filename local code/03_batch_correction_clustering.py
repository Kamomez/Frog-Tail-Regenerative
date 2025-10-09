import time
from pathlib import Path
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, rand_score, silhouette_score
import numpy as np
import pandas as pd
import scanpy as sc

def add_clusters(adata_obj, use_rep=None, suffix=""):
    """Add multiple clustering results to AnnData object"""
    # Compute neighbor graph for graph-based clustering
    sc.pp.neighbors(adata_obj, use_rep=use_rep, n_neighbors=15)
    
    # Leiden clustering
    sc.tl.leiden(adata_obj, resolution=1.0, key_added=f"leiden{suffix}")
    
    # Louvain clustering
    try:
        sc.tl.louvain(adata_obj, resolution=1.0, key_added=f"louvain{suffix}")
    except (ImportError, ModuleNotFoundError) as e:
        print(f"Warning: Louvain clustering skipped ({e}). Using Leiden only (Leiden is an improved version of Louvain).")

    # K-means clustering
    emb = adata_obj.obsm.get("X_pca") if use_rep is None else adata_obj.obsm.get(use_rep)
    if emb is not None:
        km = KMeans(n_clusters=12, random_state=7, n_init="auto").fit_predict(emb[:, :30])
        adata_obj.obs[f"kmeans{suffix}"] = km.astype(str)

def compute_metrics(adata_obj, cluster_key, embedding_key="X_pca"):
    """Calculate clustering evaluation metrics"""
    if embedding_key not in adata_obj.obsm:
        embedding_key = "X_pca"
    emb = adata_obj.obsm.get(embedding_key)
    labels = adata_obj.obs[cluster_key].astype(str)

    if emb is None or labels.nunique() < 2:
        return {"silhouette": np.nan, "ari": np.nan, "rand": np.nan}

    try:
        sil = float(silhouette_score(emb, labels))
    except ValueError:
        sil = np.nan

    if "cluster" in adata_obj.obs:
        ref = adata_obj.obs["cluster"].astype(str)
        ari = float(adjusted_rand_score(ref, labels))
        rand = float(rand_score(ref, labels))
    else:
        ari = np.nan
        rand = np.nan

    return {"silhouette": sil, "ari": ari, "rand": rand}

def main():
    """Main function"""
    print("Batch Correction and Clustering Analysis")

    # Create necessary directories
    Path("results").mkdir(exist_ok=True)
    Path("figures").mkdir(exist_ok=True)

    # Check if denoised data exists
    baseline_file = Path("results/baseline_denoised.h5ad")
    if not baseline_file.exists():
        print("Error: run 02_data_denoising.py first")
        return

    # Load data
    print("Loading data...")
    datasets = {}

    # Load baseline data
    adata_hvg = sc.read_h5ad(baseline_file)
    datasets[""] = (adata_hvg, None)

    # Load MAGIC data (if exists)
    magic_file = Path("results/magic_denoised.h5ad")
    if magic_file.exists():
        datasets["magic"] = (sc.read_h5ad(magic_file), None)

    # Load ALRA data (if exists)
    alra_file = Path("results/alra_denoised.h5ad")
    if alra_file.exists():
        datasets["alra"] = (sc.read_h5ad(alra_file), "X_alra")

    # Check batch information
    batch_key = 'DaysPostAmputation' if 'DaysPostAmputation' in adata_hvg.obs.columns else None
    if batch_key:
        print(f"Using batch key: {batch_key}")
        print(f"Batch distribution: {adata_hvg.obs[batch_key].value_counts()}")
    else:
        print("No batch information found")

    # 3.1 Check batch correction library availability
    print("\n3.1 Check batch correction library availability")

    try:
        from scanpy.external.pp import bbknn
        BBKNN_AVAILABLE = True
        print("[OK] BBKNN available")
    except ImportError:
        BBKNN_AVAILABLE = False
        print("[FAIL] BBKNN unavailable")

    try:
        from scanpy.external.pp import harmony_integrate
        HARMONY_AVAILABLE = True
        print("[OK] Harmony available")
    except ImportError:
        HARMONY_AVAILABLE = False
        print("[FAIL] Harmony unavailable")

    try:
        import scanorama
        SCANORAMA_AVAILABLE = True
        print("[OK] Scanorama available")
    except ImportError:
        SCANORAMA_AVAILABLE = False
        print("[FAIL] Scanorama unavailable")

    # 3.2 BBKNN Batch Correction
    print("\n3.2 BBKNN Batch Correction")

    bbknn_data = None
    if BBKNN_AVAILABLE and batch_key:
        # Create BBKNN data copy
        bbknn_data = adata_hvg.copy()

        print("Starting BBKNN batch correction...")
        print("Note: BBKNN may be slow on large datasets. Consider using Harmony instead.")
        start_time = time.time()

        try:
            # Apply BBKNN with optimized parameters
            bbknn(bbknn_data, batch_key=batch_key, neighbors_within_batch=3, n_pcs=20)

            # Calculate UMAP
            sc.tl.umap(bbknn_data, random_state=0, n_components=2, maxiter=100)

            elapsed_time = time.time() - start_time
            print(f"BBKNN batch correction completed, time taken: {elapsed_time:.2f} seconds")
            datasets["bbknn"] = (bbknn_data, None)
        except KeyboardInterrupt:
            print("BBKNN interrupted by user, skipping...")
            bbknn_data = None
        except Exception as e:
            print(f"BBKNN batch correction failed: {e}")
            bbknn_data = None
    else:
        print("Skipping BBKNN batch correction")

    # 3.3 Scanorama Batch Correction
    print("\n3.3 Scanorama Batch Correction")

    corrected_scanorama = None
    if SCANORAMA_AVAILABLE and batch_key:
        try:
            from scanpy.external.pp import scanorama_integrate
            
            # Create Scanorama data copy and sort by batch (required for scanorama_integrate)
            scanorama_data = adata_hvg.copy()
            scanorama_data = scanorama_data[scanorama_data.obs[batch_key].argsort()].copy()

            print("Starting Scanorama batch correction...")
            print(f"Data shape: {scanorama_data.shape}")
            print(f"Number of batches: {scanorama_data.obs[batch_key].nunique()}")
            start_time = time.time()

            # Use scanpy's scanorama_integrate
            # This method is much faster and more stable than scanorama.correct()
            scanorama_integrate(scanorama_data, key=batch_key, basis='X_pca', adjusted_basis='X_scanorama')
            
            corrected_scanorama = scanorama_data
            
            # Calculate UMAP using the integrated representation
            sc.pp.neighbors(corrected_scanorama, n_neighbors=10, use_rep='X_scanorama')
            sc.tl.umap(corrected_scanorama, random_state=0, n_components=2, maxiter=100)

            elapsed_time = time.time() - start_time
            print(f"Scanorama batch correction completed, time taken: {elapsed_time:.2f} seconds")
            datasets["scanorama"] = (corrected_scanorama, 'X_scanorama')

        except Exception as e:
            print(f"Scanorama batch correction failed: {e}")
            import traceback
            traceback.print_exc()
            corrected_scanorama = None
    else:
        print("Skipping Scanorama batch correction")

    # 3.4 Simplified Batch Correction (Linear Regression)
    print("\n3.4 Simplified Batch Correction (Linear Regression)")
    
    linear_corrected = None
    if batch_key:
        print("Starting simplified batch correction (linear regression)...")
        #print("Note: This is a basic method that corrects batch effects using linear regression.")
        start_time = time.time()

        try:
            from sklearn.linear_model import LinearRegression

            linear_corrected = adata_hvg.copy()
            X = linear_corrected.X.toarray() if hasattr(linear_corrected.X, 'toarray') else linear_corrected.X

            print(f"Correcting batch effects for {X.shape[1]} genes...")
            
            # Correct batch effects for each gene
            batch_indicators = pd.get_dummies(linear_corrected.obs[batch_key]).values
            
            # Use first batch as reference
            reference_batch = batch_indicators[:, 0] == 1

            X_corrected = np.zeros_like(X)

            # Correct all genes
            for gene_idx in range(X.shape[1]):
                if gene_idx % 500 == 0:
                    print(f"  Progress: {gene_idx}/{X.shape[1]} genes corrected")
                
                gene_expr = X[:, gene_idx]

                # Build linear model: gene expression = batch effect + residual
                lr = LinearRegression(fit_intercept=True)
                lr.fit(batch_indicators, gene_expr)

                # Remove batch effects (preserve reference batch level)
                predicted_batch_effects = lr.predict(batch_indicators)
                reference_mean = np.mean(gene_expr[reference_batch]) if np.sum(reference_batch) > 0 else 0

                X_corrected[:, gene_idx] = gene_expr - predicted_batch_effects + reference_mean

            linear_corrected.X = X_corrected

            # Recalculate PCA and UMAP
            print("Recalculating PCA and UMAP...")
            sc.pp.scale(linear_corrected, max_value=10)
            sc.tl.pca(linear_corrected, n_comps=50)
            sc.pp.neighbors(linear_corrected, n_neighbors=10, n_pcs=20)
            sc.tl.umap(linear_corrected, random_state=0, n_components=2, maxiter=100)

            elapsed_time = time.time() - start_time
            print(f"Simplified batch correction completed, time taken: {elapsed_time:.2f} seconds")
            datasets["linear"] = (linear_corrected, None)

        except Exception as e:
            print(f"Simplified batch correction failed: {e}")
            import traceback
            traceback.print_exc()
            linear_corrected = None
    else:
        print("Skipping simplified batch correction (no batch information)")

    # 3.5 Harmony Batch Correction
    print("\n3.5 Harmony Batch Correction")

    harmony_data = None
    if HARMONY_AVAILABLE and batch_key:
        # Create Harmony data copy
        harmony_data = adata_hvg.copy()

        print("Starting Harmony batch correction...")
        start_time = time.time()

        try:
            # IMPORTANT: Harmony requires batch key to be string or categorical, not int
            # Convert batch key to string to avoid 'unique' error
            harmony_data.obs[batch_key] = harmony_data.obs[batch_key].astype(str)
            print(f"Converted {batch_key} to string type for Harmony compatibility")
            
            # Apply Harmony batch correction
            harmony_integrate(harmony_data, key=batch_key, max_iter_harmony=20)

            # Calculate UMAP
            sc.pp.neighbors(harmony_data, n_neighbors=10, n_pcs=20)
            sc.tl.umap(harmony_data, random_state=0, n_components=2, maxiter=100)

            elapsed_time = time.time() - start_time
            print(f"Harmony batch correction completed, time taken: {elapsed_time:.2f} seconds")
            datasets["harmony"] = (harmony_data, None)

        except Exception as e:
            print(f"Harmony batch correction failed: {e}")
            import traceback
            traceback.print_exc()
            harmony_data = None

    else:
        print("Skipping Harmony batch correction")

    # 3.6 Execute Clustering Analysis
    print("\n3.6 Execute Clustering Analysis")

    for suffix, (dataset, rep) in datasets.items():
        print(f"Executing clustering on {suffix or 'baseline'} dataset...")
        add_clusters(dataset, use_rep=rep, suffix=suffix)

    # 3.7 Calculate clustering evaluation metrics
    print("\n3.7 Calculate clustering evaluation metrics")

    clustering_results = []
    for suffix, (dataset, rep) in datasets.items():
        dataset_name = suffix or "baseline"
        embedding_key = rep or "X_pca"

        for method in ["kmeans"]:
            cluster_key = f"{method}{suffix}"
            if cluster_key in dataset.obs:
                metrics = compute_metrics(dataset, cluster_key, embedding_key)
                clustering_results.append({
                    "dataset": dataset_name,
                    "method": method,
                    "silhouette": metrics["silhouette"],
                    "ari": metrics["ari"],
                    "rand": metrics["rand"],
                    "n_clusters": dataset.obs[cluster_key].nunique()
                })

    clustering_df = pd.DataFrame(clustering_results)
    print("Clustering evaluation results:")
    print(clustering_df)

    # 3.8 UMAP Visualization (Reproduce Figure 1B)
    print("\n3.8 UMAP Visualization")

    # Calculate UMAP for baseline dataset
    baseline_dataset = adata_hvg
    sc.pp.neighbors(baseline_dataset, n_neighbors=10, n_pcs=20)
    sc.tl.umap(baseline_dataset, random_state=0, n_components=2, maxiter=100)

    # Create clustering comparison visualization (only show baseline)
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    sc.pl.umap(baseline_dataset, color="kmeans", show=False, ax=ax,
              legend_loc="on data", size=10, title="K-means Clustering (Baseline)")

    plt.tight_layout()
    plt.savefig("figures/umap_clustering_comparison.png", dpi=150, bbox_inches='tight')
    plt.close()

    print("UMAP clustering comparison plot saved")

    # 3.8.1 Create comprehensive UMAP clustering metrics comparison
    print("\n3.8.1 Create comprehensive UMAP clustering metrics comparison")

    if len(datasets) > 1:
        # Calculate UMAP for all datasets that don't have it
        for suffix, (dataset, rep) in datasets.items():
            if 'X_umap' not in dataset.obsm:
                try:
                    sc.pp.neighbors(dataset, n_neighbors=10, n_pcs=20, use_rep=rep)
                    sc.tl.umap(dataset, random_state=0, n_components=2, maxiter=100)
                    print(f"Calculated UMAP for {suffix or 'baseline'} dataset")
                except Exception as e:
                    print(f"Failed to calculate UMAP for {suffix or 'baseline'} dataset: {e}")

        # Create a figure comparing different batch correction methods
        n_methods = len(datasets)
        fig, axes = plt.subplots(1, n_methods, figsize=(6*n_methods, 5))

        if n_methods == 1:
            axes = [axes]

        for i, (suffix, (dataset, rep)) in enumerate(datasets.items()):
            method_name = suffix or "baseline"
            cluster_key = f"kmeans{suffix}"

            if cluster_key in dataset.obs and 'X_umap' in dataset.obsm:
                sc.pl.umap(dataset, color=cluster_key, show=False, ax=axes[i],
                          legend_loc="on data" if i == 0 else None,
                          size=8, title=f"{method_name.upper()} Clustering")

        plt.tight_layout()
        plt.savefig("figures/umap_clustering_metrics_comparison.png", dpi=150, bbox_inches='tight')
        plt.close()
        print("UMAP clustering metrics comparison plot saved")
    else:
        print("Skipping clustering metrics comparison (only one dataset available)")

    # 3.8.2 Create UMAP with time point coloring (Figure 1B replication)
    print("\n3.8.2 Create UMAP with time point coloring (Figure 1B replication)")

    # Use the batch key (DaysPostAmputation) as time information
    time_column = batch_key or 'time'
    if time_column and time_column in baseline_dataset.obs.columns:
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))

        # Color by time points
        sc.pl.umap(baseline_dataset, color=time_column, show=False, ax=ax,
                  size=10, title="UMAP by Time Points (Figure 1B Replication)")

        plt.tight_layout()
        plt.savefig("figures/umap_umap_figure_1b_replication.png", dpi=150, bbox_inches='tight')
        plt.close()
        print("UMAP Figure 1B replication plot saved")
    else:
        print("Skipping Figure 1B replication (no time information)")

    # 3.8.3 Create UMAP showing marker gene expression
    print("\n3.8.3 Create UMAP showing marker gene expression")

    # Identify top marker genes for the largest cluster
    if "kmeans" in baseline_dataset.obs:
        # Find the largest cluster
        largest_cluster = baseline_dataset.obs["kmeans"].value_counts().index[0]

        # Get marker genes for this cluster
        try:
            sc.tl.rank_genes_groups(baseline_dataset, "kmeans", method="wilcoxon",
                                   key_added="marker_genes")

            # Get top 6 marker genes
            top_markers = []
            for i in range(min(6, len(baseline_dataset.uns["marker_genes"]["names"][largest_cluster]))):
                gene = baseline_dataset.uns["marker_genes"]["names"][largest_cluster][i]
                if gene in baseline_dataset.var_names:
                    top_markers.append(gene)

            if top_markers:
                # Create UMAP plots for top marker genes
                n_markers = len(top_markers)
                n_cols = min(3, n_markers)
                n_rows = (n_markers + n_cols - 1) // n_cols

                fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))

                if n_rows == 1 and n_cols == 1:
                    axes = [axes]
                elif n_rows == 1:
                    axes = axes.flatten()
                else:
                    axes = axes.flatten()

                for i, gene in enumerate(top_markers):
                    if i < len(axes):
                        sc.pl.umap(baseline_dataset, color=gene, show=False, ax=axes[i],
                                  size=8, title=f"{gene} Expression", cmap="viridis")

                # Hide unused subplots
                for i in range(len(top_markers), len(axes)):
                    axes[i].set_visible(False)

                plt.tight_layout()
                plt.savefig("figures/umap_marker_gene_expression.png", dpi=150, bbox_inches='tight')
                plt.close()
                print(f"UMAP marker gene expression plot saved ({len(top_markers)} genes)")
            else:
                print("No valid marker genes found for expression visualization")

        except Exception as e:
            print(f"Failed to create marker gene expression plot: {e}")
    else:
        print("Skipping marker gene expression plot (no clustering available)")

    # 3.9 Save Results
    print("\n3.9 Save Results")

    # Save clustering metrics
    clustering_df.to_csv("results/clustering_metrics.csv", index=False)

    # Save batch-corrected data
    for suffix, (dataset, rep) in datasets.items():
        if suffix:  # Only save batch-corrected data
            dataset.write_h5ad(f"results/{suffix}_corrected.h5ad")

    # Save baseline dataset with clustering results
    baseline_dataset.write_h5ad("results/baseline_clustered.h5ad")

    print("All completed")

if __name__ == "__main__":
    main()