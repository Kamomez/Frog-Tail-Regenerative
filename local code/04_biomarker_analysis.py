from pathlib import Path
from sklearn.linear_model import LogisticRegression
from scipy.stats import ranksums
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os
import seaborn as sns

def find_markers_logistic_regression(adata, cluster_key, cluster_of_interest):
    # Identify marker genes using logistic regression
    # Prepare data
    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    y = (adata.obs[cluster_key] == cluster_of_interest).astype(int)

    # Train logistic regression model
    clf = LogisticRegression(random_state=42, max_iter=1000)
    clf.fit(X, y)

    # Get feature importance
    feature_importance = np.abs(clf.coef_[0])

    # Create results DataFrame
    markers_df = pd.DataFrame({
        'gene': adata.var_names,
        'importance': feature_importance
    })

    # Sort by importance
    markers_df = markers_df.sort_values('importance', ascending=False)

    return markers_df

def find_markers_wilcoxon(adata, cluster_key, cluster_of_interest):
    # Identify marker genes using Wilcoxon rank-sum test
    cluster_mask = adata.obs[cluster_key] == cluster_of_interest
    other_mask = adata.obs[cluster_key] != cluster_of_interest

    p_values = []
    log_fold_changes = []

    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X

    for i in range(X.shape[1]):
        cluster_expr = X[cluster_mask, i]
        other_expr = X[other_mask, i]

        # Calculate log fold change
        cluster_mean = np.mean(cluster_expr)
        other_mean = np.mean(other_expr)
        lfc = np.log2(cluster_mean + 1) - np.log2(other_mean + 1)
        log_fold_changes.append(lfc)

        # Wilcoxon rank-sum test
        try:
            stat, p_val = ranksums(cluster_expr, other_expr)
            p_values.append(p_val)
        except ValueError:
            p_values.append(1.0)

    # Create results DataFrame
    markers_df = pd.DataFrame({
        'gene': adata.var_names,
        'p_value': p_values,
        'log_fold_change': log_fold_changes,
        'neg_log_p': -np.log10(np.array(p_values) + 1e-10)
    })

    # Sort by p-value
    markers_df = markers_df.sort_values('p_value', ascending=True)

    return markers_df

def main():
    print("Biomarker Identification and Functional Analysis")

    # Create necessary directories
    Path("results").mkdir(exist_ok=True)
    Path("figures").mkdir(exist_ok=True)

    # Check if clustering data exists
    baseline_file = Path("results/baseline_denoised.h5ad")
    if not baseline_file.exists():
        print("Error: Clustering data file not found, please run 03_batch_correction_clustering.py first")
        return

    # Load data (use clustered baseline data)
    print("Loading data...")
    baseline_file = "results/baseline_clustered.h5ad"
    if not os.path.exists(baseline_file):
        print(f"Warning: {baseline_file} not found, using preprocessed data")
        baseline_file = "results/preprocessed_data.h5ad"

    adata_hvg = sc.read_h5ad(baseline_file)

    # Check batch information
    batch_key = 'DaysPostAmputation' if 'DaysPostAmputation' in adata_hvg.obs.columns else None

    # 4.1 Identify potential ROC clusters (with leiden+kmeans support)
    print("\n4.1 Identify potential ROC clusters")

    potential_roc_clusters = {}

    # Support both leiden and kmeans - automatically select the best one
    for method in ["leiden", "kmeans"]:
        cluster_key = f"{method}"
        if cluster_key in adata_hvg.obs:
            print(f"Analyzing {method.upper()} clustering results:")

            if batch_key and batch_key in adata_hvg.obs.columns:
                cluster_time_dist = pd.crosstab(adata_hvg.obs[cluster_key], adata_hvg.obs[batch_key], normalize='index')
                early_timepoints = cluster_time_dist.columns[:2]
                early_enrichment = cluster_time_dist[early_timepoints].sum(axis=1)

                roc_candidates = early_enrichment[early_enrichment > early_enrichment.mean() + early_enrichment.std()].index.tolist()

                if roc_candidates:
                    potential_roc_clusters[f"baseline_{method}"] = roc_candidates
                    print(f"Potential ROC clusters: {roc_candidates}")
                    print(f"Early time point enrichment: {early_enrichment[roc_candidates].values}")
                else:
                    print("No significantly enriched clusters found")

    print(f"Potential ROC clusters summary: {potential_roc_clusters}")

    # Select the best ROC cluster (highest early enrichment)
    best_roc_key = None
    best_enrichment = 0
    best_roc_cluster = None

    if potential_roc_clusters:
        for roc_key, clusters in potential_roc_clusters.items():
            dataset_name, method = roc_key.rsplit('_', 1)
            cluster_key = f"{method}"
            
            if cluster_key in adata_hvg.obs and batch_key in adata_hvg.obs.columns:
                cluster_time_dist = pd.crosstab(adata_hvg.obs[cluster_key], adata_hvg.obs[batch_key], normalize='index')
                early_timepoints = cluster_time_dist.columns[:2]
                early_enrichment = cluster_time_dist[early_timepoints].sum(axis=1)
                
                for cluster in clusters:
                    enrichment = early_enrichment.loc[cluster]
                    if enrichment > best_enrichment:
                        best_enrichment = enrichment
                        best_roc_key = roc_key
                        best_roc_cluster = cluster
        
        if best_roc_key:
            print(f"   Best ROC cluster selected: {best_roc_cluster} from {best_roc_key}")
            print(f"   Early enrichment: {best_enrichment:.1%}")

    # 4.2 Marker gene identification
    print("\n4.2 Marker gene identification")

    cluster_key = "kmeans"  # Default to kmeans clustering

    if best_roc_cluster:
        dataset_name, method = best_roc_key.rsplit('_', 1)
        cluster_key = f"{method}"

        if cluster_key in adata_hvg.obs:
            roc_cluster = best_roc_cluster
            print(f"Selected ROC cluster: {roc_cluster} (from {best_roc_key})")

            # Logistic regression
            print("Using Logistic regression to identify marker genes...")
            markers_logistic = find_markers_logistic_regression(adata_hvg, cluster_key, roc_cluster)
            markers_logistic.to_csv("results/roc_markers_logistic.csv", index=False)

            # Wilcoxon test
            print("Using Wilcoxon test to identify marker genes...")
            markers_wilcoxon = find_markers_wilcoxon(adata_hvg, cluster_key, roc_cluster)
            markers_wilcoxon.to_csv("results/roc_markers_wilcoxon.csv", index=False)

            print("Logistic regression - top 10 marker genes:")
            print(markers_logistic.head(10))

            print("Wilcoxon test - top 10 marker genes:")
            print(markers_wilcoxon.head(10))
    else:
        print("No ROC clusters identified, proceeding with general marker gene analysis")

    # 4.2.1 Marker gene visualization - Always run for kmeans
    print("\n4.2.1 Marker gene visualization")

    cluster_key = "kmeans"  # Use kmeans clustering for visualization
    if cluster_key in adata_hvg.obs:
        # Compute dendrogram first to remove warning
        try:
            sc.tl.dendrogram(adata_hvg, groupby=cluster_key)
        except:
            pass  # Skip if dendrogram computation fails

        # Compute marker genes using Wilcoxon test (more robust for sparse data)
        print("Computing marker genes for visualization...")
        sc.tl.rank_genes_groups(
            adata_hvg,
            groupby=cluster_key,
            method="wilcoxon",
            n_genes=30,
            key_added="rank_genes_clusters"
        )

        # Save marker genes results
        marker_df = sc.get.rank_genes_groups_df(adata_hvg, group=None, key="rank_genes_clusters")
        marker_df.to_csv("results/marker_genes_all_clusters.csv", index=False)
        print("Marker genes saved to: results/marker_genes_all_clusters.csv")

        # Visualization - rank genes groups plot
        print("Creating rank genes groups plot...")
        try:
            sc.pl.rank_genes_groups(adata_hvg, key="rank_genes_clusters", n_genes=5, sharey=False, show=False)
            plt.savefig('figures/rank_genes_groups_plot.png', dpi=300, bbox_inches='tight')
            plt.close()
            print("Rank genes groups plot saved to figures/rank_genes_groups_plot.png")
        except Exception as e:
            print(f"Rank genes groups plot failed: {e}")

        # Visualization - dotplot (vertical orientation)
        print("Creating dotplot...")
        try:
            sc.pl.rank_genes_groups_dotplot(adata_hvg, key="rank_genes_clusters", n_genes=5, swap_axes=True, show=False)
            plt.savefig('figures/marker_genes_dotplot.png', dpi=300, bbox_inches='tight')
            plt.close()
            print("Dotplot saved to figures/marker_genes_dotplot.png")
        except Exception as e:
            print(f"Dotplot failed: {e}")

        # Visualization - heatmap
        print("Creating heatmap...")
        try:
            sc.pl.rank_genes_groups_heatmap(adata_hvg, key="rank_genes_clusters", n_genes=10,
                                          show_gene_labels=True, show=False, swap_axes=True, figsize=(8, 20))
            plt.savefig('figures/marker_genes_heatmap.png', dpi=300, bbox_inches='tight')
            plt.close()
            print("Heatmap saved to figures/marker_genes_heatmap.png")
        except Exception as e:
            print(f"Heatmap failed: {e}")

        print("Marker gene visualization completed")

    # 4.3 Comparison with Supplementary Table 3
    print("\n4.3 Comparison with Supplementary Table 3")

    known_roc_markers = [
        "wnt11b", "fzd10", "notum", "sfrp1", "sfrp5", "dkk1", "axin2",
        "lef1", "tcf7l1", "myc", "ccnd1", "cdkn1a", "cdkn1b", "tp53"
    ]

    print("Known ROC marker genes:")
    print(known_roc_markers)

    if 'markers_logistic' in locals():
        top_markers = markers_logistic.head(50)['gene'].tolist()
        overlap = set(top_markers) & set(known_roc_markers)
        print(f"Overlap with known ROC marker genes: {len(overlap)}")
        print(f"Overlapping genes: {list(overlap)}")

        overlap_df = pd.DataFrame({
            'method': ['logistic_regression'] * len(overlap),
            'overlapping_gene': list(overlap)
        })
        overlap_df.to_csv("results/marker_overlap_summary.csv", index=False)

    # 4.3.1 Cross-method marker gene validation
    print("\n4.3.1 Cross-method marker gene validation")

    if 'markers_logistic' in locals() and 'markers_wilcoxon' in locals():
        # Get top genes from each method
        top_logistic = set(markers_logistic.head(50)['gene'].tolist())
        top_wilcoxon = set(markers_wilcoxon.head(50)['gene'].tolist())
        
        # Find overlap
        common_markers = top_logistic & top_wilcoxon
        
        print(f"\nMarker gene comparison:")
        print(f"  Logistic regression (top 50): {len(top_logistic)} genes")
        print(f"  Wilcoxon test (top 50): {len(top_wilcoxon)} genes")
        print(f"  Common markers: {len(common_markers)} genes")
        print(f"  Overlap percentage: {len(common_markers)/50*100:.1f}%")
        
        if common_markers:
            print(f"\nHigh-confidence markers (validated by both methods):")
            # Get details for common markers
            common_markers_list = sorted(list(common_markers))
            for i, gene in enumerate(common_markers_list[:10], 1):  # Show top 10
                logistic_rank = markers_logistic[markers_logistic['gene'] == gene].index[0] + 1
                wilcoxon_rank = markers_wilcoxon[markers_wilcoxon['gene'] == gene].index[0] + 1
                logistic_score = markers_logistic[markers_logistic['gene'] == gene]['importance'].values[0]
                print(f"  {i}. {gene}")
                print(f"     - Logistic: rank #{logistic_rank}, score={logistic_score:.3f}")
                print(f"     - Wilcoxon: rank #{wilcoxon_rank}")
            
            # Save common markers
            common_markers_df = pd.DataFrame({
                'gene': list(common_markers),
                'method': 'both'
            })
            common_markers_df.to_csv("results/common_markers_validated.csv", index=False)
            print(f"\nCommon markers saved to: results/common_markers_validated.csv")
            
            # Create bar chart comparing methods
            fig, ax = plt.subplots(figsize=(10, 6))
            
            comparison_data = pd.DataFrame({
                'Method': ['Logistic\nOnly', 'Both\nMethods', 'Wilcoxon\nOnly'],
                'Count': [
                    len(top_logistic - common_markers),
                    len(common_markers),
                    len(top_wilcoxon - common_markers)
                ],
                'Color': ['#377eb8', '#4daf4a', '#e41a1c']
            })
            
            bars = ax.bar(comparison_data['Method'], comparison_data['Count'], 
                         color=comparison_data['Color'], edgecolor='black', linewidth=1.5, alpha=0.8)
            
            ax.set_ylabel('Number of Marker Genes', fontsize=12, fontweight='bold')
            ax.set_title('Marker Gene Identification - Method Comparison', fontsize=14, fontweight='bold', pad=20)
            ax.grid(True, alpha=0.3, axis='y')
            
            # Add value labels on bars
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{int(height)}',
                       ha='center', va='bottom', fontweight='bold', fontsize=12)
            
            plt.tight_layout()
            plt.savefig('figures/marker_method_comparison.png', dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()
            
            print("Method comparison plot saved to: figures/marker_method_comparison.png")
        else:
            print("No common markers found between methods")
    else:
        print("Marker genes not yet identified, skipping validation")

    # 4.4 ROC marker genes heatmap visualization
    print("\n4.4 ROC marker genes heatmap visualization")

    if 'markers_logistic' in locals() and cluster_key in adata_hvg.obs:
        # Select top 20 ROC marker genes
        top_roc_markers = markers_logistic.head(20)['gene'].tolist()

        # Ensure genes exist in data
        available_markers = [gene for gene in top_roc_markers if gene in adata_hvg.var_names]
        print(f"Using {len(available_markers)} available marker genes for heatmap")

        if available_markers:
            # Calculate average expression of marker genes in each cluster
            cluster_means = []
            cluster_labels = []

            for cluster in sorted(adata_hvg.obs[cluster_key].unique()):
                cluster_mask = adata_hvg.obs[cluster_key] == cluster
                cluster_data = adata_hvg[cluster_mask]

                # Calculate average expression of marker genes in this cluster
                gene_means = []
                for gene in available_markers:
                    gene_expr = cluster_data[:, gene].X
                    if hasattr(gene_expr, 'toarray'):
                        gene_expr = gene_expr.toarray()
                    mean_expr = np.mean(gene_expr) if gene_expr.size > 0 else 0
                    gene_means.append(mean_expr)

                cluster_means.append(gene_means)
                cluster_labels.append(f"Cluster {cluster}")

            # Convert to numpy array and transpose to make it vertical
            heatmap_data = np.array(cluster_means).T

            # Create heatmap
            plt.figure(figsize=(10, 16))

            # Create heatmap using viridis colormap
            sns.heatmap(heatmap_data,
                       xticklabels=cluster_labels,
                       yticklabels=available_markers,
                       cmap='viridis',
                       annot=False,
                       cbar_kws={'label': 'Mean Expression'})

            plt.title('ROC Marker Genes Expression Across Clusters', fontsize=14, pad=20)
            plt.xlabel('Clusters', fontsize=12)
            plt.ylabel('ROC Marker Genes', fontsize=12)

            # Rotate labels for better readability
            plt.xticks(rotation=45, ha='right')
            plt.yticks(rotation=0)

            plt.tight_layout()
            plt.savefig('figures/roc_marker_heatmap.png', dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()

            print("ROC marker genes heatmap saved to figures/roc_marker_heatmap.png")
        else:
            print("Not available")
    else:
        print("Need to identify ROC marker genes and clustering results first")

    # 4.5 GO enrichment analysis
    print("\n4.5 GO enrichment analysis")

    try:
        import gseapy as gp

        if 'markers_logistic' in locals():
            roc_genes = markers_logistic.head(50)['gene'].tolist()

            print(f"Number of genes analyzed: {len(roc_genes)}")
            print(f"Gene list example: {roc_genes[:5]}")  # Show first 5 genes

            # Try with more lenient settings, but expect failure
            try:
                go_results = gp.enrichr(
                    gene_list=roc_genes,
                    gene_sets=['KEGG_2021_Human'],  # Use human KEGG pathways
                    organism='human',  # Use human genome as reference
                    outdir='results/go_analysis',
                    cutoff=0.5  # Further lower threshold
                )

                # Display results
                print("\nKEGG pathway analysis results:")
                if not go_results.results.empty:
                    print(go_results.results.head(10)[['Term', 'Adjusted P-value', 'Genes']])
                    # Save results
                    go_results.results.to_csv("results/go_enrichment_results.csv", index=False)
                else:
                    print("No significantly enriched pathways found")

            except Exception as e:
                print(f"KEGG analysis also failed: {e}")

        else:
            print("Need to identify marker genes first for GO analysis")

    except Exception as e:
        print(f"GO analysis error: {e}")

    # 4.6 Save analysis results
    print("\n4.6 Save analysis results")

    # Create analysis summary
    analysis_summary = {
        "analysis_type": ["roc_cluster_identification", "marker_gene_logistic", "marker_gene_wilcoxon", "go_enrichment", "heatmap_visualization"],
        "completed": [bool(potential_roc_clusters), 'markers_logistic' in locals(), 'markers_wilcoxon' in locals(), False, 'figures/roc_marker_heatmap.png' in [f.name for f in Path('figures').glob('*')]],
        "output_files": [
            "results/roc_markers_logistic.csv, results/roc_markers_wilcoxon.csv",
            "results/roc_markers_logistic.csv",
            "results/roc_markers_wilcoxon.csv",
            "results/go_enrichment_results.csv",
            "figures/roc_marker_heatmap.png"
        ]
    }

    summary_df = pd.DataFrame(analysis_summary)
    summary_df.to_csv("results/analysis_summary.csv", index=False)

    print("All completed!")
    print("Generated files:")
    print("- results/roc_markers_logistic.csv")
    print("- results/roc_markers_wilcoxon.csv")
    print("- results/marker_overlap_summary.csv")
    print("- results/go_enrichment_results.csv (if gseapy available)")
    print("- figures/roc_marker_heatmap.png")
    print("- results/analysis_summary.csv")

if __name__ == "__main__":
    main()