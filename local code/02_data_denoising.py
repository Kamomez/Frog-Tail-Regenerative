
#Methods: MAGIC, ALRA, DCA, kNN-smoothing

import time
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors

def main():
    print("=== Frog Tail Regeneration Analysis - 2. Data Denoising Techniques Comparison ===")

    # Create necessary directories
    Path("results").mkdir(exist_ok=True)
    Path("figures").mkdir(exist_ok=True)

    # Check if preprocessed data exists
    preprocessed_file = Path("results/preprocessed_data.h5ad")
    if not preprocessed_file.exists():
        print("Error: Preprocessed data file not found, please run 01_data_preprocessing.py first")
        return

    # Load preprocessed data
    print("Loading preprocessed data...")
    adata_hvg = sc.read_h5ad(preprocessed_file)
    
    print(f"Dataset info: {adata_hvg.n_obs} cells, {adata_hvg.n_vars} genes")

    # Import denoising libraries
    try:
        from scanpy.external.pp import magic as magic_denoise
        MAGIC_AVAILABLE = True
    except ImportError:
        print("Warning: MAGIC not available, please install magic-impute")
        MAGIC_AVAILABLE = False

    ALRA_AVAILABLE = True

    try:
        import dca
        from dca.api import dca as dca_denoise
        DCA_AVAILABLE = True
        print("DCA successfully imported")
    except (ImportError, ModuleNotFoundError) as e:
        print(f"Note: DCA not available (Keras compatibility issues)")
        DCA_AVAILABLE = False

    KNN_AVAILABLE = True

    print(f"MAGIC available: {MAGIC_AVAILABLE}")
    print(f"ALRA available: {ALRA_AVAILABLE} (using sklearn implementation)")
    print(f"DCA available: {DCA_AVAILABLE}")
    print(f"kNN-smoothing available: {KNN_AVAILABLE}")

    # 3.1 MAGIC Denoising
    magic_data = None
    if MAGIC_AVAILABLE:
        try:
            # Create MAGIC data copy
            magic_data = adata_hvg.copy()

            print("Starting MAGIC denoising...")
            start_time = time.time()

            # Apply MAGIC denoising
            magic_denoise(magic_data, name_list=list(magic_data.var_names), knn=17)

            # Re-scale and PCA
            sc.pp.scale(magic_data, max_value=10)
            sc.tl.pca(magic_data, n_comps=50)

            # Store MAGIC embedding
            magic_data.obsm["X_magic"] = magic_data.obsm.get("X_pca")

            elapsed_time = time.time() - start_time
            print(f"MAGIC denoising completed, time elapsed: {elapsed_time:.2f} seconds")
            print(f"MAGIC data shape: {magic_data.shape}")
        except Exception as e:
            print(f"MAGIC denoising failed: {e}")
            print("This is often due to numba/numpy compatibility issues")
            print("Continuing with other denoising methods...")
            magic_data = None
    else:
        print("Skipping MAGIC denoising")
        magic_data = None

    # 3.2 ALRA Denoising
    alra_data = None
    if ALRA_AVAILABLE:
        print("\n3.2 ALRA Denoising")
        print("Starting ALRA denoising (Adaptively-thresholded Low Rank Approximation)...")
       
        alra_data = adata_hvg.copy()

        print(f"ALRA data shape: {alra_data.shape}")
        print(f"Data type: {type(alra_data.X)}")
        print(f"Data range: {alra_data.X.min():.1f} - {alra_data.X.max():.1f}")

        start_time = time.time()

        try:
            print("Applying low-rank approximation with adaptive thresholding...")
            X = alra_data.X.toarray() if hasattr(alra_data.X, 'toarray') else alra_data.X

            n_components = min(100, min(X.shape) - 1)

            svd = TruncatedSVD(n_components=n_components, random_state=42)
            X_transformed = svd.fit_transform(X)
            
            # Calculate explained variance ratio
            explained_variance = svd.explained_variance_
            explained_variance_ratio = svd.explained_variance_ratio_
            cumulative_variance = np.cumsum(explained_variance_ratio)

            variance_threshold = 0.90
            optimal_components = np.where(cumulative_variance >= variance_threshold)[0]
            
            if len(optimal_components) > 0:
                optimal_k = optimal_components[0] + 1
                optimal_k = max(optimal_k, 15)  # ALRA recommends at least 10-20 components
                print(f"Adaptive rank selection: k={optimal_k} (explaining {cumulative_variance[optimal_k-1]:.1%} variance)")
            else:
                optimal_k = max(20, n_components // 2)
                print(f"Using default rank: k={optimal_k}")
            
            # Step 2: Recompute SVD with optimal rank
            svd_optimal = TruncatedSVD(n_components=optimal_k, random_state=42)
            X_transformed = svd_optimal.fit_transform(X)
            X_reconstructed = svd_optimal.inverse_transform(X_transformed)
            
            # Step 3: Adaptive thresholding (key ALRA step)
            print("Applying adaptive thresholding...")

            threshold = 0.1  # Threshold for considering a value as "zero" in log-space
            imputation_mask = (X < threshold)

            X_alra = X.copy()
            X_alra[imputation_mask] = X_reconstructed[imputation_mask]
            
            # Ensure non-negative values
            X_alra = np.maximum(X_alra, 0)
            
            # Calculate imputation statistics
            n_imputed = np.sum(imputation_mask)
            n_total = imputation_mask.size
            imputation_rate = n_imputed / n_total * 100
            print(f"Imputed {n_imputed:,} / {n_total:,} values ({imputation_rate:.1f}%)")
            
            # Update AnnData object
            alra_data.X = X_alra

            # Re-scale and PCA
            sc.pp.scale(alra_data, max_value=10)
            sc.tl.pca(alra_data, n_comps=50)

            # Store ALRA embedding
            alra_data.obsm["X_alra"] = alra_data.obsm.get("X_pca")

            elapsed_time = time.time() - start_time
            print(f"ALRA denoising completed, time elapsed: {elapsed_time:.2f} seconds")
            print(f"Final rank: k={optimal_k}, variance explained: {cumulative_variance[optimal_k-1]:.1%}")

        except Exception as e:
            print(f"ALRA denoising failed: {e}")
            import traceback
            traceback.print_exc()
            alra_data = None

    else:
        print("Skipping ALRA denoising")
        alra_data = None

    # 3.3 DCA
    dca_data = None
    if DCA_AVAILABLE:
        print("\n3.3 DCA Denoising")
        print("Starting DCA denoising...")
        
        # Create DCA data copy
        dca_data = adata_hvg.copy()
        
        # DCA works best with raw counts
        if dca_data.raw is not None:
            print("Using raw counts from .raw layer")
            dca_input = sc.AnnData(
                X=dca_data.raw.X,
                obs=dca_data.obs,
                var=dca_data.raw.var
            )
        else:
            print("Note: No raw counts available, using current data")
            dca_input = dca_data.copy()
        
        print(f"DCA input shape: {dca_input.shape}")
        start_time = time.time()
        
        try:
            # Apply DCA denoising
            print("Running DCA with ZINB autoencoder...")
            dca_denoise(
                dca_input,
                mode='denoise',
                ae_type='zinb-conddisp',
                normalize_per_cell=True,
                scale=True,
                log1p=True,
                hidden_size=[64, 32, 64],
                epochs=300,
                learning_rate=0.001,
                verbose=True
            )
            
            # Copy denoised data
            dca_data.X = dca_input.X
            
            # Compute PCA
            sc.pp.scale(dca_data, max_value=10)
            sc.tl.pca(dca_data, n_comps=50)
            dca_data.obsm["X_dca"] = dca_data.obsm.get("X_pca")
            
            elapsed_time = time.time() - start_time
            print(f"DCA denoising completed, time taken: {elapsed_time:.2f} seconds")
            print(f"DCA data shape: {dca_data.shape}")
            
        except Exception as e:
            print(f"DCA denoising failed: {e}")
            print("This is expected due to Keras compatibility issues")
            dca_data = None
    
    else:
        print("\nSkipping DCA denoising (library not available or incompatible)")
        dca_data = None

    # 3.4 kNN-smoothing Denoising
    knn_data = None
    if KNN_AVAILABLE:
        print("\n3.4 kNN-smoothing Denoising")
        print("Starting kNN-smoothing denoising...")
        
        # Create kNN data copy
        knn_data = adata_hvg.copy()
        
        print(f"kNN input shape: {knn_data.shape}")
        start_time = time.time()
        
        try:
            # Get data as array
            if hasattr(knn_data.X, 'toarray'):
                X = knn_data.X.toarray()
            else:
                X = knn_data.X.copy()
            
            # Set k-nearest neighbors
            n_neighbors = 15
            print(f"Using {n_neighbors} nearest neighbors for smoothing...")
            
            # Find k nearest neighbors
            nbrs = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean', n_jobs=-1)
            nbrs.fit(X)
            distances, indices = nbrs.kneighbors(X)
            
            # Weight by inverse distance
            epsilon = 1e-10  # Avoid division by zero
            weights = 1.0 / (distances + epsilon)
            weights = weights / weights.sum(axis=1, keepdims=True)
            
            # Apply weighted averaging
            print("Applying weighted k-nearest neighbor averaging...")
            X_smoothed = np.zeros_like(X)
            for i in range(len(X)):
                neighbor_expr = X[indices[i]]
                X_smoothed[i] = (neighbor_expr * weights[i][:, np.newaxis]).sum(axis=0)
            
            # Store smoothed data
            knn_data.X = X_smoothed
            
            # Compute PCA
            sc.pp.scale(knn_data, max_value=10)
            sc.tl.pca(knn_data, n_comps=50)
            knn_data.obsm["X_knn"] = knn_data.obsm.get("X_pca")
            
            elapsed_time = time.time() - start_time
            print(f"kNN-smoothing completed, time taken: {elapsed_time:.2f} seconds")
            print(f"kNN data shape: {knn_data.shape}")
            print(f"Smoothing strength: {n_neighbors} neighbors")
            
        except Exception as e:
            print(f"kNN-smoothing failed: {e}")
            import traceback
            traceback.print_exc()
            knn_data = None
    
    else:
        print("\nSkipping kNN-smoothing")
        knn_data = None

    # Save results
    print("\nSaving results...")
    adata_hvg.write_h5ad("results/baseline_denoised.h5ad")
    if magic_data is not None:
        magic_data.write_h5ad("results/magic_denoised.h5ad")
        print("Saved: results/magic_denoised.h5ad")
    if alra_data is not None:
        alra_data.write_h5ad("results/alra_denoised.h5ad")
        print("Saved: results/alra_denoised.h5ad")
    if dca_data is not None:
        dca_data.write_h5ad("results/dca_denoised.h5ad")
        print("Saved: results/dca_denoised.h5ad")
    if knn_data is not None:
        knn_data.write_h5ad("results/knn_denoised.h5ad")
        print("Saved: results/knn_denoised.h5ad")

    print("\nDenoising analysis completed!")
    total_methods = sum([magic_data is not None, alra_data is not None, dca_data is not None, knn_data is not None])
    print(f"Successfully denoised datasets: {total_methods} / 4")

if __name__ == "__main__":
    main()