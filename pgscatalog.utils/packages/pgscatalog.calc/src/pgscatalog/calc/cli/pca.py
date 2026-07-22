import polars as pl
import numpy as np
from scipy.sparse import csr_array
from dask import array as da


# Generate per-chromsome LD matricies (with MAF information) on QC'd SNPs from reference panel

# 1. PER CHROMOSOME SELECT variants to be included in PCA:
# Choose variants that are matched (and have MAF information in the reference panel)
ld_pairs = pl.read_csv('./test_ancestry/reference/ld_ref_21.vcor.zst',
                       separator='\t',
                       new_columns=["ID_A", "ID_B",
                                    "UNPHASED_R2"])  # ToDo make the LD reading lazy so that it doesn't read variants that will just be filtered
n_vars = 10000
matched_var_maf = pl.DataFrame(
    {"ID": ld_pairs['ID_B'].unique().sample(n_vars), "MAF": np.random.uniform(0, 0.5, size=n_vars)}).with_row_index(
    name="index")

# Filter LD matrix to matched variants
ld_pairs_filt = ld_pairs.filter(pl.col("ID_A").is_in(matched_var_maf['ID']),
                                pl.col("ID_B").is_in(matched_var_maf['ID']))
ld_pairs_filt = (
    ld_pairs_filt
    # Join on ID_A to get the index for variant A
    .join(
        matched_var_maf.select(["index", "ID"]),
        left_on="ID_A",
        right_on="ID",
        how="left"
    )
    .rename({"index": "idx_A"})
    # Join on ID_B to get the index for variant B
    .join(
        matched_var_maf.select(["index", "ID"]),
        left_on="ID_B",
        right_on="ID",
        how="left"
    )
    .rename({"index": "idx_B"})
    # Select and reorder columns as needed
    .select(["idx_A", "idx_B", "UNPHASED_R2"])
)

# Convert to sparse LD matrix (symmetric)
ld_matrix = csr_array((ld_pairs_filt["UNPHASED_R2"], (ld_pairs_filt["idx_A"], ld_pairs_filt["idx_B"])),
                      shape=(n_vars, n_vars))
ld_matrix = ld_matrix + ld_matrix.T


# LD thin the matrix
def ld_thinning_maf_ordered(r2_matrix: csr_array, mafs: np.ndarray, r2_threshold: float = 0.2) -> np.ndarray:
    """Function to thin LD matrix, prioritisting higher MAF variants"""
    n_variants = r2_matrix.shape[0]
    keep_mask = np.ones(n_variants, dtype=bool)

    maf_order = np.argsort(mafs)[::-1]

    indptr = r2_matrix.indptr
    indices = r2_matrix.indices
    data = r2_matrix.data

    for i in maf_order:
        if not keep_mask[i]:
            continue

        start = indptr[i]
        end = indptr[i + 1]

        cols = indices[start:end]
        vals = data[start:end]

        prune_cols = cols[(vals > r2_threshold) & (cols != i)]  # in case symmetric

        if prune_cols.size > 0:
            keep_mask[prune_cols] = False

    return keep_mask


var_mask = ld_thinning_maf_ordered(ld_matrix, mafs=matched_var_maf['MAF'], r2_threshold=0.2)
thinned_matrix = ld_matrix[var_mask, :][:, var_mask]
selected_var_maf = matched_var_maf.filter(var_mask)

# 2. Assemble dosage array for the REFERENCE dataset and TARGET

# Generate dummy PCA data
def DummyDosgaes_TwoPop(n_variants=10000, fraction_mafdiff=0.05,
                        n_ref_samples=(3000,2000), n_target_samples=(0,1000), maf_range=(0.01, 0.5), missing_rate=0.01,
                        seed=42):
    """Simulate genotype dosages for a reference panel and a target population, with missingness"""
    # Random MAF for each variant
    rng = np.random.default_rng(seed)
    mafs = rng.uniform(0.01, 0.5, size=n_variants)
    n_maffdiff = int(n_variants*fraction_mafdiff)
    mafs_pop2 = mafs.copy()
    mafs_pop2[:n_maffdiff] = rng.uniform(maf_range[0], maf_range[1], size=n_maffdiff)

    # Simulate reference population:
    pop_ref = np.concatenate([np.zeros(n_ref_samples[0]), np.ones(n_ref_samples[1])]) + 1
    genotypes_ref = np.concatenate([
        rng.binomial(n=2, p=mafs[:, None], size=(n_variants, n_ref_samples[0])).astype(np.float32),
        rng.binomial(n=2, p=mafs_pop2[:, None], size=(n_variants, n_ref_samples[1])).astype(np.float32)],
        axis=1)
    missing_mask = np.random.rand(*genotypes_ref.shape) < missing_rate
    genotypes_ref[missing_mask] = np.nan

    pop_target = np.concatenate([np.zeros(n_target_samples[0]), np.ones(n_target_samples[1])]) + 1
    genotypes_target = np.concatenate([
        rng.binomial(n=2, p=mafs[:, None], size=(n_variants, n_target_samples[0])).astype(np.float32),
        rng.binomial(n=2, p=mafs_pop2[:, None], size=(n_variants, n_target_samples[1])).astype(np.float32)],
        axis=1)
    missing_mask = np.random.rand(*genotypes_target.shape) < missing_rate
    genotypes_target[missing_mask] = np.nan

    return (pop_ref, genotypes_ref), (pop_target, genotypes_target)

dd_ref, dd_target = DummyDosgaes_TwoPop()

def standardize(X, mean=None, std=None):
    assert np.issubdtype(X.dtype, np.floating)
    p, n = X.shape
    is_miss = np.isnan(X)
    if (mean is None) or (std is None):
        mean = np.zeros(p)
        std = np.zeros(p)
        for i in range(p):
            row_nomiss = X[i,:][~is_miss[i,:]]
            mean[i] = np.mean(row_nomiss)
            std[i] = np.std(row_nomiss)
        std[std == 0] = 1
    mean = mean.reshape((-1, 1))
    std = std.reshape((-1, 1))
    X_std = (X - mean) / std
    X_std[is_miss] = 0
    return X_std, mean, std


def svd_eigcov(XTX):
    ssq, V = np.linalg.eigh(XTX)
    V = np.squeeze(V)
    ssq = np.squeeze(ssq)
    s = np.sqrt(abs(ssq))
    V = V.T[::-1].T
    s = s[::-1]
    return s, V

def eig_ref(X):
    print('Calculating reference covariance matrix...')
    XTX = X.T @ X
    print('Eigendecomposition on reference covariance matrix...')
    s, V = svd_eigcov(XTX)
    return s, V, XTX
def source_fraposa_pca(genotypes_ref: np.array, genotypes_target= None, dim_ref=10, method='oadp', dim_stu=None, dim_online=None):
    # Setup the number of PCs we want to compute for difference methods
    dim_stu = dim_ref * 2
    assert method in ['oadp', 'adp', 'sp']
    if method in ['oadp', 'adp']:
        if dim_stu is None:
            dim_stu = dim_ref * 2
        assert dim_ref <= dim_stu
    if method in ['oadp']:
        if dim_online is None:
            dim_online = dim_stu * 2
        assert dim_stu <= dim_online

    # PCA: reference data
    X, X_mean, X_std = standardize(genotypes_ref.astype(np.float64))
    if method == 'oadp':
        s, V = eig_ref(X)[:2]
        V = V[:, :dim_online]
        U = X @ (V / s[:dim_online])
    elif method == 'sp':
        s, V = eig_ref(X)[:2]
        V = V[:, :dim_ref]
        U = X @ (V / s[:dim_ref])
    elif method == 'adp':
        s, V, XTX = eig_ref(X)
        V = V[:, :dim_ref]
    pcs_ref = V[:, :dim_ref] * s[:dim_ref]

    return pcs_ref

pca_fraposa = source_fraposa_pca(genotypes_ref=dd_ref[1])
np.savetxt('test_pca.txt', r_pca)


def standardize_dask(X, mean=None, std=None):
    """Vectorized standardization for Dask arrays."""
    if (mean is None) or (std is None):
        # Compute per-row (per-variant) stats along samples (axis 1)
        mean = da.nanmean(X, axis=1, keepdims=True)
        std = da.nanstd(X, axis=1, keepdims=True, ddof=0)

        # Handle monomorphic variants (std == 0)
        std = da.where(std == 0, 1.0, std)
    else:
        # Ensure 2D column shapes for broadcasting (p, 1)
        mean = da.asarray(mean).reshape((-1, 1))
        std = da.asarray(std).reshape((-1, 1))

    # Standardize: (p, n) - (p, 1) / (p, 1)
    X_std = (X - mean) / std

    # Impute missing values to 0
    X_std = da.where(da.isnan(X), 0.0, X_std)

    return X_std, mean, std


def svd_eigcov_dask(XTX):
    """SVD implementation aligned specifically to match `np.linalg.eigh` outputs."""
    # 1. Compute SVD on the symmetric covariance matrix XTX
    # SVD on XTX yields: ssq (eigenvalues) and V (right singular vectors)
    _, ssq, VT = da.linalg.svd(XTX)

    # 2. V is VT.T
    V = VT.T

    # 3. Ensure ssq is positive (numerical precision safety)
    ssq = da.abs(ssq)

    # 4. Standard deviations of principal components s = sqrt(eigenvalues)
    # Note: da.linalg.svd ALREADY outputs values sorted descending (largest first)
    s = da.sqrt(ssq)

    # No array reversing needed here because SVD is already largest-to-smallest!
    return s, V


def eig_ref_dask(X):
    print("Calculating reference covariance matrix (n x n)...")
    # X is (p x n) -> X.T @ X is (n x n)
    XTX = X.T @ X

    print("Eigendecomposition on covariance matrix...")
    s, V = svd_eigcov_dask(XTX)

    return s, V, XTX


def source_fraposa_pca_dask(
    genotypes_ref,
    genotypes_target=None,
    dim_ref=10,
    method="oadp",
    dim_stu=None,
    dim_online=None,
):
    # Setup dimension requirements
    dim_stu = dim_ref * 2
    assert method in ["oadp", "adp", "sp"]

    if method in ["oadp", "adp"]:
        if dim_stu is None:
            dim_stu = dim_ref * 2
        assert dim_ref <= dim_stu
    if method == "oadp":
        if dim_online is None:
            dim_online = dim_stu * 2
        assert dim_stu <= dim_online

    # Convert genotype matrix to Dask array if it isn't one already
    genotypes_ref = da.asarray(genotypes_ref, dtype=float)

    # 1. Standardize Reference
    X, X_mean, X_std = standardize_dask(genotypes_ref)

    # 2. Compute Eigendecomposition
    if method == "oadp":
        s, V, _ = eig_ref_dask(X)
        V = V[:, :dim_online]
        U = X @ (V / s[:dim_online])
    elif method == "sp":
        s, V, _ = eig_ref_dask(X)
        V = V[:, :dim_ref]
        U = X @ (V / s[:dim_ref])
    elif method == "adp":
        s, V, XTX = eig_ref_dask(X)
        V = V[:, :dim_ref]

    # Reference PCs
    pcs_ref = V[:, :dim_ref] * s[:dim_ref]

    return pcs_ref

pca_dask = source_fraposa_pca_dask(genotypes_ref=dd_ref[1]).compute()
