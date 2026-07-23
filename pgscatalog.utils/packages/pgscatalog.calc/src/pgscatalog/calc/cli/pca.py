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
def DummyDosgaes_TwoPop(n_variants=10000, fraction_mafdiff=0.1,
                        n_ref_samples=(2000,1000), n_target_samples=(0,100), maf_range=(0.01, 0.5), missing_rate=0.01,
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

def svd_online(U1, d1, V1, b, l=None):
    n, k = V1.shape
    if l is None:
        l = k
    assert U1.shape[1] == k
    assert len(d1) == k
    assert(l <= k)
    p = U1.shape[0]
    b = b.reshape((p,1)) # Make sure the new sample is a column vec
    b_tilde = b - U1 @ (U1.T @ b)
    b_tilde = b_tilde / np.sqrt(sum(np.square(b_tilde)))
    R = np.concatenate((np.diag(d1), U1.transpose() @ b), axis = 1)
    R_tail = np.concatenate((np.zeros((1,k)), b_tilde.transpose() @ b), axis = 1)
    R = np.concatenate((R, R_tail), axis = 0)
    d2, R_Vt = np.linalg.svd(R, full_matrices=False)[1:]
    V_new = np.zeros((k+1, n+1))
    V_new[:k, :n] = V1.transpose()
    V_new[k, n] = 1
    V2 = (R_Vt @ V_new).transpose()[:,:l]
    return d2, V2

def procrustes(Y_mat, X_mat, return_transformed=False):
    ''' Find the best transformation from X to Y '''
    X = np.array(X_mat, dtype=np.double, copy=True)
    Y = np.array(Y_mat, dtype=np.double, copy=True)
    n_X = X.shape[0]
    X_mean = np.mean(X, 0)
    Y_mean = np.mean(Y, 0)
    X -= X_mean
    Y -= Y_mean
    C = Y.T @ X
    U, s, VT = np.linalg.svd(C, full_matrices=False)
    trXX = np.sum(X**2)
    trS = np.sum(s)
    R = VT.T @ U.T
    rho = trS / trXX
    c = Y_mean - rho * X_mean @ R
    if return_transformed:
        X_new = X_mat @ R * rho + c
        return R, rho, c, X_new
    else:
        return R, rho, c

def ref_aug_procrustes(pcs_ref, pcs_aug):
    n_ref, p_ref = pcs_ref.shape
    n_aug, p_aug = pcs_aug.shape
    assert n_aug == n_ref + 1
    assert p_aug >= p_ref
    pcs_aug_head = pcs_aug[:-1, :]
    pcs_aug_tail = pcs_aug[-1, :].reshape((1,-1))
    R, rho, c = procrustes_diffdim(pcs_ref, pcs_aug_head)
    pcs_aug_tail_trsfed = pcs_aug_tail @ R * rho + c
    return pcs_aug_tail_trsfed.flatten()

def procrustes_diffdim(Y_mat, X_mat, n_iter_max=10000, epsilon_min=1e-6, return_transformed=False):
    X = np.array(X_mat, dtype=np.double, copy=True)
    Y = np.array(Y_mat, dtype=np.double, copy=True)
    n_X, p_X = X.shape
    n_Y, p_Y = Y.shape
    assert n_X == n_Y
    assert p_X >= p_Y
    if p_X == p_Y:
        return procrustes(Y, X, return_transformed)
    else:
        Z = np.zeros((n_X, p_X - p_Y))
        for i in range(n_iter_max):
            W = np.hstack((Y, Z))
            R, rho, c = procrustes(W, X)
            X_new = X @ R * rho + c
            Z_new = X_new[:, p_Y:]
            Z_new_mean = np.mean(Z_new, 0)
            Z_new_centered = Z_new - Z_new_mean
            Z_diff = Z_new - Z
            epsilon = np.sum(Z_diff**2) / np.sum(Z_new_centered**2)
            if(epsilon < epsilon_min):
                break
            else:
                Z = Z_new
        if return_transformed:
            return R, rho, c, X_new
        else:
            return R, rho, c

def oadp(U, s, V, b, dim_ref=4, dim_stu=None, dim_online=None):
    if dim_stu is None:
        dim_stu = dim_ref * 2
    if dim_online is None:
        dim_online = dim_stu * 2
    pcs_ref = V[:, :dim_ref] * s[:dim_ref]
    s_aug, V_aug = svd_online(U[:,:dim_online], s[:dim_online], V[:,:dim_online], b)
    s_aug, V_aug = s_aug[:dim_stu], V_aug[:, :dim_stu]
    pcs_aug = V_aug * s_aug
    pcs_stu = ref_aug_procrustes(pcs_ref, pcs_aug)
    return pcs_stu[:dim_ref]






def source_fraposa_pca(genotypes_ref: np.array, genotypes_target= None, dim_ref=10, method='oadp', dim_stu=None, dim_online=None):
    # Setup the number of PCs we want to compute for difference methods
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

    print('PCA of REFERENCE samples...')
    pcs_ref = V[:, :dim_ref] * s[:dim_ref]

    # Process target samples
    if (genotypes_target is not None) and (type(genotypes_target) in [np.ndarray, da.Array]):
        print('Processing target samples...')
        # ToDo ensure alignment of variants between reference and target
        W = genotypes_target.astype(np.float64)
        p_ref = len(X_mean)
        p_stu, n_stu = W.shape
        pcs_stu = np.zeros((n_stu, dim_ref))

        for i in range(n_stu):
            w = W[:, i].astype(np.float64).reshape((-1, 1))
            w, _, _ = standardize(w, X_mean, X_std)
            if method == 'oadp':
                pcs_stu[i, :] = oadp(U, s, V, w, dim_ref, dim_stu, dim_online)
            if method == 'sp':
                pcs_stu[i, :] = w.T @ U[:, :dim_ref]
        del W
        return pcs_ref, pcs_stu
    else:
        return pcs_ref

#pca_fraposa = source_fraposa_pca(genotypes_ref=dd_ref[1])
pca_fraposa_ref, pca_fraposa_target = source_fraposa_pca(genotypes_ref=dd_ref[1], genotypes_target=dd_target[1], method='sp')
_, oadp_fraposa_target = source_fraposa_pca(genotypes_ref=dd_ref[1], genotypes_target=dd_target[1], method='oadp')
#np.savetxt('test_pca.txt', r_pca)


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
    XTX = X.T @ X
    s, V = svd_eigcov_dask(XTX)
    return s, V, XTX



def fraposa_dask(
    genotypes_ref,
    genotypes_target=None,
    dim_ref=10,
    method="oadp",
    dim_stu=None,
    dim_online=None
):
    """Implementation of fraposa_pgsc that supports Dask arrays for large-scale genotype data."""
    # FRAPOSA: Fast and Robust Ancestry Prediction by Online singular value decomposition and Shrinkage Adjustment
    # Original Author: David (Daiwei) Zhang
    # fraposa_pgsc modifications (QC, translation to dask): Samuel Lambert

    # Setup dimension requirements
    assert method in ["oadp", "sp"]

    if method == "oadp":
        if dim_stu is None:
            dim_stu = dim_ref * 2
        assert dim_ref <= dim_stu

        if dim_online is None:
            dim_online = dim_stu * 2
        assert dim_stu <= dim_online

    # Convert genotype matrix to Dask array if it isn't one already
    genotypes_ref = da.asarray(genotypes_ref, dtype=float)

    # 1. Standardize Reference
    X, X_mean, X_std = standardize_dask(genotypes_ref)

    # 2. Run methods
    if method == "oadp":
        X = X.compute()  # materialise the genotypes
        s, V = eig_ref(X)[:2]
        V = V[:, :dim_online]
        U = X @ (V / s[:dim_online])

        print('PCA of REFERENCE samples...')
        pcs_ref = V[:, :dim_ref] * s[:dim_ref]
        if (genotypes_target is not None) and (type(genotypes_target) in [np.ndarray, da.Array]):
            W = da.asarray(genotypes_target, dtype=float)
            W_std, _, _ = standardize_dask(W, mean=X_mean, std=X_std)
            W_std = W_std.compute()
            p_stu, n_stu = W_std.shape
            pcs_stu = np.zeros((n_stu, dim_ref))

            for i in range(n_stu):
                w = W_std[:, i].astype(np.float64).reshape((-1, 1))
                pcs_stu[i, :] = oadp(U, s, V, w, dim_ref, dim_stu, dim_online)
            return pcs_ref, pcs_stu
        else:
            return pcs_ref
    elif method == "sp":
        s, V, _ = eig_ref_dask(X) # This implementation uses a different solver, not compatiple with OADP!
        V = V[:, :dim_ref]
        U = X @ (V / s[:dim_ref])

        pcs_ref = V[:, :dim_ref] * s[:dim_ref]
        print('PCA of REFERENCE samples...')
        pcs_ref = pcs_ref.compute()

        if (genotypes_target is not None) and (type(genotypes_target) in [np.ndarray, da.Array]):
            print('Processing TARGET samples...')
            # ToDo ensure alignment of variants between reference and target
            W = da.asarray(genotypes_target, dtype=float)
            W_std, _, _ = standardize_dask(W, mean=X_mean, std=X_std)

            if method == 'sp':
                pcs_stu = W_std.T @ U[:, :dim_ref]
                pcs_stu = pcs_stu.compute()
            return pcs_ref, pcs_stu
        else:
            return pcs_ref

pca_dask_ref, pca_dask_target = fraposa_dask(genotypes_ref=dd_ref[1], genotypes_target=dd_target[1], method='sp')
oadp_dask_ref, oadp_dask_target = fraposa_dask(genotypes_ref=dd_ref[1], genotypes_target=dd_target[1], method='oadp')
