import logging
import pandas as pd
import numpy as np
from sklearn.covariance import MinCovDet, EmpiricalCovariance
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LinearRegression, GammaRegressor
from scipy.stats import chi2, percentileofscore, mannwhitneyu
from scipy import optimize as opt
import json
import gzip


logger = logging.getLogger(__name__)

####
# Methods for comparing sample ancestry to a reference dataset using PCA data & population labels
###

comparison_method_threshold = {
    "Mahalanobis": 1e-10,
    "RandomForest": 0.5,
}  # default p-value thresholds to define low-confidence population matching
_mahalanobis_methods = ["MinCovDet", "EmpiricalCovariance"]


def choose_pval_threshold(args):
    set_threshold = comparison_method_threshold[args.method_compare]  # method default
    if args.pThreshold is not None:
        if (args.pThreshold > 0) and (args.pThreshold < 1):
            set_threshold = args.pThreshold
        else:
            logging.warning(
                "p-value threshold out of range, assigning as method default: {}"
            )
    return set_threshold


def get_covariance_method(method_name):
    match method_name:
        case "MinCovDet":
            covariance_model = MinCovDet()
        case "EmpiricalCovariance":
            covariance_model = EmpiricalCovariance()
        case _:
            assert False, "Invalid covariance method"

    return covariance_model


def compare_ancestry(
    ref_df: pd.DataFrame,
    ref_pop_col: str,
    target_df: pd.DataFrame,
    ref_train_col=None,
    n_pcs=4,
    method="RandomForest",
    covariance_method="EmpiricalCovariance",
    p_threshold=None,
):
    """
    Function to compare target sample ancestry to a reference panel with PCA data
    :param ref_df: reference dataset
    :param ref_pop_col: training labels for population assignment in reference dataset
    :param target_df: target dataset
    :param ref_train_col: column name of TRUE/FALSE labels for inclusion in training ancestry assignments (e.g. unrelated)
    :param n_pcs: number of PCs to use in ancestry assignment
    :param method: One of Mahalanobis or RandomForest
    :param covariance_method: Used to calculate Mahalanobis distances One of EmpiricalCovariance or MinCovDet
    :param p_threshold: used to define LowConfidence population assignments
    :return: dataframes for reference (predictions on training set) and target (predicted labels) datasets
    """
    # Check that datasets have the correct columns
    assert (
        method in comparison_method_threshold.keys()
    ), "comparison method parameter must be Mahalanobis or RF"
    if method == "Mahalanobis":
        assert (
            covariance_method in _mahalanobis_methods
        ), "covariance estimation method must be MinCovDet or EmpiricalCovariance"

    cols_pcs = ["PC{}".format(x + 1) for x in range(0, n_pcs)]
    assert all(
        [col in ref_df.columns for col in cols_pcs]
    ), "Reference Dataset (ref_df) is missing some PC columns for ancestry comparison (max:{})".format(
        n_pcs
    )
    assert all(
        [col in target_df.columns for col in cols_pcs]
    ), "Target Dataset (target_df) is missing some PC columns for ancestry comparison (max:{})".format(
        n_pcs
    )
    assert (
        ref_pop_col in ref_df.columns
    ), "Population label column ({}) is missing from reference dataframe".format(
        ref_pop_col
    )
    ref_populations = ref_df[ref_pop_col].unique()

    # Extract columns for analysis
    ref_df = ref_df[cols_pcs + [ref_pop_col, ref_train_col]].copy()
    target_df = target_df[cols_pcs].copy()

    # Create Training dfs
    if ref_train_col:
        assert (
            ref_train_col in ref_df.columns
        ), "Training index column({}) is missing from reference dataframe".format(
            ref_train_col
        )
        ref_train_df = ref_df.loc[ref_df[ref_train_col],]
    else:
        ref_train_df = ref_df

    # Check outlier-ness of target with regard to the reference PCA space
    compare_info = {}
    pop = "ALL"
    ref_covariance_model = get_covariance_method(covariance_method)
    ref_covariance_fit = ref_covariance_model.fit(ref_train_df[cols_pcs])
    colname_dist = "Mahalanobis_dist_{}".format(pop)
    colname_pval = "Mahalanobis_P_{}".format(pop)
    target_df[colname_dist] = ref_covariance_fit.mahalanobis(target_df[cols_pcs])
    target_df[colname_pval] = chi2.sf(target_df[colname_dist], n_pcs - 1)
    compare_info["Mahalanobis_P_ALL"] = dict(target_df[colname_pval].describe())
    logger.info(
        "Mahalanobis Probability Distribution (train: all reference samples): {}".format(
            compare_info["Mahalanobis_P_ALL"]
        )
    )

    ## Check if PCs only capture target/reference stratification
    if target_df.shape[0] >= 20:
        for col_pc in cols_pcs:
            mwu_pc = mannwhitneyu(ref_train_df[col_pc], target_df[col_pc])
            compare_info[col_pc] = {"U": mwu_pc.statistic, "pvalue": mwu_pc.pvalue}
            if mwu_pc.pvalue < 1e-4:
                logger.warning(
                    "{} *may* be capturing target/reference stratification (Mann-Whitney p-value={}), "
                    "use visual inspection of PC plot to confirm".format(
                        col_pc, mwu_pc.pvalue
                    )
                )

    # Run Ancestry Assignment methods
    if method == "Mahalanobis":
        logger.debug("Calculating Mahalanobis distances")
        # Calculate population distances
        pval_cols = []
        for pop in ref_populations:
            logger.debug("Fitting Mahalanobis distances: {}".format(pop))
            # Fit the covariance matrix for the current population
            colname_dist = "Mahalanobis_dist_{}".format(pop)
            colname_pval = "Mahalanobis_P_{}".format(pop)

            covariance_model = get_covariance_method(covariance_method)

            covariance_fit = covariance_model.fit(
                ref_train_df.loc[ref_train_df[ref_pop_col] == pop, cols_pcs]
            )

            # Caclulate Mahalanobis distance of each sample to that population
            # Reference Samples
            ref_df[colname_dist] = covariance_fit.mahalanobis(ref_df[cols_pcs])
            ref_df[colname_pval] = chi2.sf(ref_df[colname_dist], n_pcs - 1)
            # Target Samples
            target_df[colname_dist] = covariance_fit.mahalanobis(target_df[cols_pcs])
            target_df[colname_pval] = chi2.sf(target_df[colname_dist], n_pcs - 1)

            pval_cols.append(colname_pval)

        # Assign population (maximum probability)
        logger.debug("Assigning Populations (max Mahalanobis probability)")
        ref_assign = ref_df[pval_cols].copy()
        ref_assign = ref_assign.assign(MostSimilarPop=ref_assign.idxmax(axis=1))

        target_assign = target_df[pval_cols].copy()
        target_assign = target_assign.assign(
            MostSimilarPop=target_assign.idxmax(axis=1)
        )

        ref_assign["MostSimilarPop_LowConfidence"] = np.nan
        target_assign["MostSimilarPop_LowConfidence"] = np.nan

        if p_threshold:
            logger.debug(
                "Comparing Population Similarity to p-value threshold (p < {})".format(
                    p_threshold
                )
            )
            ref_assign["MostSimilarPop_LowConfidence"] = [
                (ref_assign[x][i] < p_threshold)
                for i, x in enumerate(ref_assign["MostSimilarPop"])
            ]
            target_assign["MostSimilarPop_LowConfidence"] = [
                (target_assign[x][i] < p_threshold)
                for i, x in enumerate(target_assign["MostSimilarPop"])
            ]

        # Cleanup variable names
        ref_assign["MostSimilarPop"] = [
            x.split("_")[-1] for x in ref_assign["MostSimilarPop"]
        ]
        target_assign["MostSimilarPop"] = [
            x.split("_")[-1] for x in target_assign["MostSimilarPop"]
        ]

    elif method == "RandomForest":
        # Assign SuperPop Using Random Forest (PCA loadings)
        logger.debug("Training RandomForest classifier")
        clf_rf = RandomForestClassifier(random_state=32)
        clf_rf.fit(
            ref_train_df[cols_pcs],  # X (training PCs)
            ref_train_df[ref_pop_col].astype(str),
        )  # Y (pop label)

        # Predict most similar population using RF classifier
        logger.debug("Find most similar Populations (max RF probability)")
        ref_assign = pd.DataFrame(
            clf_rf.predict_proba(ref_df[cols_pcs]),
            index=ref_df.index,
            columns=["RF_P_{}".format(x) for x in clf_rf.classes_],
        )
        ref_assign["MostSimilarPop"] = clf_rf.predict(ref_df[cols_pcs])

        target_assign = pd.DataFrame(
            clf_rf.predict_proba(target_df[cols_pcs]),
            index=target_df.index,
            columns=["RF_P_{}".format(x) for x in clf_rf.classes_],
        )
        target_assign["MostSimilarPop"] = clf_rf.predict(target_df[cols_pcs])

        # Define confidence using p-value thresholds
        ref_assign["MostSimilarPop_LowConfidence"] = np.nan
        target_assign["MostSimilarPop_LowConfidence"] = np.nan

        if p_threshold:
            logger.debug(
                "Comparing Population Similarity to p-value threshold (p < {})".format(
                    p_threshold
                )
            )
            ref_assign["MostSimilarPop_LowConfidence"] = [
                (ref_assign["RF_P_{}".format(x)].iloc[i] < p_threshold)
                for i, x in enumerate(clf_rf.predict(ref_df[cols_pcs]))
            ]
            target_assign["MostSimilarPop_LowConfidence"] = [
                (target_assign["RF_P_{}".format(x)].iloc[i] < p_threshold)
                for i, x in enumerate(clf_rf.predict(target_df[cols_pcs]))
            ]

    return ref_assign, target_assign, compare_info


####
# Methods for adjusting/reporting polygenic score results that account for genetic ancestry
####
normalization_methods = ["empirical", "mean", "mean+var"]


def pgs_adjust(
    ref_df,
    target_df,
    scorecols: list,
    ref_pop_col,
    target_pop_col,
    use_method: list,
    norm2_2step=False,
    ref_train_col=None,
    n_pcs=4,
    norm_centerpgs=True,
    std_pcs=True,
):
    """
    Function to adjust PGS using population references and/or genetic ancestry (PCs)
    :param ref_df: reference dataset
    :param target_df: datagframe with target samples
    :param scorecols: [list of columns containing PGS that should be adjusted]
    :param ref_pop_col: ref_df column with population labels
    :param target_pop_col: target_df column with population labels that will be matched to ref_df population distributions
    :param use_method: list of ["empirical", "mean", "mean+var"]
    :param norm2_2step: boolean (default=False) whether to use the two-step model vs. the full-fit
    :param ref_train_col: column name with true/false labels of samples that should be included in training PGS methods
    :param n_pcs: number of genetic PCs that will be used for PGS-adjustment
    :return: [results_ref:df , results_target:df , results_models: dict] adjusted dfs for reference and target populations, and a dictionary with model fit/parameters.
    """
    # Check that datasets have the correct columns
    ## Check that score is in both dfs
    assert all(
        [x in ref_df.columns for x in scorecols]
    ), "Reference Dataset (ref_df) is missing some PGS column(s)".format()
    assert all(
        [x in target_df.columns for x in scorecols]
    ), "Target Dataset (target_df) is missing some PGS column(s)".format()

    ## Check that PCs is in both dfs
    cols_pcs = ["PC{}".format(x + 1) for x in range(0, n_pcs)]
    assert all(
        [col in ref_df.columns for col in cols_pcs]
    ), "Reference Dataset (ref_df) is missing some PC columns for PCA adjustment (max:{})".format(
        n_pcs
    )
    assert all(
        [col in target_df.columns for col in cols_pcs]
    ), "Target Dataset (target_df) is missing some PC columns for PCA adjustment (max:{})".format(
        n_pcs
    )
    assert (
        ref_pop_col in ref_df.columns
    ), "Population label column ({}) is missing from reference dataframe".format(
        ref_pop_col
    )
    ref_populations = ref_df[ref_pop_col].unique()
    assert (
        target_pop_col in target_df.columns
    ), "Population label column ({}) is missing from target dataframe".format(
        target_pop_col
    )

    ## Create Training dfs
    if ref_train_col:
        assert (
            ref_train_col in ref_df.columns
        ), "Training index column({}) is missing from reference dataframe".format(
            ref_train_col
        )
        ref_train_df = ref_df.loc[ref_df[ref_train_col],].copy()
    else:
        ref_train_df = ref_df.copy()

    ## Create results structures
    results_ref = {}
    results_target = {}
    results_models = {}  # used to store regression information
    scorecols_drop = set()
    for c_pgs in scorecols:
        # Makes melting easier later
        sum_col = "SUM|{}".format(c_pgs)
        results_ref[sum_col] = ref_df[c_pgs]
        results_target[sum_col] = target_df[c_pgs]
        results_models = {}

        # Check that PGS has variance (e.g. not all 0)
        if 0 in [np.var(results_ref[sum_col]), np.var(results_target[sum_col])]:
            scorecols_drop.add(c_pgs)
            logger.warning(
                "Skipping adjustment: {} has 0 variance in PGS SUM".format(c_pgs)
            )

    # Report PGS values with respect to distribution of PGS in the most similar reference population
    if "empirical" in use_method:
        logger.debug(
            "Adjusting PGS using most similar reference population distribution."
        )
        results_models["dist_empirical"] = {}

        for c_pgs in scorecols:
            # Initialize Output
            percentile_col = "percentile_MostSimilarPop|{}".format(c_pgs)
            results_ref[percentile_col] = pd.Series(index=ref_df.index, dtype="float64")
            results_target[percentile_col] = pd.Series(
                index=target_df.index, dtype="float64"
            )
            z_col = "Z_MostSimilarPop|{}".format(c_pgs)
            results_ref[z_col] = pd.Series(index=ref_df.index, dtype="float64")
            results_target[z_col] = pd.Series(index=target_df.index, dtype="float64")

            if c_pgs not in scorecols_drop:
                r_model = {}

                # Adjust for each population
                for pop in ref_populations:
                    r_pop = {}
                    i_ref_pop = ref_df[ref_pop_col] == pop
                    i_target_pop = target_df[target_pop_col] == pop

                    # Reference Score Distribution
                    c_pgs_pop_dist = ref_train_df.loc[
                        ref_train_df[ref_pop_col] == pop, c_pgs
                    ]

                    # Calculate Percentile
                    results_ref[percentile_col].loc[i_ref_pop] = percentileofscore(
                        c_pgs_pop_dist, ref_df.loc[i_ref_pop, c_pgs]
                    )
                    results_target[percentile_col].loc[
                        i_target_pop
                    ] = percentileofscore(
                        c_pgs_pop_dist, target_df.loc[i_target_pop, c_pgs]
                    )
                    r_pop["percentiles"] = np.percentile(
                        c_pgs_pop_dist, range(0, 101, 1)
                    )

                    # Calculate Z
                    r_pop["mean"] = c_pgs_pop_dist.mean()
                    r_pop["std"] = c_pgs_pop_dist.std(ddof=0)

                    results_ref[z_col].loc[i_ref_pop] = (
                        ref_df.loc[i_ref_pop, c_pgs] - r_pop["mean"]
                    ) / r_pop["std"]
                    results_target[z_col].loc[i_target_pop] = (
                        target_df.loc[i_target_pop, c_pgs] - r_pop["mean"]
                    ) / r_pop["std"]

                    r_model[pop] = r_pop

                results_models["dist_empirical"][c_pgs] = r_model
                # ToDo: explore handling of individuals who have low-confidence population labels
                #  -> Possible Soln: weighted average based on probabilities? Small Mahalanobis P-values will complicate this
    # PCA-based adjustment
    if any([x in use_method for x in ["mean", "mean+var"]]):
        logger.debug("Adjusting PGS using PCA projections")
        results_models["adjust_pcs"] = {"PGS": {}}

        # Make copies of ref/target dfs for normalizing
        normcols = scorecols + cols_pcs
        ref_norm = ref_df[normcols].copy()
        target_norm = target_df[normcols].copy()

        if std_pcs:
            pcs_norm = {}
            for pc_col in cols_pcs:
                # Calculate norm factors
                pc_mean = ref_train_df[pc_col].mean()
                pc_std = ref_train_df[pc_col].std(ddof=0)
                pcs_norm[pc_col] = {"mean": pc_mean, "pc_std": pc_std}
                results_models["adjust_pcs"]["norm_pcs"] = pcs_norm

                # Normalize data
                ref_train_df[pc_col] = (ref_train_df[pc_col] - pc_mean) / pc_std
                ref_norm[pc_col] = (ref_norm[pc_col] - pc_mean) / pc_std
                target_norm[pc_col] = (target_norm[pc_col] - pc_mean) / pc_std

        for c_pgs in scorecols:
            if c_pgs in scorecols_drop:
                # fill the output with NAs
                adj_cols = ["Z_norm1|{}".format(c_pgs)]
                if "mean+var" in use_method:
                    adj_cols.append("Z_norm2|{}".format(c_pgs))
                for adj_col in adj_cols:
                    results_ref[adj_col] = pd.Series(
                        index=ref_df.index, dtype="float64"
                    )  # fill na
                    results_target[adj_col] = pd.Series(
                        index=target_df.index, dtype="float64"
                    )  # fill na
            else:
                results_models["adjust_pcs"]["PGS"][c_pgs] = {}
                if norm_centerpgs:
                    pgs_mean = ref_train_df[c_pgs].mean()
                    ref_train_df[c_pgs] = ref_train_df[c_pgs] - pgs_mean
                    ref_norm[c_pgs] = ref_norm[c_pgs] - pgs_mean
                    target_norm[c_pgs] = target_norm[c_pgs] - pgs_mean
                    results_models["adjust_pcs"]["PGS"][c_pgs]["pgs_offset"] = pgs_mean

                # Method 1 (Khera et al. Circulation (2019): normalize mean (doi:10.1161/CIRCULATIONAHA.118.035658)
                adj_col = "Z_norm1|{}".format(c_pgs)
                # Fit to Reference Data
                pcs2pgs_fit = LinearRegression().fit(
                    ref_train_df[cols_pcs], ref_train_df[c_pgs]
                )
                ref_train_pgs_pred = pcs2pgs_fit.predict(ref_train_df[cols_pcs])
                ref_train_pgs_resid = ref_train_df[c_pgs] - ref_train_pgs_pred
                ref_train_pgs_resid_mean = ref_train_pgs_resid.mean()
                ref_train_pgs_resid_std = ref_train_pgs_resid.std(ddof=0)

                ref_pgs_resid = ref_norm[c_pgs] - pcs2pgs_fit.predict(
                    ref_norm[cols_pcs]
                )
                results_ref[adj_col] = ref_pgs_resid / ref_train_pgs_resid_std
                # Apply to Target Data
                target_pgs_pred = pcs2pgs_fit.predict(target_norm[cols_pcs])
                target_pgs_resid = target_norm[c_pgs] - target_pgs_pred
                results_target[adj_col] = target_pgs_resid / ref_train_pgs_resid_std
                results_models["adjust_pcs"]["PGS"][c_pgs][
                    "Z_norm1"
                ] = package_skl_regression(pcs2pgs_fit)

                if "mean+var" in use_method:
                    # Method 2 (Khan et al. Nature Medicine (2022)): normalize variance (doi:10.1038/s41591-022-01869-1)
                    # Normalize based on residual deviation from mean of the distribution [equalize population sds]
                    # (e.g. reduce the correlation between genetic ancestry and how far away you are from the mean)
                    # USE gamma distribution for predicted variance to constrain it to be positive (b/c using linear
                    # regression we can get negative predictions for the sd)
                    adj_col = "Z_norm2|{}".format(c_pgs)
                    pcs2var_fit_gamma = GammaRegressor(max_iter=1000).fit(
                        ref_train_df[cols_pcs],
                        (ref_train_pgs_resid - ref_train_pgs_resid_mean) ** 2,
                    )
                    if norm2_2step:
                        # Return 2-step adjustment
                        results_ref[adj_col] = ref_pgs_resid / np.sqrt(
                            pcs2var_fit_gamma.predict(ref_norm[cols_pcs])
                        )
                        results_target[adj_col] = target_pgs_resid / np.sqrt(
                            pcs2var_fit_gamma.predict(target_norm[cols_pcs])
                        )
                        results_models["adjust_pcs"]["PGS"][c_pgs][
                            "Z_norm2"
                        ] = package_skl_regression(pcs2var_fit_gamma)
                    else:
                        # Return full-likelihood adjustment model
                        # This jointly re-fits the regression parameters from the mean and variance prediction to better
                        # fit the observed PGS distribution. It seems to mostly change the intercepts. This implementation is
                        # adapted from https://github.com/broadinstitute/palantir-workflows/blob/v0.14/ImputationPipeline/ScoringTasks.wdl,
                        # which is distributed under a BDS-3 license.
                        params_initial = np.concatenate(
                            [
                                [pcs2pgs_fit.intercept_],
                                pcs2pgs_fit.coef_,
                                [pcs2var_fit_gamma.intercept_],
                                pcs2var_fit_gamma.coef_,
                            ]
                        )
                        pcs2full_fit = fullLL_fit(
                            df_score=ref_train_df,
                            scorecol=c_pgs,
                            predictors=cols_pcs,
                            initial_params=params_initial,
                        )

                        results_ref[adj_col] = fullLL_adjust(
                            pcs2full_fit, ref_norm, c_pgs
                        )
                        results_target[adj_col] = fullLL_adjust(
                            pcs2full_fit, target_norm, c_pgs
                        )

                        if pcs2full_fit["params"]["success"] is False:
                            logger.warning(
                                "{} full-likelihood: {} {}".format(
                                    c_pgs,
                                    pcs2full_fit["params"]["status"],
                                    pcs2full_fit["params"]["message"],
                                )
                            )
                        results_models["adjust_pcs"]["PGS"][c_pgs][
                            "Z_norm2"
                        ] = pcs2full_fit
    # Only return results
    logger.debug("Outputting adjusted PGS & models")
    results_ref = pd.DataFrame(results_ref)
    results_target = pd.DataFrame(results_target)
    return results_ref, results_target, results_models


def f_var(df, beta):
    """Predict the result of a Gamma regression (log link) w/ intercept"""
    return np.exp(beta[0] + np.inner(beta[1:], df))


def f_mu(df, beta):
    """Predict the result of a Gaussian regression w/ intercept"""
    return beta[0] + np.inner(beta[1:], df)


def nLL_mu_and_var(theta, df, c_score, l_predictors):
    """Negative log-likelihood for regression that fits the expected mean and variance based on a
    set of predictors. In this case fitting the PGS as a function of geneitc ancestry (PCA loadings).
    Adapted from https://github.com/broadinstitute/palantir-workflows/blob/v0.14/ImputationPipeline/ScoringTasks.wdl,
    which is distributed under a BDS-3 license."""
    i_split = int(1 + (len(theta) - 2) / 2)
    theta_mu = theta[0:i_split]
    theta_var = theta[i_split:]
    x = df[c_score]
    return sum(
        np.log(np.sqrt(f_var(df[l_predictors], theta_var)))
        + (1 / 2)
        * (x - f_mu(df[l_predictors], theta_mu)) ** 2
        / f_var(df[l_predictors], theta_var)
    )


def grdnt_mu_and_var(theta, df, c_score, l_predictors):
    """Gradient used to optimize the nLL_mu_and_var fit function.
    Adapted from https://github.com/broadinstitute/palantir-workflows/blob/v0.14/ImputationPipeline/ScoringTasks.wdl,
    which is distributed under a BDS-3 license."""
    i_split = int(1 + (len(theta) - 2) / 2)
    beta_mu = theta[0:i_split]
    beta_var = theta[i_split:]
    x = df[c_score]

    pred_var = f_var(df[l_predictors], beta_var)  # current prediction of variance

    mu_coeff = -(x - f_mu(df[l_predictors], beta_mu)) / f_var(
        df[l_predictors], beta_var
    )
    sig_coeff = 1 / (2 * f_var(df[l_predictors], beta_var)) - (1 / 2) * (
        x - f_mu(df[l_predictors], beta_mu)
    ) ** 2 / (f_var(df[l_predictors], beta_var) ** 2)

    grad = np.concatenate(
        [
            [sum(mu_coeff * 1)],
            [sum(mu_coeff * df[x]) for x in l_predictors],
            [sum(sig_coeff * (1 * pred_var))],
            [sum(sig_coeff * (df[x] * pred_var)) for x in l_predictors],
        ]
    )

    return grad


def fullLL_fit(df_score, scorecol, predictors, initial_params):
    """Fit the regression of values based on predictors that effect the mean and variance of the distribution"""
    fit_result = opt.minimize(
        fun=nLL_mu_and_var,
        x0=initial_params,
        method="BFGS",
        jac=grdnt_mu_and_var,
        args=(df_score, scorecol, predictors),
        options={"maxiter": 1000},
    )

    # package result for output and use in prediction
    rp = dict(fit_result)
    rp["initial_params"] = initial_params
    x = rp.pop("x")  # fitted coefficients

    i_split = int(1 + (len(x) - 2) / 2)
    x_mu = x[0:i_split]
    x_var = x[i_split:]
    return {
        "params": rp,
        "mu_intercept": x_mu[0],
        "mu_coef": dict(zip(predictors, x_mu[1:])),
        "var_intercept": x_var[0],
        "var_coef": dict(zip(predictors, x_var[1:])),
    }


def fullLL_adjust(fullLL_model, df_score, scorecol):
    """Function to adjust PGS based on the full likelihood fit"""
    predictors = fullLL_model["mu_coef"].keys()
    mu_coef = np.concatenate(
        [[fullLL_model["mu_intercept"]], list(fullLL_model["mu_coef"].values())]
    )
    var_coef = np.concatenate(
        [[fullLL_model["var_intercept"]], list(fullLL_model["var_coef"].values())]
    )
    return (df_score[scorecol] - f_mu(df_score[predictors], mu_coef)) / np.sqrt(
        f_var(df_score[predictors], var_coef)
    )


def package_skl_regression(model):
    """Extract relevant details from sklearn regression model"""
    return {
        "model": str(type(model)),
        "params": model.get_params(),
        "_intercept": model.intercept_,
        "_coef": dict(zip(model.feature_names_in_, model.coef_)),
    }


class NumpyEncoder(json.JSONEncoder):
    """Special json encoder for numpy types (taken from:
    https://stackoverflow.com/questions/26646362/numpy-array-is-not-json-serializable)"""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def write_model(d: dict, outname):
    """Use numpy encoder to write json file for models"""
    logger.debug("Writing PGS adjustment models to: {}".format(outname))
    if outname.endswith(".gz"):
        outfile = gzip.open(outname, "wt")
    else:
        outfile = open(outname, "w")

    outfile.write(json.dumps(d, indent=2, cls=NumpyEncoder))
    outfile.close()
