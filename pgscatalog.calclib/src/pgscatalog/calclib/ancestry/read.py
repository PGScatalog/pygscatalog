import logging
import pandas as pd


logger = logging.getLogger(__name__)


def read_pcs(loc_pcs: list[str], dataset: str, loc_related_ids=None, nPCs=None):
    """
    Read the .pc file outputs of the fraposa_pgsc projection
    :param loc_pcs: list of locations for .pcs files
    :param dataset: name of the dataset being read (used for index)
    :param loc_related_ids: path to newline-delimited list of IDs for related samples that can be used to filter
    :return: pandas dataframe with PC information
    """
    proj = pd.DataFrame()

    for i, path in enumerate(loc_pcs):
        logger.debug("Reading PCA projection: {}".format(path))
        df = pd.read_csv(path, sep="\t", converters={"IID": str}, header=0)
        df["sampleset"] = dataset
        df.set_index(["sampleset", "IID"], inplace=True)

        if i == 0:
            logger.debug("Initialising combined DF")
            proj = df.copy()
        else:
            logger.debug("Appending to combined DF")
            proj = pd.concat([proj, df])

    # Drop PCs
    if nPCs:
        logger.debug("Filtering to relevant PCs")
        dropcols = []
        for x in proj.columns:
            if int(x[2:]) > nPCs:
                dropcols.append(x)
        proj = proj.drop(dropcols, axis=1)

    # Read/process IDs for unrelated samples (usually reference dataset)
    if loc_related_ids:
        logger.debug("Flagging related samples with: {}".format(loc_related_ids))
        proj["Unrelated"] = True
        with open(loc_related_ids, "r") as infile:
            IDs_related = [x.strip() for x in infile.readlines()]
        proj.loc[
            proj.index.get_level_values(level=1).isin(IDs_related), "Unrelated"
        ] = False
    else:
        # if unrelated is all nan -> dtype is float64
        # if unrelated is only true / false -> dtype is bool
        # if unrelated contains None, dtype stays bool, and pd.concat warning disappears
        proj["Unrelated"] = None

    return proj


def extract_ref_psam_cols(
    loc_psam, dataset: str, df_target, keepcols=("SuperPop", "Population")
):
    psam = pd.read_csv(loc_psam, sep="\t", header=0)

    match psam.columns[0]:
        # handle case of #IID -> IID (happens when #FID is present)
        case "#IID":
            psam.rename({"#IID": "IID"}, axis=1, inplace=True)
        case "#FID":
            psam.drop(["#FID"], axis=1, inplace=True)
        case _:
            assert False, "Invalid columns"
    psam["sampleset"] = dataset
    psam.set_index(["sampleset", "IID"], inplace=True)

    return pd.merge(df_target, psam[keepcols], left_index=True, right_index=True)


def read_pgs(loc_aggscore):
    """
    Function to read the output of aggreagte_scores
    :param loc_aggscore: path to aggregated scores output
    :return:
    """
    logger.debug("Reading aggregated score data: {}".format(loc_aggscore))
    df = pd.read_csv(
        loc_aggscore,
        sep="\t",
        index_col=["sampleset", "IID"],
        converters={"IID": str},
        header=0,
    ).pivot(columns=["PGS"], values=["SUM", "AVG"])
    # join column levels ({PGS}_{VALUE})
    df.columns = [f"{j}_{i}" for i, j in df.columns]

    return df
