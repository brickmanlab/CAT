import argparse
import warnings

from . import __version__
from .cat import run
from .utils import init_logger

warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=DeprecationWarning)
warnings.simplefilter(action="ignore", category=RuntimeWarning)


def main():
    parser = argparse.ArgumentParser(description=f"Cluster Alignment Tool (CAT)")
    parser.add_argument(
        "--ds1",
        action="store",
        type=str,
        help="Processed dataset (h5/h5ad)",
    )
    parser.add_argument(
        "--ds1_name",
        type=str,
        help="Dataset name",
    )
    parser.add_argument(
        "--ds1_cluster",
        type=str,
        help="Column name for comparison",
    )
    parser.add_argument(
        "--ds1_genes",
        type=str,
        default=None,
        help="Gene column, using `index` as default",
    )
    parser.add_argument(
        "--ds2",
        action="store",
        type=str,
        help="Processed dataset (h5/h5ad)",
    )
    parser.add_argument(
        "--ds2_name",
        type=str,
        help="Dataset name",
    )
    parser.add_argument(
        "--ds2_cluster",
        type=str,
        help="Column name for comparison",
    )
    parser.add_argument(
        "--ds2_genes",
        type=str,
        default=None,
        help="Gene column, using `index` as default",
    )
    parser.add_argument(
        "--features",
        type=str,
        default=None,
        help="File containing list of genes on new lines",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="./results",
        help="Output location",
    )
    parser.add_argument(
        "--distance", type=str, default="euclidean", help="Distance measurement"
    )
    parser.add_argument(
        "--sigma", type=float, default=1.6, help="Sigma cutoff (1.6 => p-value: 0.05)"
    )
    parser.add_argument(
        "--n_iter",
        type=int,
        default=1_000,
        help="Number of bootstraps (default: 1,000)",
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="Number of threads (default: 1)"
    )
    parser.add_argument(
        "--format",
        type=str,
        default="excel",
        choices=["excel", "html"],
        help="Report output format",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Verbose mode",
    )
    parser.add_argument("--version", action="version", version=f"CAT v{__version__}")

    args = parser.parse_args()
    init_logger(args.verbose)
    run(args)


if __name__ == "__main__":
    main()
