import argparse
import logging
import sys

from rich_argparse import RichHelpFormatter

from cat import run
from cat.constants import N_ITERATIONS, SIGMA, __version__


def init_logger(verbose: bool):
    """Initialize logger

    Parameters
    ----------
    verbose : bool
        Verbose mode
    """
    from rich.logging import RichHandler

    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(message)s",
        datefmt="%d-%m-%Y %H:%M:%S",
        handlers=[RichHandler()],
    )


def parse_args() -> argparse.Namespace:
    """Parse cli arguments

    Returns
    -------
    argparse.Namespace
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Cluster Alignment Tool (CAT)",
        formatter_class=RichHelpFormatter,
    )
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
        "--sigma",
        type=float,
        default=SIGMA,
        help=f"Sigma cutoff ({SIGMA} => p-value: 0.05)",
    )
    parser.add_argument(
        "--n_iter",
        type=int,
        default=N_ITERATIONS,
        help="Number of bootstraps, default 1,000",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Verbose mode",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"cat-python v{__version__}",
    )

    return parser.parse_args(args=None if sys.argv[1:] else ["--help"])


def main():
    """Main function"""
    args = parse_args()
    init_logger(args.verbose)
    run(args)


if __name__ == "__main__":
    main()
