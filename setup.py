from setuptools import find_packages, setup

__version__ = "1.0"

with open("README.md") as f:
    README = f.read()

with open("requirements.txt") as f:
    INSTALL_REQUIREMENTS = f.read().splitlines()

setup(
    name="cat-python",
    version=__version__,
    python_requires=">=3.7.*",
    description="Cluster Alignment Tool",
    long_description=README,
    long_description_content_type="text/markdown",
    author="Alex Valentin Nielsen, Martin Proks, Ala Trusina",
    author_email="alexander.nielsen@nbi.ku.dk",
    url="https://github.com/brickmanlab/CAT",
    license="MIT",
    packages=find_packages(),
    install_requires=INSTALL_REQUIREMENTS,
    zip_safe=False,
    entry_points={"console_scripts": ["catcli=cat.__main__:main"]},
    classifiers=[
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
