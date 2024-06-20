from setuptools import find_packages, setup

setup(
    name="Alleleome",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas == 2.2.2",
        "numpy==1.26.4",
        "biopython == 1.83",
    ],
    package_data={"Alleleome": ["sample_data/Oenococcus_oeni/*"]},
    entry_points={"console_scripts": ["alleleome=Alleleome.cli:main"]},
)
