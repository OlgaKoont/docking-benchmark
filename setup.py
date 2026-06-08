"""Install ToxAffinity docking pipeline package."""

from pathlib import Path

from setuptools import find_packages, setup

readme = Path(__file__).parent / "README.md"
long_description = readme.read_text(encoding="utf-8") if readme.exists() else ""

setup(
    name="toxaffinity",
    version="1.0.0",
    description="ToxDock-Bench: toxicity-oriented docking benchmark pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.9",
    install_requires=[
        "pandas>=1.3.0",
        "numpy<2.0",
        "matplotlib>=3.3.0",
        "seaborn>=0.11.0",
        "pyyaml>=5.4.0",
        "scipy>=1.7.0",
        "meeko>=0.5.0",
        "biopython>=1.79",
    ],
    entry_points={
        "console_scripts": [
            "toxdock-pipeline=docking_benchmark2.cli.run_benchmark:main",
        ],
    },
)
