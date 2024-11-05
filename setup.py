from setuptools import setup, find_packages
import os

# Fonction pour lire le fichier requirements.txt
def load_requirements(filename="requirements.txt"):
    with open(filename, "r") as f:
        return f.read().splitlines()

setup(
    name="vIPer",
    version="0.1.0",
    description="python code relative to the paper 'Estimating the vertical ionization potential of single-stranded DNA molecules' (M.Rooman and F. Pucci, submitted).",
    author="Benoitdw",
    url="https://github.com/Benoitdw/vIPer",
    packages=find_packages(),
    install_requires=load_requirements(),
    python_requires=">=3.6",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
