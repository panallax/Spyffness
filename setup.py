from setuptools import setup, find_packages

setup(
    name="Spyffness",
    version="1.0.0",
    description="Direct Stiffness Matrix Solver and Post-Optimization for 3D truss structures",
    packages=find_packages(),
    install_requires=[
        "networkx",
        "scipy >= 1.11.0",
        "numpy",
        "matplotlib"
    ]
)