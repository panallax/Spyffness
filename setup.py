from setuptools import setup, find_packages

setup(
    name="Spyffness",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "networkx",
        "scipy",
    ],
    py_modules=["config"],
    author="Alex Bernadí",
    description="A Python library for structural analysis",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/panallax/Spyffness",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)