import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="osyris",
    version="1.0.5",
    author="Neil Vaytet",
    author_email="neil.vaytet@esss.se",
    description="A package to visualize AMR data from the RAMSES code",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nvaytet/osyris",
    packages=setuptools.find_packages("src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "matplotlib>=2.0.0",
        "numpy>=1.0.0",
        "scipy>=0.1.0",
        "ipyvolume>=0.5",
    ],
)