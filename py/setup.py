# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

VERSION = '0.2.2'

setup(
    name="easySCFpy",  # package name
    version=VERSION,  # package version
    author="Lei Cui",
    author_email="cuilei798@qq.com",
    maintainer="Lei Cui",
    maintainer_email="cuilei798@qq.com",
    license="MIT License",
    platforms=["linux"],
    url="https://github.com/xleizi/easySCFpy",
    description="The easySCFpy is a Python package for transformating single-cell data.",
    long_description=open("README.md").read(),
    packages=find_packages(),
    zip_safe=False,
    # entry_points={
    #     "console_scripts": [
    #         "easyBio=easyBio.easyBio:main",
    #     ]
    # },
    install_requires=["anndata", "h5py", "scipy", "numpy"],
    package_data={"Utils": ["Utils/*"]},
    classifiers=[
        "Operating System :: OS Independent",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    python_requires=">=3",
)
