#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='ProtFlex',
    version='0.1.0',
    description='ProtFlex: a python standalone program to predict and assess the flexibility of proteins',
    author ='Maria Artigues-Lleixà and Pau Torrén-Duran',
    author_email='maria.artigues01@estudiant.upf.edu, pau.torren01@estudiant.upf.edu',
    url='https://github.com/mariaartlle/PYT-SBI-project-ProtFlex',
    packages=find_packages(),
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3"
    ],
)
