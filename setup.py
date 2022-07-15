from setuptools import setup, find_packages

setup(
    name="pyvaporation",
    packages=find_packages(),
    version="1.1.2",
    license="Apache license 2.0",
    description="Set of tools for modelling pervaporation processes",
    author="Denis Sapegin, Aleksei Chekmachev",
    author_email="a.checkmachev@gmail.com",
    url="https://github.com/Membrizard/PyVaporation",
    download_url="https://github.com/Membrizard/PyVaporation/archive/refs/tags/v1.1.2.tar.gz",
    long_description=open("README.md", "r").read(),
    long_description_content_type="text/markdown",
    keywords=[
        "pervaporation",
        "membrane",
        "chemistry",
        "modelling",
    ],
    install_requires=[
        "joblib>=1.1.0",
        "matplotlib>=3.5.2",
        "pandas>=1.3.5",
        "scipy>=1.7.3",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Build Tools",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
)
