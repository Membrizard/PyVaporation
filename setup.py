from distutils.core import setup

setup(
    name="pyvaporation",
    packages=["pyvaporation"],
    version="0.2.0",
    license="Apache license 2.0",
    description="Set of tools for modelling pervaporation processes",
    author="Denis Sapegin, Aleksei Chekmachev",
    author_email="a.checkmachev@gmail.com",
    url="https://github.com/Membrizard/PyVaporation",
    download_url="https://github.com/Membrizard/PyVaporation/archive/refs/tags/v0.2.0.tar.gz",
    long_description_content_type="text/markdown",
    keywords=[
        "pervaporation",
        "membrane",
        "chemistry",
        "modelling",
    ],
    install_requires=[
        "joblib>=1.1.0",
        "matplotlib>=3.3.4",
        "pandas>=1.1.5",
        "scipy>=1.5.4",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Build Tools",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
)
