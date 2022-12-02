from setuptools import find_packages, setup

setup(
    name="pyvaporation",
    packages=find_packages(),
    version="1.1.10",
    license="Apache license 2.0",
    description="Set of tools for modelling pervaporation processes",
    author="Denis Sapegin, Aleksei Chekmachev",
    author_email="dasapegin@sr-systems.ru",
    url="https://github.com/Membrizard/PyVaporation",
    download_url="https://github.com/Membrizard/PyVaporation/archive/refs/tags/v1.1.5.tar.gz",
    long_description=open("README.md", "r").read(),
    long_description_content_type="text/markdown",
    keywords=[
        "pervaporation",
        "membrane",
        "chemistry",
        "modelling",
        "chemical-engineering",
        "scientific"
    ],
    install_requires=[
        "joblib==1.1.0",
        "matplotlib==3.5.2",
        "pandas==1.3.5",
        "scipy==1.7.3",
        "attr==0.3.1",
        "attrs==21.4.0",
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
