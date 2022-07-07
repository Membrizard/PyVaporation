from distutils.core import setup

setup(
    name="pyvaporation",
    packages=["pyvaporation"],
    version="0.1.6",
    license="Apache license 2.0",
    description="Set of tools for modelling membrane pervaporations",
    author="Denis Sapegin, Aleksei Chekmachev",
    author_email="a.checkmachev@gmail.com",
    url="https://github.com/Membrizard/PyVaporation",
    download_url="https://github.com/Membrizard/PyVaporation/archive/refs/tags/v0.1.6.tar.gz",
    long_description_content_type='text/markdown',
    keywords=[
        "pervaporation",
        "membrane",
        "chemistry",
        "modelling",
    ],
    install_requires=[
        'joblib>=1.1.0',
        'matplotlib>=3.5.2',
        'pandas>=1.4.2',
        'scipy>=1.8.0',
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Build Tools",
        "License :: OSI Approved :: Apache license 2.0",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
)
