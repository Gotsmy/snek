import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="snek",
    version="0.01",
    author="Mathias Gotsmy",
    author_email="mathias.gotsmy@univie.ac.at",
    description="Convenience functions for COBRApy.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Gotsmy/snek",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNUv3 License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "."},
    packages=setuptools.find_packages(where="."),
    python_requires=">=3.7",
#    include_package_data=True,
#    package_data={'': ['data/*.csv']},
)
