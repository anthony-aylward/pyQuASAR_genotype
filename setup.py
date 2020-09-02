import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='pyQuASAR_genotype',
    version='0.1.7',
    author='Anthony Aylward',
    author_email='aaylward@eng.ucsd.edu',
    description='Command line tool for genotyping with QuASAR',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/anthony-aylward/pyQuASAR_genotype.git',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
    install_requires=['pyQuASAR', 'seqalign', 'wasp_map'],
    entry_points={
        'console_scripts': ['pyQuASAR-genotype=pyQuASAR_genotype.genotype:main']
    }
)
