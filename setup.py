# Licensed under BSD-3-Clause License - see LICENSE

from setuptools import setup, find_packages

with open('README.md', 'r') as f:
    long_description = f.read()

version = {}
with open('GC_formation_model/version.py') as fp:
    exec(fp.read(), version)

setup(
    name = 'GC_formation_model',
    packages = find_packages(),
    version = version['__version__'],
    url = 'https://github.com/ybillchen/GC_formation_model',
    license = 'BSD-3-Clause',
    author = 'Bill Chen',
    author_email = 'ybchen@umich.edu',
    description = 'A post-processing model for globular cluster formation in cosmological simulations.',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    install_requires = ['numpy', 'scipy', 'h5py'],
    python_requires = '>=3.8',
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
    ],
)
