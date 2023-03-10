# GC-formation-model

[![version](https://img.shields.io/badge/version-v0.0-blue.svg?style=flat)](https://github.com/EnthalpyBill/GC-formation-model)
[![license](https://img.shields.io/github/license/EnthalpyBill/GC-formation-model?style=flat)](LICENSE)
[![workflows](https://img.shields.io/github/actions/workflow/status/ybillchen/GC-formation-model/build.yaml?logo=github)](https://github.com/ybillchen/GC-formation-model/actions/workflows/build.yaml)

A post-processing model for globular cluster formation in cosmological simulations.

If you use this code for a publication, we kindly request you to acknowledge the following original papers.

- [Muratov \& Gnedin (2010)](https://ui.adsabs.harvard.edu/abs/2010ApJ...718.1266M/abstract)
- [Li \& Gnedin (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJ...796...10L/abstract)
- [Choksi, Gnedin, \& Li (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.2343C/abstract)
- [Chen \& Gnedin (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.514.4736C/abstract)
- [Chen \& Gnedin (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv230108218C/abstract)

## Table of Contents

- [Background](#background)
- [Install](#install)
- [Usage](#usage)
- [Maintainers](#maintainers)
- [License](#license)

## Background

This is a long developing project starting from [Muratov \& Gnedin (2010)](https://ui.adsabs.harvard.edu/abs/2010ApJ...718.1266M/abstract). Some important milestones are [Li \& Gnedin (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJ...796...10L/abstract), [Choksi, Gnedin, \& Li (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.2343C/abstract), and [Chen \& Gnedin (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.514.4736C/abstract). This version is the most recent model by [Chen \& Gnedin (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv230108218C/abstract). 

## Install

you can `git clone` the source package from [GitHub](https://github.com/ybillchen/GC-formation-model):
```shell
$ git clone https://github.com/ybillchen/GC-formation-model.git
```
To build and install `GC_formation_model`, `cd` the folder and `pip install` it:
```shell
$ cd GC_formation_model/
$ pip install -e .
```
The `-e` command allows you to make changes to the code.

## Usage

To use the package, just import it as
```python
>>> import GC_formation_model as gm
```
To start with, let's run the model with default settings
```python
>>> gm.run()
```

## Maintainers

- [@Bill Chen](https://github.com/ybillchen)
- [@Oleg Gnedin](https://github.com/ognedin)

## License

[BSD 3-Clause License](LICENSE) &copy; Bill Chen
