# GC_formation_model

[![version](https://img.shields.io/badge/version-0.1-blue.svg)](https://github.com/ybillchen/GC_formation_model)
[![license](https://img.shields.io/github/license/ybillchen/GC_formation_model)](LICENSE)
[![workflows](https://img.shields.io/github/actions/workflow/status/ybillchen/GC_formation_model/build.yaml?logo=github)](https://github.com/ybillchen/GC_formation_model/actions/workflows/build.yaml)

A post-processing model for globular cluster formation in cosmological simulations.

The code is open source under a [BSD 3-Clause License](LICENSE), which allows you to redistribute and modify the code with moderate limitations. If you use this code for a publication, we kindly request you to cite the following original papers.

- [A. L. Muratov \& O. Y. Gnedin (2010)](https://ui.adsabs.harvard.edu/abs/2010ApJ...718.1266M/abstract), ApJ, **718**, 1266
- [H. Li \& O. Y. Gnedin (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJ...796...10L/abstract), ApJ, **796**, 10
- [N. Choksi, O. Y. Gnedin, \& H. Li (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.2343C/abstract), MNRAS, **480**, 2343
- [Y. Chen \& O. Y. Gnedin (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.514.4736C/abstract), MNRAS, **514**, 4736
- [Y. Chen \& O. Y. Gnedin (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv230108218C/abstract), MNRAS submitted (arXiv:2301.08218)

## Background

This is a long developing project starting from [Muratov \& Gnedin (2010)](https://ui.adsabs.harvard.edu/abs/2010ApJ...718.1266M/abstract). Some important milestones are [Li \& Gnedin (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJ...796...10L/abstract), [Choksi, Gnedin, \& Li (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.2343C/abstract), and [Chen \& Gnedin (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.514.4736C/abstract). This version is the most recent model by [Chen \& Gnedin (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv230108218C/abstract). 

## Install

We have tested `GC_formation_model` on `python >= 3.8`. However, lower versions may also work. The prerequisites of this package are
```
numpy
scipy
h5py
```

To download the packge, `git clone` the source code from [GitHub](https://github.com/ybillchen/GC_formation_model):
```shell
$ git clone https://github.com/ybillchen/GC_formation_model.git
```
Next, `cd` the folder and use `pip` to install it:
```shell
$ cd GC_formation_model/
$ pip install -e .
```
The `-e` command allows you to make changes to the code.

## Usage

To start with, let's run the model with default parameters
```python
>>> from GC_formation_model import run
>>> from params_example import params
>>> run(params)
```
You may want to use your own paramters. Then simply replace `params_example` with the name of your paramter file.


## Contribute

Feel free to dive in! [Raise an issue](https://github.com/ybillchen/GC_formation_model/issues/new) or submit pull requests.

## Maintainers

- [@Yingtian (Bill) Chen](https://github.com/ybillchen)
- [@Oleg Gnedin](https://github.com/ognedin)
