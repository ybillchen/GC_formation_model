# GC_formation_model

[![version](https://img.shields.io/badge/version-0.2-blue.svg)](https://github.com/ybillchen/GC_formation_model)
[![license](https://img.shields.io/github/license/ybillchen/GC_formation_model)](LICENSE)
[![workflows](https://img.shields.io/github/actions/workflow/status/ybillchen/GC_formation_model/build.yaml?logo=github)](https://github.com/ybillchen/GC_formation_model/actions/workflows/build.yaml)

A post-processing model for globular cluster formation in cosmological simulations.

The code is open source under a [BSD 3-Clause License](LICENSE), which allows you to redistribute and modify the code with moderate limitations. If you use this code for a publication, we kindly request you to cite the following original papers.

- [A. L. Muratov \& O. Y. Gnedin (2010)](https://ui.adsabs.harvard.edu/abs/2010ApJ...718.1266M/abstract), ApJ, **718**, 1266
- [H. Li \& O. Y. Gnedin (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJ...796...10L/abstract), ApJ, **796**, 10
- [N. Choksi, O. Y. Gnedin, \& H. Li (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.2343C/abstract), MNRAS, **480**, 2343
- [Y. Chen \& O. Y. Gnedin (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.514.4736C/abstract), MNRAS, **514**, 4736
- [Y. Chen \& O. Y. Gnedin (2023)](https://ui.adsabs.harvard.edu/abs/2023MNRAS.522.5638C/abstract), MNRAS, **522**, 5638
- [Y. Chen \& O. Y. Gnedin (2024a)](https://ui.adsabs.harvard.edu/abs/2024MNRAS.527.3692C/abstract), MNRAS, **527**, 3692
- [Y. Chen \& O. Y. Gnedin (2024b)](https://ui.adsabs.harvard.edu/abs/2024OJAp....7E..23C/abstract), OJAp, **7**, 23

We also provide a toolkit for parallelization. If you are planning to use the toolkit, please contact us for access!

## Background

This is a long developing project starting from [Muratov \& Gnedin (2010)](https://ui.adsabs.harvard.edu/abs/2010ApJ...718.1266M/abstract). Some important milestones are [Li \& Gnedin (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJ...796...10L/abstract), [Choksi, Gnedin, \& Li (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.2343C/abstract), [Chen \& Gnedin (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.514.4736C/abstract), and [Chen \& Gnedin (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv230108218C/abstract). This version is the most recent model by [Chen \& Gnedin (2024a)](https://ui.adsabs.harvard.edu/abs/2024MNRAS.527.3692C/abstract). 

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

### Pull request protocol

We recommend you to contribute code to `GC_formation_model` following [GitHub flow](https://docs.github.com/en/get-started/quickstart/github-flow). To summarize, you submit a pull request via the following steps:

1. Clone the repository.
2. Create and checkout a new branch. For example, a new branch called `new_feature`.
3. Make changes on `new_feature` and never touch the `main` branch again until you are ready to merge.
4. When you feel ready, submit a pull request on GitHub.
5. There may be conflicts. If so, you need to 
	1. Checkout the `main` branch and pull from `origin`.
	2. Rebase `new_feature` on `main` and address the conflicts (recommended).
	3. Alternatively, you can compare `new_feature` with `main` and fix all conflicts.
	4. Your pull request will update automatically.
6. If your pull request is approved, we will squash and merge your commits. 
7. We will delete `new_feature` on GitHub when it's merged. You can choose to delete it loacally as well. 

**_NOTE:_** Any slight modification may entirely change the random number generation! To keep repeatability of the model, please construct a new random generator for the need of new random numbers

## Maintainers

- [@Yingtian (Bill) Chen](https://github.com/ybillchen)
- [@Oleg Gnedin](https://github.com/ognedin)
