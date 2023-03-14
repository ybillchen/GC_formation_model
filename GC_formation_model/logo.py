# Licensed under BSD-3-Clause License - see LICENSE

from packaging.version import Version

from .version import __version__

def print_logo():
    print(r'    ___    ___                         _        _ ')
    print(r'   / _ \  / __\  _ __ ___    ___    __| |  ___ | |')
    print(r'  / /_\/ / /    | `_ ` _ \  / _ \  / _` | / _ \| |')
    print(r' / /_\\ / /___  | | | | | || (_) || (_| ||  __/| |')
    print(r' \____/ \____/  |_| |_| |_| \___/  \__,_| \___||_|')
    print(r'                                                  ')

def print_version():
    print('version: %s'%Version(__version__))

def print_papers():
    print(r' - A. L. Muratov & O. Y. Gnedin (2010), ApJ, 718, 1266 (arXiv:1002.1325)')
    print(r' - H. Li & O. Y. Gnedin (2014), ApJ, 796, 10 (arXiv:1405.0763)')
    print(r' - N. Choksi, O. Y. Gnedin, & H. Li (2018), MNRAS, 480, 2343 (arXiv:1801.03515)')
    print(r' - Y. Chen & O. Y. Gnedin (2022), MNRAS, 514, 4736 (arXiv:2203.00599)')
    print(r' - Y. Chen & O. Y. Gnedin (2023), MNRAS submitted (arXiv:2301.08218)')