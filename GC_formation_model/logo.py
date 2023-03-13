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