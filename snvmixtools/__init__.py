"""
SNVMixtools: Various tools for the analysis and manipulation of SNVMix files.


Copyright (c) 2013 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2013 Jeroen F.J. Laros <j.f.j.laros@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""

__version_info__ = ('0', '0', '1')


__version__ = '.'.join(__version_info__)
__author__ = 'LUMC, Jeroen F.J. Laros'
__contact__ = 'j.f.j.laros@lumc.nl'
__homepage__ = 'https://git.lumc.nl/j.f.j.laros/snvmixtools'

usage = __doc__.split("\n\n\n")

def docSplit(func):
    return func.__doc__.split("\n\n")[0]

def version(name):
    return "%s version %s\n\nAuthor   : %s <%s>\nHomepage : %s" % (name,
        __version__, __author__, __contact__, __homepage__)
