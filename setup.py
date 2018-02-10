import sys
from setuptools import setup

if sys.version_info < (2, 6):
    raise Exception('fastools requires Python 2.6 or higher.')

requires = ['matplotlib', 'pybedtools', 'scipy', 'wiggelen']

# Python 2.6 does not include the argparse module.
try:
    import argparse
except ImportError:
    requires.append('argparse')

import snvmixtools as distmeta

setup(
    name='snvmixtools',
    version=distmeta.__version__,
    description='SNVMix analysis toolkit.',
    long_description=distmeta.__doc__,
    author=distmeta.__author__,
    author_email=distmeta.__contact__,
    url=distmeta.__homepage__,
    license='MIT License',
    platforms=['any'],
    packages=['snvmixtools'],
    install_requires=requires,
    entry_points = {
        'console_scripts': [
            'snvmixtools = snvmixtools.cli:main'
        ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
    ],
    keywords='bioinformatics'
)
