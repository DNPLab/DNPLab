import setuptools
from distutils.core import setup

setup(
        name = 'dnpLab',
        packages = ['dnpLab'],
        version = '1.0',
        license = 'MIT',
        description = 'dnpLab - A NMR Processing Library for ODNP Experiments',
        author = 'Timothy Keller',
        author_email = 'tkeller@bridge12.com',
        url = '',
        download_url = '',
        keywords = ['ODNP','DNP','NMR'],
        install_requires = [
            'numpy',
            'scipy',
            'matplotlib',
            'h5py',
            ],
        classifiers = [
            'Development Status :: 3 - Alpha',
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 3',
            ]

        )
