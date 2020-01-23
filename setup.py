from distutils.core import setup

setup(
        name = 'odnpLab',
        packages = ['odnpLab'],
        version = '1.0',
        license = 'MIT',
        description = 'odnpLab - A NMR Processing Library for ODNP Experiments',
        author = 'Timothy Keller',
        author_email = 'tkeller@bridge12.com',
        url = '',
        download_url = '',
        keywords = ['ODNP','NMR'],
        install_requires = [
            'numpy',
            'scipy',
            'matplotlib',
            'h5py',
            ]
        classifiers = [
            'Development Status :: 3 - Alpha',
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 3',
            ]

        )
