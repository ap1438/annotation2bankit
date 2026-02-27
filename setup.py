from setuptools import setup, find_packages

setup(
    name             = 'augustus_bankit_toolkit',
    version          = '1.0.0',
    description      = (
        'Convert AUGUSTUS gene predictions and GenBank/GFF annotations '
        'to NCBI BankIt-compatible submission files.'
    ),
    long_description = open('README.md').read(),
    long_description_content_type = 'text/markdown',
    author           = 'Akash Parida',
    author_email     = 'akash.parida@ruhr-uni-bochum.de',
    url              = 'https://github.com/yourusername/augustus_bankit_toolkit',
    license          = 'MIT',
    python_requires  = '>=3.6',
    packages         = find_packages(),
    py_modules       = ['augustus2bankit', 'genbank2bankit', 'validate_bankit'],
    entry_points     = {
        'console_scripts': [
            'augustus2bankit = augustus2bankit:main',
            'genbank2bankit  = genbank2bankit:main',
            'validate_bankit = validate_bankit:main',
        ],
    },
    install_requires = [],   # standard library only
    extras_require   = {
        'dev': ['pytest>=7.0'],
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Operating System :: OS Independent',
    ],
    keywords = (
        'bioinformatics AUGUSTUS GenBank GFF GFF3 GTF NCBI BankIt '
        'gene annotation submission feature-table genomics'
    ),
)
