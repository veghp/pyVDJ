import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='pyvdj',
    version='0.1.1',
    author='Peter Vegh',
    author_email='peter.vegh@newcastle.ac.uk',
    description='V(D)J sequencing data analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/veghp/pyVDJ',
    packages=setuptools.find_packages(),
    install_requires=[
        'pandas',
        'anndata',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: OS Independent',
    ],
)
