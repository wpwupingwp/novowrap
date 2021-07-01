#!/usr/bin/python3

import setuptools

with open('README.md', 'r') as _:
    long_description = _.read()

with open('requirements.txt', 'r') as _:
    requires = [i.strip() for i in _.readlines()]

setuptools.setup(
    author='Ping Wu',
    author_email='wpwupingwp@outlook.com',
    description='Assembly Plastid Genome',
    install_requires=requires,
    license='GNU AGPL v3',
    long_description=long_description,
    long_description_content_type='text/markdown',
    name='novowrap',
    packages=setuptools.find_packages(),
    # f-string and Path-like
    python_requires='>=3.7',
    url='https://github.com/wpwupingwp/novowrap',
    version='1.13',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
)
