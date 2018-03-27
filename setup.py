"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
import os

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
	long_description = f.read()

setup(name='oligo_melting',
	version='2.0.1.post2',
	description='''A Python3 package for melting temperature calculation of
		oligonucleotides hybridization and secondary structures.''',
	long_description=long_description,
	long_description_content_type='text/markdown',
	url='https://github.com/ggirelli/oligo-melting',
	author='Gabriele Girelli',
	author_email='gabriele.girelli@scilifelab.se',
	license='MIT',
	classifiers=[
		'Development Status :: 5 - Production/Stable',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 3 :: Only',
	],
	keywords='DNA chemistry melting temperature modeling RNA salt denaturant',
	packages=["oligo_melting"],
	install_requires=[],
	scripts=["bin/melt_duplex", "bin/melt_secstr"],
	test_suite="nose.collector",
	tests_require=["nose"],
)
