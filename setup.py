
from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='fieldgenerator',
    version='1.0',
    author='Patrick Melix',
    author_email='chemistry@melix.me',
    url='https://github.com/patrickmelix/fieldgenerator',
    download_url='https://github.com/patrickmelix/fieldgenerator',
    license='LGPLv3',
    description='Python Library to create DL_POLY FIELD files and update them with new values',
    long_description=readme(),
    classifiers=[
            'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.5',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    keywords=[],
    install_requires=[
        'numpy',
        'ase',
        'simpleeval'
    ],
    packages=['fieldgenerator'],
    include_package_data=True,
    scripts=['bin/fieldgenerator']
)
