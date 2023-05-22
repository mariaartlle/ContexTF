from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='ContexTF',
    version='1.0',
    description='ContexTF: functional annotation of prokaryotic regulatory elements',
    author ='Maria Artigues-Lleixà',
    author_email='maria.artigues01@estudiant.upf.edu',
    url='https://github.com/mariaartlle/ContexTF',
    packages=find_packages(),
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
    'Development Status :: 5 - Production/Stable',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    entry_points={  
    'console_scripts': 'ContexTF=ContexTF.ContexTF:main'}
)
