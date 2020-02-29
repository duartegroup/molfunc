from setuptools import setup

setup(
    name='ffunc',
    version='1.0.0',
    packages=['ffunc'],
    license='MIT',
    author='Tom Young',
    url='https://github.com/duartegroup/ffunc',
    entry_points={'console_scripts': ['ffunc = ffunc.ffunc:main']},
    author_email='tom.young@chem.ox.ac.uk',
    description='fast molecular functionalisation'
)
