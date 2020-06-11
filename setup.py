from setuptools import setup

setup(name='molfunc',
      version='1.0.0a0',
      packages=['molfunc'],
      license='MIT',
      author='Tom Young',
      url='https://github.com/duartegroup/molfunc',
      entry_points={'console_scripts': ['molfunc = molfunc.molfunc:main']},
      author_email='tom.young@chem.ox.ac.uk',
      description='Fast molecular functionalisation')
