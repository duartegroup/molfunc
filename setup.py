from setuptools import setup
from Cython.Build import cythonize
from setuptools.extension import Extension
import os

# Try to build the extension from the Cython generated C
if os.path.exists('molfunc/ext/molfunc_ext.c'):
    ext = 'c'
else:
    ext = 'pyx'

extensions = [Extension('molfunc_ext', [f'molfunc/ext/molfunc_ext.{ext}'])]

setup(name='molfunc',
      version='1.0.0a3',
      packages=['molfunc'],
      license='MIT',
      author='Tom Young',
      url='https://github.com/duartegroup/molfunc',
      entry_points={'console_scripts': ['molfunc = molfunc.molfunc:main']},
      ext_modules=cythonize(extensions, language_level="3"),
      author_email='tom.young@chem.ox.ac.uk',
      description='Fast molecular functionalisation')
