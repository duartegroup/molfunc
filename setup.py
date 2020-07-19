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
      version='1.0.0b0',
      packages=['molfunc'],
      package_data={'': ['fragments_lib/*']},
      license='MIT',
      author='Tom Young',
      url='https://github.com/duartegroup/molfunc',
      entry_points={'console_scripts': ['molfunc = molfunc.molfunc:main']},
      ext_modules=cythonize(extensions, language_level="3"),
      author_email='tom.young@chem.ox.ac.uk',
      description='Fast molecular functionalisation',
      platforms='any',
      long_description='molfunc enables adding molecular fragments e.g. '
                       'Me, Ph etc. to existing 3D structures.')
