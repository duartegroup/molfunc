from setuptools import setup
from Cython.Build import cythonize
from setuptools.extension import Extension

extensions = [Extension('molfunc_ext', ['molfunc/ext/molfunc_ext.pyx'])]

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
