from setuptools import setup
from Cython.Build import cythonize
from setuptools.extension import Extension

extensions = [Extension('molfunc_ext',
                        [f'molfunc/molfunc_ext.pyx'],
                        include_dirs=['molfunc/include'],
                        language='c++',
                        compiler_directives={'language_level': '3'},
                        extra_compile_args=["-std=c++17", "-Wno-missing-braces"],
                        extra_link_args=["-std=c++17"]
                        )]

setup(name='molfunc',
      version='2.0.0',
      packages=['molfunc'],
      package_data={'': ['src/species/data/*']},
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
