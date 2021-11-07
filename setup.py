import sys
import os
from setuptools import setup
from Cython.Build import cythonize
from setuptools.extension import Extension
from setuptools.command.install import install as _install

here = os.path.dirname(os.path.abspath(__file__))


def _post_install():
    from subprocess import call
    call([sys.executable,
          'molfunc/scripts/compile_fragments.py',
          os.path.join(here, 'molfunc')],
         cwd=here)


class Install(_install):
    def run(self):
        self.execute(_post_install,
                     args=(),
                     msg="Compiling fragments into source")
        _install.run(self)


extensions = [Extension('molfunc_ext',
                        [f'molfunc/molfunc_ext.pyx'],
                        include_dirs=['molfunc/include'],
                        language='c++',
                        extra_compile_args=["-std=c++17", "-Wno-missing-braces", "-O3"],
                        extra_link_args=["-std=c++17"])
              ]

setup(name='molfunc',
      version='2.0.0b0',
      packages=['molfunc'],
      license='MIT',
      package_data={"molfunc": ["include/*.h", "src/*.cpp", "include/species/*.h",
                                "src/species/*.cpp", "scripts/*.py", "*.pyx"]},
      build_requires=["Cython"],
      author='Tom Young',
      url='https://github.com/duartegroup/molfunc',
      cmdclass={'install': Install},
      entry_points={'console_scripts': ['molfunc = molfunc.molfunc:main']},
      ext_modules=cythonize(extensions, language_level="3"),
      author_email='tom.young@chem.ox.ac.uk',
      description='Fast molecular functionalisation',
      platforms='any',
      long_description='molfunc enables adding molecular fragments e.g. '
                       'Me, Ph etc. to existing 3D structures.')
