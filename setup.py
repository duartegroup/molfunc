from setuptools import setup
from Cython.Build import cythonize
from setuptools.extension import Extension

extensions = [Extension('esp_gen',
                        ['cgbind/ext/esp_gen.pyx'])]

setup(
    name='cgbind',
    version='1.0.0',
    packages=['ffunc'],
    ext_modules=cythonize(extensions, language_level="3", annotate=True),
    url='',
    license='MIT',
    author='Tom Young',
    author_email='tom.young@chem.ox.ac.uk',
    description='fast molecular functionalisation'
)
