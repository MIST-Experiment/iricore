from numpy.distutils.core import Extension
from numpy.distutils.core import setup
from setuptools import find_packages
import os

fsources = ['iri_fcore.pyf', 'iri_fcore.f90', 'cira.for', 'igrf.for', 'iridreg.for', 'iriflip.for',
            'irifun.for', 'irisub.for', 'iritec.for']

here = os.path.dirname(os.path.realpath(__file__))

iri_fcore = Extension(
    'iricore.iri_fcore',
    sources=[os.path.join(here, 'src/iricore/src/') + fs for fs in fsources],
    extra_f90_compile_args=['-fPIC', '-O3', '-w', '-Wtabs'],
    extra_f77_compile_args=['-fPIC', '-O3', '-w', '-Wtabs'],
    f2py_options=['--static'],
)

if __name__ == "__main__":
    setup(name='iricore',
          version="1.0.1",
          author="lap1dem",
          author_email="vadym.bidula@gmail.com",
          url='https://github.com/lap1dem/iricore',
          license='LICENSE',
          description='A python interface to iri2016 Fortran code',
          long_description=open('README.md').read(),
          long_description_content_type="text/markdown",
          install_requires=[
              "numpy >= 1.21",
          ],
          package_dir={"": "src"},
          packages=find_packages(where="src"),
          ext_modules=[iri_fcore],
          # include_package_data=True,
          package_data={
              "iricore": [
              "data/ccir/*.asc",
              "data/igrf/*.dat",
              "data/index/*.dat",
              "data/mcsat/*.dat",
              "data/ursi/*.asc",
              ],
          },
          python_requires=">=3.7, <3.10",
          classifiers=[
              # How mature is this project? Common values are
              #   3 - Alpha
              #   4 - Beta
              #   5 - Production/Stable
              "Development Status :: 4 - Beta",
              "License :: OSI Approved :: MIT License",
              "Programming Language :: Python :: 3",
              "Programming Language :: Python :: 3.7",
              "Programming Language :: Python :: 3.8",
              "Programming Language :: Python :: 3.9",
          ],
          )