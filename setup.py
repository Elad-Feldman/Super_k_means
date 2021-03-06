from setuptools import setup, Extension

##   python setup.py build_ext --inplace

"""
A minimalist setup is shown.
"""
module = Extension("spkmeans",
                   sources=[
                       'spkmeans.c',
                       'spkmeansmodule.c'
                   ])
setup(name='spkmeans',
      version='1.0',
      description='Python wrapper for custom C extension',
      ext_modules=[module])