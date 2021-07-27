from setuptools import setup, Extension

"""
A minimalist setup is shown.
"""
module = Extension("spkmeans",
                   sources=[
                       'spkmeans.c',
                       'matrix_op.c',
                       'spkmeansmodule.c'
                   ])
setup(name='spkmeans',
      version='1.0',
      description='Python wrapper for custom C extension',
      ext_modules=[module])