from setuptools import setup, Extension

"""
A minimalist setup is shown.
"""
module = Extension("spkmeans",
                   sources=[
                       'spkmeansmodule.c'
                   ])
setup(name='spkmeans',
      version='1.0',
      description='Python wrapper for custom C extension',
      ext_modules=[module])