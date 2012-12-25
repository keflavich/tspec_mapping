
from distutils.core import setup

with open('README') as file:
    long_description = file.read()

setup(name='tspec_mapping',
      version='0.1',
      description='TripleSpec slit mapping',
      author='Adam Ginsburg',
      author_email='adam.g.ginsburg@gmail.com',
      url='https://github.com/keflavich/tspec_mapping',
      packages=['tspec_mapping'], 
     )
