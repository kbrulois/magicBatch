import os
import sys
import shutil
from warnings import warn
from distutils.core import setup


setup(name='magicBatch',
      version='1.0.4',
      description='MAGICbatch',
      author='Kevin Brulois',
      author_email='kevin.brulois@gmail.com',
      package_dir={'': 'src'},
      packages=['magicBatch'],
      url='https://github.com/kbrulois/magicBatch',
      download_url='https://github.com/kbrulois/magicBatch/archive/refs/tags/v1.0.0.tar.gz',
      install_requires=[
          'numpy>=1.10.0',
          'scipy>=0.14.0',
          'seaborn',
          'scikit-learn',
          'networkx',
          'statsmodels',
          'datatable'
      ],
      scripts=['src/magicBatch/MAGIC.py'],
      )


# get location of setup.py
setup_dir = os.path.dirname(os.path.realpath(__file__))
