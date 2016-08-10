"""
Prediction of protein beta-beta contacts at the residue level using direct coupling patterns

bbcontacts is a Python program predicting residue-level contacts between beta-strands by detecting patterns 
in matrices of predicted couplings. 
bbcontacts can make use of a secondary structure assignment or a secondary structure prediction.
"""

from setuptools import setup

# ===================================
DOCSTRING = __doc__.replace("\n", "")
VERSION = 1.0
# ===================================

setup(
    name='bbcontacts',
    description=DOCSTRING,
    version=VERSION,
    author='Jessica Andreani & Johannes Soeding',
    author_email='jessica.andreani@mpibpc.mpg.de',
    license='GNU Affero General Public License v3',
    url='https://github.com/fsimkovic/bbcontacts',
    package_dir={'bbcontacts': 'bbcontacts'},
    packages=['bbcontacts', 'bbcontacts.hmm'],
    package_data={'bbcontacts': ['bbcontacts.conf', 'paramfiles/*']},
    platforms=['Linux', 'Mac OS-X', 'Unix', 'Windows'],
    install_requires=['numpy', 'scipy'],
    entry_points={'console_scripts': ['bbcontacts = bbcontacts.main:main']},
    include_package_data=True,
    zip_safe=False,
)

