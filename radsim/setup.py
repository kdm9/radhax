from setuptools import setup
from setuptools.command.test import test as TestCommand
import versioneer


# Inspired by the example at https://pytest.org/latest/goodpractises.html
class NoseCommand(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        # Run nose ensuring that argv simulates running nosetests directly
        import nose
        nose.run_exit(argv=['nosetests'])

versioneer.VCS = 'git'
versioneer.versionfile_source = 'radsim/_version.py'
versioneer.versionfile_build = 'radsim/_version.py'
versioneer.tag_prefix = 'radsim-'
versioneer.parentdir_prefix = 'radsim-'

desc = """
RADsim: simulate reduced-represention sequencing platforms based on restriction
digestion, including RADseq and GBS.
"""

setup_requires = [
    'nose>=1.3',
]

install_requires = [
    'biopython==1.65',
    'docopt>=0.6',
    'screed>=0.8',
    'six>=1.9',
    'progressbar2>=2.7',
]

test_requires = [
    'pep8>=1.6,<1.7',
    'pylint>=1.4,<1.5',
    'coverage>=3.7',
]

command_classes=versioneer.get_cmdclass()
command_classes['test'] =  NoseCommand

setup(
    name="radsim",
    packages=['radsim', ],
    scripts=[
    ],
    version=versioneer.get_version(),
    cmdclass=command_classes,
    install_requires=install_requires,
    tests_require=test_requires,
    setup_requires=setup_requires,
    description=desc,
    author="Kevin Murray",
    author_email="spam@kdmurray.id.au",
    url="https://github.com/kdmurray91/radhax",
    keywords=[
        "bioinformatics",
        "SRA",
        "Next-gen Sequencing",
    ],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later " +
            "(GPLv3+)",
    ],
)
