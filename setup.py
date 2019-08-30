from distutils.core import setup
from setuptools import find_packages


setup(
    name="seqhighway",
    version="0.1.0",
    packages=find_packages(),
    install_requires=open("requirements.txt").readlines(),
    python_requires='>=3.5',
    long_description=open("readme.md").read(),
    description="HTML writer for annotated biological sequences.",
)
