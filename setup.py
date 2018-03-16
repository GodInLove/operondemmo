from setuptools import setup, find_packages

self_version = "0.0"

try:
    LONG_DESCRIPTION = open("README.rst", "rb").read().decode("utf-8")
except IOError:
    LONG_DESCRIPTION = "an independent demmo of KNOWN operon predict method"

setup(
    name='operondemmo',
    version=self_version,
    packages=find_packages(),
    url='https://github.com/GodInLove/operondemmo',
    license='GPL-3.0',
    author='yaodongliu',
    author_email='yd.liu.scu@gmail.com',
    description=LONG_DESCRIPTION,
)
