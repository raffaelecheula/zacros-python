import setuptools

with open('requirements.txt', 'r') as f:
    requirements = f.readlines()

with open('README.md', 'r') as f:
    readme = f.read()

setuptools.setup(
    name="zacros-python",
    version="0.1.1",
    url="https://github.com/raffaelecheula/zacros-python.git",

    author="Raffaele Cheula",
    author_email="cheula.raffaele@gmail.com",

    description="Python wrapper for Zacros KMC program.",
    long_description=readme,
    license='GPL-3.0',

    packages=[
        'zacros',
    ],
    package_dir={
        'zacros': 'zacros'
    },
    install_requires=requirements,
    python_requires='>=3.5, <4',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python',
    ],
)
