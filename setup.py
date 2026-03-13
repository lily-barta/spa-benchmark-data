from setuptools import setup, find_packages

setup(
    name="spa-benchmark",
    version="0.1",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "tequila-basic==1.9.10.dev0",
        "pandas",
        "matplotlib", 
        "qulacs",
        "pyscf",
        "openfermion"
    ],
    url="https://github.com/lily-barta/spa-benchmark",
    author="Lily Barta",
)
