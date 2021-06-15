from setuptools import setup

setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name='Measurements',
    url='https://github.com/jladan/package_demo',
    author='Tobia Diggelmann',
    author_email='dtobia@ethz.ch',
    # Needed to actually package something
    packages=['platepy'],
    # Needed for dependencies
    install_requires=['numpy', 'scipy', 'pandas', 'gmsh', 'matplotlib'],
    # *strongly* suggested for sharing
    version='0.1',
    # The license can be anything you like
    license='MIT',
    description='package to compute plates using finite elements',
    # We will also need a readme eventually (there will be a warning)
    # long_description=open('README.txt').read(),
)