from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='feastruct',
      version='0.1',
      description='A python package for structural finite element analysis.',
      long_description=readme(),
      long_description_content_type='text/markdown',
      url='https://github.com/robbievanleeuwen/feastruct',
      author='Robbie van Leeuwen',
      author_email='robbie.vanleeuwen@gmail.com',
      license='MIT',
      packages=[
          'feastruct', 'feastruct.fea', 'feastruct.fea.elements', 'feastruct.post',
          'feastruct.pre', 'feastruct.solvers',
          'feastruct.examples', 'feastruct.tests.implementation',
          'feastruct.tests.mastan2_textbook'
      ],
      install_requires=[
          'numpy', 'scipy', 'matplotlib'
      ],
      include_package_data=True,
      zip_safe=False)
