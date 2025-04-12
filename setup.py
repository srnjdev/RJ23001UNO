# setup.py
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="RJ23001UNO",  
    version="1.0.0",
    author="Sergio Ramírez",
    author_email="rj23001@ues.edu.sv",
    description="Una librería para resolver sistemas de ecuaciones lineales y no lineales. Usando los métodos de: Eliminación de Gauss, Gauss-Jordan, Crammer, Jacobi, Gauss-Seidel, Bisección",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/srnjdev/RJ23001UNO",  
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
