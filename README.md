# TensorNetworkManyBodyCavityQED


Tensor Network methods for simulating spin-1/2 system interacting via a single-mode (possibly leaky) cavity mode.

The tensor network codes are written in C++, with taylored functions writte on top of the ITensor v3 library for tensor network methods.

Exact diagonalization calculation are written in Python, based on the quimb Python package. These codes are mostly used as benchmarks for tensor newtork codes.

## Description

Contains codes for simulating one-dimensional many-body cavity-QED systems.

- rydberg_in_cavity: contains codes for simulating 1d array of Rydberg atoms interacting with a (also leaky) bosonic mode (i.e. cavity field)

- collective_light_matter: contains codes for simulating paradigmatic collective photon-matter models, namely the Dicke and Tavis-Cummings model with tensor-networks.

- library_cpp: contains functions needed for tensor network

- library_python: contains methods useful for Exact Diagonalization calculation.

## Prerequisites

C++ ITensor v3 library: You need the ITensor v3 library to be installed. See the ITensor Installation Guide for details https://itensor.org/docs.cgi?vers=cppv3&page=install

C++17: This library requires a C++17-compliant compiler.

quimb: needed the Python package quimb for Exact Diagonalization (written in Python). Install it via the command
```bash
pip install quimb
```

## Installation
Clone the repository:
```bash
git clone https://github.com/riccardovalencia/TensorNetworkManyBodyCavityQED.git
```

## Usage
Modify path to your ITensor library in order to compile via 'make' command the .cpp codes.

## License

[MIT](https://choosealicense.com/licenses/mit/)
