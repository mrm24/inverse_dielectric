# Spglib

```{toctree}
---
maxdepth: 1
hidden:
---
   How to install spglib C-API <install>
   C-APIs <api>
   Spglib dataset <dataset>
   Magnetic spglib dataset (Experimental) <magnetic_dataset>
   Variables <variable>
   Spglib for Python <python-interface>
   Rust, Fortran, and Ruby interfaces <interface>
   Definitions and conventions <definition>
   References <references>
   Summary of releases <releases>
   API Reference <api/index>
   Developer Documentation <development/develop>
```

Spglib is a library for finding and handling crystal symmetries
written in C. There are the [python](python-interface.md),
{ref}`fortran, rust, and ruby <interfaces>` interfaces.

The detailed algorithm of spglib is described the following text:

- "$\texttt{Spglib}$: a software library for crystal symmetry search", Atsushi
  Togo, Kohei Shinohara, and Isao Tanaka, Sci. Technol. Adv. Mater., Meth. 4,
  2384822-2384836 (2024)
  \[[link-sg](https://doi.org/10.1080/27660400.2024.2384822)\]
- "Algorithms for magnetic symmetry operation search and identification of
  magnetic space group from magnetic crystal structure", Kohei Shinohara,
  Atsushi Togo, and Isao Tanaka, Acta Cryst. A 79, 390-398 (2023)
  \[[link-msg](https://doi.org/10.1107/S2053273323005016)\]

## Requirements

- C standard: As of version 2.1.0, Spglib explicitly enforces a minimum standard of
  [C11](https://en.cppreference.com/w/c/11)

## Features

- Find symmetry operations
- Identify space-group type
- Wyckoff position assignment
- Refine crystal structure
- Find a primitive cell
- Search irreducible k-points

## For more information

- Repository: <https://github.com/spglib/spglib>
- License: BSD-3-Clause

## Links

- Documentation (this page): <https://spglib.readthedocs.io/>
- Repository <https://github.com/spglib/spglib>
- Conda: <https://anaconda.org/conda-forge/spglib>
- PyPI: <https://pypi.org/project/spglib/>

## Contributors

Spglib thrives from the contributions of the community, and we would like to thank each and every one of [the contributors][contributors] as well as everyone participating in the [issue] and [discussion] forums.

The current main developers of Spglib are:

- {user}`Atsushi Togo <atztogo>` (National Institute for Materials Science)
- {user}`Kohei Shinohara <lan496>`
- {user}`Cristian Le <LecrisUT>`

## Acknowledgments

Spglib development is supported by

- National Institute for Materials Science

Spglib project acknowledges Paweł T. Jochym for deployment and packaging, Yusuke
Seto for the Crystallographic database, Jingheng Fu for layer group
implementation, Juan Rodriguez-Carvajal for the Magnetic space-group database,
Dimitar Pashov for the Fortran interface, and many other contributors.

[contributors]: https://github.com/spglib/spglib/graphs/contributors
[discussion]: https://github.com/spglib/spglib/discussions
[issue]: https://github.com/spglib/spglib/issues
