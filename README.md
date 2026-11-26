# insertsize

[![GitHub release](https://img.shields.io/github/v/release/grenaud/insertsize)](https://github.com/grenaud/insertsize/releases/latest)
[![License](https://img.shields.io/github/license/grenaud/insertsize)](LICENSE)
[![How to cite](https://img.shields.io/badge/how%20to-cite-blue.svg)](#how-to-cite)

This program produces the insert size for a BAM file.

To build the software simply tyoe:

    cd src/
    make

To check the build against the small test set in `test/`:

    cd src/
    make test

See [test/README.md](test/README.md) for what the test data covers.

## How to cite

insertsize has no accompanying paper. If you use it, please cite this repository:

> Renaud G. insertsize. https://github.com/grenaud/insertsize

BibTeX:

```bibtex
@misc{Renaud_insertsize,
  author       = {Renaud, Gabriel},
  title        = {{insertsize}: insert size from {BAM} files},
  year         = {2026},
  howpublished = {\url{https://github.com/grenaud/insertsize}}
}
```

## License

Copyright Gabriel Renaud. insertsize is free software, released under the GNU General Public License, version 3 or later (GPL-3.0-or-later). See [LICENSE](LICENSE).
