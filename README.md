# Tacopy

Taxonomic Compositions in Python

## Installation

```bash
pip install tacopy
```

## Usage

Read in data from tidytacos format (e.g. a folder with a table for the counts, a table for the samples and a table for the taxa.): 

```python
import tacopy as tp

tt = tp.read_tacos('tests/data/urt')
ttg = tp.aggregate_taxa(tt, rank="genus")

```