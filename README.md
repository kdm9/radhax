RADHAX
======

Scripts for analysis of reduced-representation sequencing (GBS, RADseq)

`radsim`
========

Simulates GBS/RAD digests.

### Installation

    pip install radsim

### Usage

There are three commands in radsim


##### ``radsim-digest``

Digitally digests the reference genome, returning GBS fragments.

##### ``radsim-rebed``

Returns a BED file of restriction sites. Optionally, one can return a window
around each RE site, to simulate digest-then-random-fragment libraries like the
original RADseq (`--length`).

