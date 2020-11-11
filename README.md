# thread

`thread` is a software package for reconstructing ancestral genomes, given a
pedigree structure and genotype data from recent individuals. See our [preprint](https://www.biorxiv.org/content/10.1101/2020.01.15.908459v2) for more details about the method. Before using `thread`, genotype or
sequence data should be phased (recommended: `SHAPEIT2`) and then IBD segments
should be called (recommended: `GERMLINE`). Then `thread` can be applied using two steps:

1) Produce a `json` file (creates a dictionary to hold IBD information from `GERMLINE`):

```
$ python3 match2json.py -g example/toy_germline.match -s example/toy_pedigree.txt -m example/toy.map -p example/toy_geno.ped -j example/toy.json
```

2) Run the reconstruction pipeline, optionally producing a `ped` file of sequence information
for the reconstructed individuals:

```
$ python3 thread.py -g example/toy_germline.match -s example/toy_pedigree.txt -m example/toy.map -j example/toy.json -p example/recon.ped
```

This `recon.ped` file can be compared with the true ancestral haplotypes (known in
    this simulated example): `toy_anc.ped`.

We have implemented two different algorithms for selecting the source of each IBD segment ("min path" and "max prob"). "min path"
is the default. To use "max prob" instead, use the `-x` flag:

```
$ python3 thread.py -g example/toy_germline.match -s example/toy_pedigree.txt -m example/toy.map -j example/toy.json -p example/recon.ped -x
```

To create reproducible results, use the `PYTHONHASHSEED`. For example:

```
$ PYTHONHASHSEED=42 python3 thread.py -g example/toy_germline.match -s example/toy_pedigree.txt -m example/toy.map -j example/toy.json -p example/recon.ped
```

Note that this software is under active development. Please contact Sara
Mathieson (smathieson [at] haverford [dot] edu) with any questions.


Contributors:

* Kelly Finke `kellyfinke`
* Michael Kourakos `mkourak1` / `MikeyManiac`
* Gabriela Brown `gabcbrown`
* Sara Mathieson `saramathieson`
