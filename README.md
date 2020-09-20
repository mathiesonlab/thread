# thread

`thread` is a software package for reconstructing ancestral genomes, given a
pedigree structure and genotype data from recent individuals. Genotype or
sequence data should be phased (recommended: `SHAPEIT2`) and then IBD segments
should be called (recommended: `GERMLINE`).

1) Produce a `json` file (creates a dictionary to hold IBD information from GERMLINE):

```
$ python3 match2json.py -g example/toy_germline.match -s example/toy_pedigree.txt -m example/toy.map -p example/toy_geno.ped -j example/toy.json
```

2) Run the reconstruction pipeline, optionally producing a `ped` file of sequence information
for the reconstructed individuals:

```
$ python3 thread.py -g example/toy_germline.match -s example/toy_pedigree.txt -m example/toy.map -j example/toy.json -p example/recon.ped
```

This `recon.ped` file can be compared with the true ancestral haplotypes (known in
    this simulated) example: `toy_anc.ped`.

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
