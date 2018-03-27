oligo-melting v2.0.1.post2
===

A Python3 package for melting temperature calculation of oligonucleotides hybridization and secondary structures.

#### Features

* Handles DNA:DNA, RNA:RNA and DNA/RNA hybridizations.
* Corrects for both salt and chemicals.
* Compatible with UNAfold (OligoArrayAux) output for direct salt and chemicals (denaturants) correction.
* Produces melting curves for the provided sequence.
* Input either as single sequence or FASTA file.
* Slightly faster than BioPython.SeqUtils.MeltingTemp.
* Provides dG, dS and dH alongside melting temperature.

#### Limitations

* Does not handle mismatches or dangling ends
* Does not handle ambiguous bases

Installation
---

To **install**, run the following:

```
git clone http://github.com/ggirelli/gpseq-img-py
cd gpseq-img-py
sudo -H pip3 install .
```

To **uninstall** run the following from within the repository folder:

```
sudo -H pip3 uninstall oligo_melting
```

To **update**, first uninstall, and then run the following from within the repository folder.

```
git pull
sudo -H pip3 uninstall oligo_melting
sudo -H pip3 install .
```

Usage
---

### From command line

#### Duplexes

The `melt_duplex` command allows to calculate the melting temperature of a nucleic acid duplex, provided the sequence of one of the two strands.

The hybridization delta free energy calculation is based on the N-N thermodynamic values in literature and is available for DNA:DNA[3], RNA:RNA[1] and DNA:RNA[2] duplexes. The melting temperature calculation is based on Santalucia, 1998[4]. Sodium and cagnesium concentration correction is based on the work of Owczarzy et al[5][6]. Formamide correction can be performed based on two different published models[7][8].

* Use the `-t` option to specify the **type of nucleic acid duplex**.
* Use `--fa-mode` to switch between linear melting temperature **formamide-based correction**[7] and linear &Delta;G formamide-based correction[8].
* Use `--fa-mvalue` together with `--fa-mode wright` to specify the **m-value** for the formamide-based correction.
* Provide the path to a **fasta file** instead of a single sequence to calculate the melting temperature of every sequence in the file.
* Use the `-v` option to trigger the verbose mode, which provides more **details** for every single sequence.
* Use `-C` for the temperature in **degree Celsius** instead of Kelvin.
* Use `--out-curve` to specify a file where to save estimated single-sequence **melting curves** with temperature range and step around the melting temperature as defined with `--t-curve`.

#### Secondary structure

The `melt_secstr` script allows to correct the melting temperature of a nucleic acid secondary structure, previously calculated with OligoArrayAux, and to produce the corresponding melting curves.

### As a library

Import the package and use the corresponding functions.

```python3
import oligo_melting as OligoMelt

seq = "CAGTCAGTCGATC"

# Calculate melting temperature for 25uM oligos
(name, g, h, s, tm, seq) = OligoMelt.Duplex.calc_tm(seq, oligo_conc = 25e-6)
print(tm)

# Adjust for 300 mM [Na+]
tm = OligoMelt.Duplex.adj_ions(tm, 0.3, 0, seq)
print(tm)
```

The `Duplex` module contains functions for duplex hybridization and melting temperature calculation, while the `SecStr` module contains similar methods for evaluating secundary structure melting temperatures.

References
---

* [1]: Freier et al, PNAS(83), 1986;
* [2]: Sugimoto et al, Biochemistry(34), 1995.
* [3]: Allawi & Santalucia, Biochemistry(36), 1997;
* [4]: SantaLucia, PNAS(95), 1998;
* [5]: Owczarzy et al, Biochemistry(43), 2004;
* [6]: Owczarzy et al, Biochemistry(47), 2008;
* [7]: McConaughy et al, Biochemistry(8), 1969;
* [8]: Wright et al, Appl. env. microbiol.(80), 2014.

License
---

```
MIT License
Copyright (c) 2017 Gabriele Girelli
```

---

This project comes from the [potpourri](https://github.com/ggirelli/potpourri) sandbox.  \\( ﾟヮﾟ)/