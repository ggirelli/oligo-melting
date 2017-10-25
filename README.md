oligo-melting
===

The `oligomelt.py` script, implemented in Python, allows to calculate the melting temperature of a nucleic acid duplex, provided the sequence of one of the two strands.

The hybridization delta free energy calculation is based on the N-N thermodynamic values in literature and is available for DNA:DNA[^3], RNA:RNA[^1] and DNA:RNA[^2] duplexes. The melting temperature calculation is based on Santalucia, 1998[^4]. Sodium and cagnesium concentration correction is based on the work of Owczarzy et al[^5][^6]. Formamide correction can be performed based on two different published models[^7][^8].

* Use the `-t` option to specify the **type of nucleic acid duplex**.
* Use `--fa-mode` to switch between linear melting temperature **formamide-based correction**[^7] and linear &delta;G formamide-based correction[^8].
* Use `--fa-mvalue` together with `--fa-mode wright` to specify the **m-value** for the formamide-based correction.
* Use the `-F` option and providing the path to a **fasta file** instead of a single sequence, the melting temperature is automatically calculated for every sequence in the file.
* Use the `-v` option to trigger the verbose mode, which provides more **details** for every single sequence.
* Use `-C` for the temperature in **degree Celsius** instead of Kelvin.
* Use `--out-curve` to specify a file where to save estimated single-sequence **melting curves** with temperature range and step around the melting temperature as defined with `--t-curve`.

### Help page

```
usage: oligomelt.py [-h] [-t {DNA:DNA,RNA:RNA,RNA:DNA,DNA:RNA}]
                    [-o oligo_conc] [-n na_conc] [-m mg_conc]
                    [-f fa_conc] [--fa-mode fa_mode] [--fa-mvalue m]
                    [--t-curve range step] [--out-curve outname]
                    [-C] [-F] [-v]
                    seq

Calculate melting temeprature of a DNA duplex at provided [oligo],
[Na+], [Mg2+]. Either provide an oligo sequence or a file with one oligo
per line (and use -F option). References:
 [1] Freier et al, PNAS(83), 1986;
 [2] Sugimoto et al, Biochemistry(34), 1995.
 [3] Allawi & Santalucia, Biochemistry(36), 1997;
 [4] SantaLucia, PNAS(95), 1998;
 [5] Owczarzy et al, Biochemistry(43), 2004;
 [6] Owczarzy et al, Biochemistry(47), 2008;
 [7] McConaughy et al, Biochemistry(8), 1969;
 [8] Wright et al, Appl. env. microbiol.(80), 2014.

positional arguments:
  seq                   DNA duplex sequence (one strand only) or path to a
                        FASTA file (use with -F).

optional arguments:
  -h, --help            show this help message and exit
  -t {DNA:DNA,RNA:RNA,RNA:DNA,DNA:RNA}, --type {DNA:DNA,RNA:RNA,RNA:DNA,DNA:RNA}
                        Duplex type. Possible values: DNA:DNA (based on ref.3,
                        default), RNA:RNA (based on ref.1), DNA:RNA (based on
                        ref.2., given DNA sequence) or RNA:DNA (based on
                        ref.2, given RNA sequence). The first nucleic acid
                        type indicates the provided sequence.
  -o oligo_conc, --oconc oligo_conc
                        Oligonucleotide concentration [M]. Default: 0.25e-6 M
  -n na_conc, --naconc na_conc
                        Na+ concentration [M]. Default: 50e-3 M
  -m mg_conc, --mgconc mg_conc
                        Mg2+ concentration [M]. Note: Mg2+ correction
                        overwrites Na+ correction. Default: 0 M
  -f fa_conc, --faconc fa_conc
                        Formamide concentration in %(v,v).
  --fa-mode fa_mode     Mode of formamide correction. "mcconaughy" for
                        classical -0.72%FA correction from ref. 7, "wright"
                        for single reaction model correction from ref.8
                        (default).
  --fa-mvalue m         Specify the formamide m-value to be used with the
                        wright correction model. Use either a single value "x"
                        or two values "xL+y" where L is the probe length.
                        Default: 0.1734
  --t-curve range step  Temperature range and step for melting curve
                        generation. Use --make-curve option to generate the
                        curve. Default: 20 degC range and 0.5 degC step.
  --out-curve outname   Path to output table containing tabulated curves.
  -C, --celsius         Output temperature in Celsius degrees. Default: Kelvin
  -F, --usefile         Use when a file path is provided instead of a single
                        sequence.
  -v, --verbose         Verbose output.
```

### References

* [^1]: Freier et al, PNAS(83), 1986;
* [^2]: Sugimoto et al, Biochemistry(34), 1995.
* [^3]: Allawi & Santalucia, Biochemistry(36), 1997;
* [^4]: SantaLucia, PNAS(95), 1998;
* [^5]: Owczarzy et al, Biochemistry(43), 2004;
* [^6]: Owczarzy et al, Biochemistry(47), 2008;
* [^7]: McConaughy et al, Biochemistry(8), 1969;
* [^8]: Wright et al, Appl. env. microbiol.(80), 2014.

---

<small>This project comes from the [potpourri](https://github.com/ggirelli/potpourri) sandbox.</small>
