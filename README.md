<p align="center"><img src="images/logo_transparent.png" alt="Polypolish" width="70%"></p>

Polypolish is a tool for polishing genome assemblies with short reads. Unlike other polishers, Polypolish uses SAM files where each read has been aligned to _all_ possible locations (not just a single best location). This allows it to repair errors in repeat regions that other alignment-based polishers cannot fix. Polypolish is also a conservative polisher, so it is very unlikely to introduce new errors during polishing.

For installation instructions, usage, deeper explanations and more, head over to the [Polypolish wiki](https://github.com/rrwick/Polypolish/wiki)!

This is the original paper describing Polypolish:<br>
[Wick RR, Holt KE. Polypolish: short-read polishing of long-read bacterial genome assemblies. PLOS Computational Biology. 2022. doi:10.1371/journal.pcbi.1009802.](https://doi.org/10.1371/journal.pcbi.1009802)

And this is a follow-up paper that describes Polypolish v0.6.0:<br>
[Bouras G, Judd LM, Edwards RA, Vreugde S, Stinear TP, Wick RR. How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies. Microbial Genomics. 2024. doi:10.1099/mgen.0.001254.](https://doi.org/10.1099/mgen.0.001254)

[![License GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
