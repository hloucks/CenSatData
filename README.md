## Centromeric Satellite Data Repository
CenSat Track Generated with [CenSat Annotation Workflow](https://github.com/kmiga/alphaAnnotation/blob/CenSat/cenSatAnnotation/centromereAnnotation.wdl)  

### Annotation Bin Overview 

Alpha-Satellites - Annotated with [Fedor Ryabov’s HumAS-HMMER](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL) 
- Active alpha (active_hor) 
- diverged HORS (dhor) 
- monomeric HORs (mon)
- mixed alpha (mixedAlpha)

Human Satellites 2 and 3 - Annotated with [Nick Altemose’s HSAT2/3 script](https://github.com/altemose/chm13_hsat)
	
Other Centromeric Satellite annotations - Annotated with RepeatMasker
- HSAT1A - SAR in DFAM
- HSAT1B - HSAT1 in DFAM
- Gamma - includes all GSAT and TAR1 in DFAM
- Beta - BSR, LSAU, and BSAT in DFAM
- CenSat - other centromeric satellites CER, SATR, SST1, ACRO, HSAT4, HSAT5, TAF11  

Centromere Transition (ct)  
Centromeres are defined by merging all above satellite annotations within 2MB (bedtools merge) and then identifying the region containing the active array. Any stretch of sequence not annotated within this region is marked "ct"


Please contact hloucks@ucsc.edu with any questions 