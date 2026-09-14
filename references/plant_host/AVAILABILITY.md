# Plant host genome availability (NCBI assembly, checked 2026-09-08)

Host panel from collaborator's rating table (31 species). Only the 9 sampled host codes
matter for mapping: CC54, CC547, Eu635, Geu182, HC83, M26, P46, T29, T49 — code->species
mapping UNCONFIRMED, ask collaborator. Working hypothesis: code prefix = host
(Eu=Euonymus, Geu=Golden Euonymus, T=Tomato/Tobacco, P=Pepper/Poncirus, HC=Hibiscus,
CC=Carrizo Citrus?, M=Mustard?).
Update 2026-09-14: samples.tsv now holds the confirmed hosts. CC = Carrizo citrange
(C. sinensis x P. trifoliata F1) and HC = Hamlin sweet orange; see CITRUS_HOST_PLAN.md §1.1.

## Chromosome-level assembly available (good to use directly)
- Helianthus annuus (sunflower), Brassica oleracea, Brassica juncea, Carica papaya,
  Quercus virginiana, Hydrangea macrophylla, Persea americana (avocado),
  Ficus benjamina, Citrus sinensis, Bergera koenigii (curry tree),
  Poncirus trifoliata, Salix babylonica, Capsicum annuum, Nicotiana benthamiana,
  Solanum lycopersicum, Vitis rotundifolia (muscadine)

## Species missing -> nearest usable relative
- Euonymus japonicus (both euonymus entries!) -> Euonymus europaeus (chromosome-level)
- Hibiscus rosa-sinensis / acetosella -> Hibiscus yunnanensis (7 chrom-level in genus)
- Citrus paradisi (grapefruit) / macrophylla / medica -> other Citrus (many)
- Carrizo citrange (C. sinensis x P. trifoliata) -> concatenate both parents' genomes
- Vitis labrusca (scaffold-only) -> V. rotundifolia or V. vinifera
- Rubus fruticosus -> Rubus spectabilis etc. (9 in genus)

## No genome anywhere in genus (would need de novo or distant relative)
- Schinus terebinthifolia, Verbesina virginica, Kalanchoe delagoensis,
  Acacia auriculiformis (genus-level not checked beyond), Euphorbia pulcherrima
  (1 scaffold-level only), Vaccinium virgatum (congeners exist, e.g. V. corymbosum)

## Next step
Once codes are confirmed: download the 9 needed assemblies + annotation via NCBI
datasets API (works on compute nodes through the proxy), place one subdir per host here,
prefix contigs `plant_<code>_` via scripts/convert_refs.py conventions / seqkit.
