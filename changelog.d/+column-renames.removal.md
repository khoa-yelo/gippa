**GIPPA 1.1 renames ten table columns** (one name per quantity). This is a one-time exception to the rule that column changes wait for a major version. Old names are not kept as aliases; `gippa report` and `gippa compare` still read 1.0 output folders by renaming on read. The renames:

- `tables/novel_orfs.tsv`: `ostir_tir` → `tir_ostir`
- `tables/determinants_features.tsv`: `mw` → `mw_da`
- `tables/determinants_features.tsv`: `promotech_max` → `promoter_score`
- `tables/determinants_features.tsv`: `promotech_tu` → `promoter_score_tu`
- `mapping/peptide_hits.tsv`: `peptide_len` → `peptide_len_aa`
- `tables/pairs.tsv`: `class` → `pair_class`
- `tables/pairs.tsv`: `promoter_between` → `promoter_score_between`
- `tables/pairs.tsv`: `terminator_between` → `terminator_score_between`
- `tables/operon_flags.tsv`: `uniop_class` → `pair_class`
- `tables/operon_flags.tsv`: `flag` → `flags`

Feature names held as values change the same way (`promotech_max` → `promoter_score`, `promotech_tu` → `promoter_score_tu`, `mw` → `mw_da`) in `tables/determinants_model.tsv` › `feature_or_block` and in the `feature` column of `plots/support/determinants_pd.tsv` and `determinants_shap.tsv`; the partial-dependence axis of the determinants plot now reads "promoter score". `pairs.tsv` and `operon_flags.tsv` now share `pair_class` and `flags`, so a join on `pair_id` needs suffixes. The output schema version is now 1.3. See "Changes in 1.1" on the outputs page.
