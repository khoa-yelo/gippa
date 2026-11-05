# ruff: noqa: RUF001  (the wording is Unicode: en dashes, ≤, ≥)
"""The one plain-text dictionary for user-facing prose (M12-04 PR 0, § 1 rule 2).

Flag tokens, skip reasons, feature columns and detection tiers map to plain words here, so
figure titles, captions, cards and notices never show a raw token. Pure: constants and one
lookup function, no file I/O.

A key ending in ``:`` is a prefix template: ``layer_missing:<layer>`` with no exact entry is
looked up by ``layer_missing:`` and its parameter fills the ``{}`` of the sentence. Values hold
no underscore, no section sign and no internal reference (AC 0.3).
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from types import MappingProxyType

__all__ = [
    "FEATURE_LABEL",
    "FLAG_TEXT",
    "SCOPED_TEXT",
    "SKIP_TEXT",
    "TIER_NAME",
    "TIER_TEXT",
    "feature_label",
    "plain",
]

_log = logging.getLogger(__name__)

#: flag token (or ``prefix:`` template) -> short sentence; table flags and run flags
FLAG_TEXT: Mapping[str, str] = MappingProxyType(
    {
        # genes
        "partial": "gene or ORF runs off a contig end",
        "crosses_origin": "gene spans the origin of a circular replicon",
        "pseudogene": "annotated as a pseudogene",
        "multi_identical": "peptide evidence shared by identical proteins",
        "shared_evidence": "evidence shared with another gene",
        "two_proteins_one_gene": "two input proteins map to this one gene",
        "cross_strain_protein": "protein matched from another strain",
        "imputed_abundance": "abundance imputed, not measured",
        "abundance_imputed": "abundance imputed, not measured",
        "tolerant_evidence": "has evidence from mismatch-tolerant matches",
        "strain_difference": "evidence points to a sequence difference from the reference strain",
        "nterm_extension_evidence": "peptides suggest an N-terminal extension",
        "cai_reference": "in the codon adaptation reference set (ribosomal proteins); not used to "
        "fit the model",
        "top3_lt3": "Top3 abundance from fewer than three peptides",
        "detectability_constant": "detectability model not fitted; a constant rate is used",
        "abundance_not_explained_by_sequence": "abundance not explained by sequence features",
        # placement and peptides
        "same_db": "peptide searched against this genome's own database",
        "fs": "placement crosses a frameshift",
        "st": "placement reads through a stop codon",
        "engine_protein_mismatch": "search engine protein differs from the mapped gene",
        "origin_spanning": "peptide spans the origin of a circular replicon",
        "initiator_m": "peptide starts with the initiator methionine at an annotated non-ATG "
        "start codon",
        "low_complexity": "low-complexity peptide",
        "multi": "peptide maps to more than one place",
        "pyl_site": "peptide reads a TAG codon as pyrrolysine",
        "sec_site": "peptide reads a TGA codon as selenocysteine",
        "readthrough": "peptide reads through a stop codon",
        # regions and context
        "small_region": "region too small for a reliable estimate",
        "dosage_plus_regulation": "gene dosage and regulation are not separated",
        "descriptive": "descriptive only; no test",
        "no_oric": "replication origin not found",
        "uncalibrated": "test not calibrated for this setting",
        "uncalibrated:": "{} is not calibrated in this run",
        "uncalibrated:operon_distance": "the operon distance model is not calibrated in this run",
        "uncalibrated:uniop_pair_threshold": "the UniOp pair threshold is not calibrated in "
        "this run",
        "uncalibrated:terminatornet_yes": "the TerminatorNet call threshold is not calibrated in "
        "this run",
        "depth_adjusted": "adjusted for plasmid copy number",
        "few_shifts": "too few shifts for a significant result",
        "few_values": "too few values to scan",
        "split_contigs": "tested separately on each contig",
        "hl_variance_nonpositive": "variance estimate not positive; no p-value",
        "no_reference_replicon": "no reference replicon to compare with",
        "not_estimable": "statistic cannot be estimated",
        "ci_few_units": "fewer than 10 resampling units; the interval is unreliable or not "
        "computed",
        "ci_unestimated": "confidence interval not estimated",
        "tu_block_tu_error": "the transcription unit null assumes exact unit boundaries",
        "tu_source_distance_model": "transcription units from the distance model",
        # operons
        "rna_between": "an RNA gene lies between the pair",
        "relative_scale": "abundance on a relative scale; not tested",
        "internal_promoter_supported": "a promoter inside the operon explains higher downstream "
        "abundance",
        "internal_terminator_supported": "a terminator inside the operon explains lower downstream "
        "abundance",
        "unexplained_discordance": "discordant abundance inside an operon, not explained by a "
        "promoter or terminator",
        "conflict_terminator_in_operon": "a terminator lies inside an operon with concordant "
        "abundance",
        "novel_operon_candidate": "co-abundant neighbours not predicted as an operon",
        "coabundant_despite_terminator": "co-abundant despite a terminator between them",
        "other": "another evidence combination",
        "no_detection_limit_model": "no detection limit model; one-detected pairs are not tested",
        "pmiss_covariates_missing": "covariates for the detection model are missing",
        "too_few_no_pairs": "too few non-operon pairs for a null",
        "too_few_op_pairs": "too few operon pairs to fit the model",
        "unestimated": "not estimated",
        # enrichment
        "elim_empty": "no genes left after removing more specific terms",
        "null_degenerate": "the detection null is degenerate; no p-value",
        "null_uncalibrated": "the permutation p-value is not calibrated "
        "(it assumes operons are exchangeable)",
        "too_few_rows": "too few genes to fit the model",
        "odds_out_of_range": "odds ratio out of range; no p-value",
        "ora_uncalibrated": "over-representation test not calibrated",
        "p_detect_clipped": "detection probabilities of 0 or 1 were clipped",
        "tu_recal_skipped": "transcription unit recalibration skipped",
        "tu_recal_reduced": "transcription unit recalibration reduced",
        "detection_confounded": "the result may reflect detectability, not abundance",
        "genomic_cluster": "set members cluster along the genome",
        "log2err_unestimated": "precision of the p-value not estimated",
        "size_filtered": "set too small or too large to test",
        "too_few_ranked_tus": "too few ranked transcription units",
        "no_eligible_genes": "no eligible genes",
        # comparison
        "both_empty": "both sets are empty",
        "missing_in:": "absent from run {}",
        "scale_differs": "the runs use different abundance measures or scales",
        "too_few_shared": "fewer than three shared genes",
        # variants
        "mass_equivalent": "variant has the same mass as the reference residue",
        "snv": "single nucleotide change",
        "multi_nt": "change of more than one nucleotide",
        "database_derived": "supported by a peptide from another database",
        "start_codon": "variant at a start codon",
        "stop_codon": "variant at a stop codon",
        "vcf_unrepresentable": "spans the origin; not written to the VCF",
        # novel ORFs
        "junction_spanned": "a peptide spans the pseudogene junction",
        "junction_unknown": "pseudogene junction not located",
        "homopolymer_junction": "junction at a homopolymer run",
        "class_fdr_unestimated": "class FDR not estimated",
        "rbs_supported": "supported by a ribosome binding site",
        "promoter_supported": "supported by a promoter",
        "terminator_supported": "supported by a terminator",
        "denovo_only_unsupported": "de novo peptides only, without sequence support",
        "no_start_codon": "no start codon found",
        "ostir_missing": "no translation initiation rate computed",
        "rbs_threshold_unestimated": "ribosome binding site threshold not estimated",
        "readthrough_candidate": "possible stop codon readthrough",
        # run
        "anti_sd_default": "the E. coli anti-Shine–Dalgarno sequence was used",
        "c_assumed": "de novo precision assumed, not measured",
        "unreported": "de novo precision not estimable in this score bin; not reported",
        "bakta_db_light": "annotated with the light Bakta database",
        "operon_source_distance_model": "operons from the distance model",
        "input_prefiltered": "input was filtered before GIPPA",
        "cai_reference_missing": "no reference gene set for codon adaptation",
        "foreign_protein_ids": "proteins come from a related annotation; identical matches are not "
        "flagged per gene",
        "layer_missing:": "the {} layer was not computed",
        "layer_missing:abundance": "abundance not available",
        "layer_missing:antismash": "biosynthetic gene clusters not predicted",
        "layer_missing:codon": "codon usage not computed",
        "layer_missing:cog": "COG categories not assigned",
        "layer_missing:features": "gene features not computed",
        "layer_missing:genomad": "mobile elements not predicted",
        "layer_missing:go": "GO terms not assigned",
        "layer_missing:hits_tolerant": "mismatch-tolerant matching not run",
        "layer_missing:kegg": "KEGG pathways not assigned",
        "layer_missing:ostir": "translation initiation rates not computed",
        "layer_missing:promoters": "promoters not predicted",
        "layer_missing:rbs": "ribosome binding sites not predicted",
        "layer_missing:regions": "genomic regions not computed",
        "layer_missing:terminators": "terminators not predicted",
        "layer_missing:tx_rates": "transcription rates not computed",
    }
)

#: table -> flag token -> sentence, where a token means something else in that table
SCOPED_TEXT: Mapping[str, Mapping[str, str]] = MappingProxyType(
    {
        "pairs": MappingProxyType(
            {"abundance_imputed": "abundance of a gene of the pair imputed, not measured"}
        ),
        "determinants_model": MappingProxyType(
            {"null_degenerate": "the permutation null is degenerate; no p-value"}
        ),
    }
)

#: skip-reason token (or ``prefix:`` template) -> sentence
SKIP_TEXT: Mapping[str, str] = MappingProxyType(
    {
        "module_off:": "the {} module was not run",
        "input_absent:": "an input this view needs is missing from this run",
        "input_absent:autocorrelation": "no autocorrelation table in this run",
        "input_absent:detection_summary": "no detection summary in this run",
        "input_absent:determinants_model": "no determinants model in this run",
        "input_absent:genes_tsv": "no gene table in this run",
        "input_absent:genes": "no gene table in this run",
        "input_absent:mapping_fdr": "no mapping FDR table in this run",
        "input_absent:novel_orfs": "no novel ORF table in this run",
        "input_absent:pairs": "no gene pair table in this run",
        "input_absent:peptide_hits": "no peptide table in this run",
        "input_absent:preflight": "no input check in this run",
        "input_absent:unmapped": "no unmapped peptide table in this run",
        "input_absent:variants": "no variant table in this run",
        "not_registered": "this plot is not available in this GIPPA version",
    }
)

#: feature column -> human name
FEATURE_LABEL: Mapping[str, str] = MappingProxyType(
    {
        # identity columns of determinants_features.tsv
        "locus_tag": "locus tag",
        "tu_id": "transcription unit",
        "detected": "detection status",
        "abundance_log2": "log2 abundance",
        "abundance_imputed": "abundance imputed",
        "ribosomal_protein": "ribosomal protein",
        # features
        "promotech_max": "promoter score",
        "promotech_tu": "promoter score of the transcription unit",
        "promoter_score": "promoter score",
        "promoter_score_tu": "promoter score of the transcription unit",
        "tx_rate_max_tu": "highest transcription rate in the transcription unit",
        "tx_rate_sum_tu": "summed transcription rate of the transcription unit",
        "n_promoters_tu_to_gene": "promoters between unit start and gene",
        "tu_position": "position in the transcription unit",
        "dist_to_tu_start_bp": "distance to the transcription unit start (bp)",
        "terminator_conf_downstream": "downstream terminator confidence",
        "tir_ostir": "translation initiation rate (OSTIR)",
        "sd_dg": "Shine–Dalgarno binding energy",
        "spacer_bp": "spacer length (bp)",
        "start_codon": "start codon",
        "mrna_dg": "mRNA folding energy",
        "cai": "codon adaptation index",
        "tai": "tRNA adaptation index",
        "enc": "effective number of codons",
        "fop": "frequency of optimal codons",
        "gc3": "GC content at third codon positions",
        "cpb": "codon pair bias",
        "at_2_10": "AT content, codons 2–10",
        "rare_codons_2_10": "rare codons, codons 2–10",
        "length_aa": "protein length (aa)",
        "mw": "molecular weight",
        "mw_da": "molecular weight (Da)",
        "pi": "isoelectric point",
        "gravy": "hydropathy (GRAVY)",
        "nterm_residue": "N-terminal residue",
        "tm_segments": "transmembrane segments",
        "n_observable_peptides": "observable peptides",
        "replicon_type": "replicon type",
        "replichore_position": "position on the replichore",
        "in_prophage": "in a prophage",
        "in_bgc": "in a biosynthetic gene cluster",
        "gc_local_dev": "local GC deviation",
    }
)

#: detection tier -> sentence without the tier name (the caller adds it); ``{alpha}`` is the
#: run's FDR level
TIER_TEXT: Mapping[str, str] = MappingProxyType(
    {
        "decoy_picked": "a gene counts as detected when its gene-level q-value "
        "(target–decoy competition) is ≤ {alpha}.",
        "mayu": "a gene counts as detected when it has an accepted peptide at q ≤ {alpha}; "
        "the gene-level FDR of that list is estimated with Mayu.",
        "two_peptide": "a gene counts as detected when ≥ 2 distinct accepted peptides map to it; "
        "no FDR is estimated.",
        "protein_list": "a gene counts as detected when a protein in your list maps to it; "
        "no FDR is estimated.",
    }
)

#: detection tier -> its name in prose
TIER_NAME: Mapping[str, str] = MappingProxyType(
    {
        "decoy_picked": "decoy-picked",
        "mayu": "Mayu",
        "two_peptide": "two-peptide",
        "protein_list": "protein list",
    }
)

#: a template key looked up bare (no parameter) -> sentence without the placeholder
_BARE: Mapping[str, str] = MappingProxyType(
    {
        "layer_missing:": "a layer was not computed",
        "uncalibrated:": "a setting not calibrated in this run",
        "missing_in:": "absent from another run",
        "module_off:": "a module was not run",
    }
)
#: template prefixes whose parameter is user text (a run name), filled in verbatim
_VERBATIM = frozenset({"missing_in"})


def _input_key(param: str) -> str:
    """``tables/variants.tsv`` -> ``variants``: an output path reduced to its table key."""
    name = param.rsplit("/", 1)[-1]
    for ext in (".gz", ".tsv", ".json", ".vcf"):
        name = name.removesuffix(ext)
    return name


def _fill(template: str, prefix: str, param: str) -> str:
    """``template`` with its ``{}`` filled by ``param``: underscores become spaces (a parameter is
    a token, prose holds none) except for a verbatim prefix; an empty parameter takes the bare
    sentence."""
    if prefix not in _VERBATIM:
        param = " ".join(param.replace("_", " ").split())
    if not param:
        return _BARE.get(f"{prefix}:", template.replace("{}", "").strip())
    return template.replace("{}", param)


def _lookup(token: str) -> str | None:
    if token in _BARE:
        return _BARE[token]
    for table in (FLAG_TEXT, SKIP_TEXT):
        if token in table:
            return table[token]
    if ":" not in token:
        return None
    prefix, param = token.split(":", 1)
    if not param:
        return None
    if prefix == "input_absent":
        exact = SKIP_TEXT.get(f"input_absent:{_input_key(param)}")
        if exact is not None:
            return exact
    for table in (FLAG_TEXT, SKIP_TEXT):
        template = table.get(f"{prefix}:")
        if template is not None:
            return _fill(template, prefix, param)
    return None


def plain(token: str, table: str | None = None) -> str:
    """The plain-words sentence for a flag token or skip reason; ``table`` (a table name such as
    ``pairs``) selects the sentence of a token that means something else there and otherwise
    changes nothing. An unknown token comes back unchanged with a logged warning. Never raises."""
    token = str(token)
    if table is not None:
        scoped = SCOPED_TEXT.get(str(table), {}).get(token)
        if scoped is not None:
            return scoped
    text = _lookup(token)
    if text is None:
        _log.warning("no plain-text entry for %r", token)
        return token
    return text


def feature_label(col: str) -> str:
    """The human name of a determinants column: a ``FEATURE_LABEL`` entry, or a one-hot column
    ``<feature>_<level>`` as ``"<feature name>: <level>"`` (``start_codon_ATG`` ->
    ``"start codon: ATG"``); any other column with its underscores turned into spaces."""
    col = str(col)
    label = FEATURE_LABEL.get(col)
    if label is not None:
        return label
    for feat in sorted(FEATURE_LABEL, key=len, reverse=True):
        if col.startswith(f"{feat}_") and len(col) > len(feat) + 1:
            level = col[len(feat) + 1 :].replace("_", " ")
            return f"{FEATURE_LABEL[feat]}: {level}"
    return col.replace("_", " ")
