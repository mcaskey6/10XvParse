from __future__ import annotations

import logging
from pathlib import Path

import yaml

from .classes import RunSettings, AnalysisConfig, IndexConfig, TenXPaths, HashtagsPaths, ParsePaths
from . import utils, parse_config, io


def _load_index_config(settings: RunSettings, config: AnalysisConfig) -> IndexConfig:
    index_config_file = settings.root_dir / "Configs" / "indexes" / f"{config.species}.yaml"
    return IndexConfig.from_yaml(index_config_file)


def _download_reference(settings: RunSettings, paths: TenXPaths, index_config: IndexConfig, logger: logging.Logger) -> None:
    '''Download reference genome. Builds a barnyard (dual-species) reference when index_config.is_barnyard.'''
    ref_exists = paths.genome_file.is_file() and paths.gtf_file.is_file()
    if not ref_exists or settings.overwrite:
        if index_config.is_barnyard:
            utils.make_barnyard_reference(paths, index_config, logger)
        else:
            utils.get_reference(paths.genome_file, paths.gtf_file, index_config.fasta_url, index_config.gtf_url, logger)
    else:
        logger.info("Reference files already exist. Skipping reference download.")


def _run_core_pipeline(settings: RunSettings, paths, config: AnalysisConfig, assay: str, logger: logging.Logger) -> None:
    '''Route to ERA, SRA, or local-files pipeline based on what accessions the config provides.'''
    if config.era:
        utils.era_core_pipeline(settings, paths, config, assay, logger)
    elif config.sra:
        utils.core_pipeline(settings, paths, config, assay, logger)
    else:
        utils.local_pipeline(settings, paths, config, assay, logger)


def load_10x(settings: RunSettings, config_file: str, assay: str, logger: logging.Logger) -> None:
    '''Get H5AD files for 10X data. Supports SRA and ERA sources, single-species and barnyard references.'''
    config = AnalysisConfig.from_yaml(config_file, assay)
    paths = TenXPaths.build(settings, config, assay)
    paths.ensure_dirs(logger)

    index_config = _load_index_config(settings, config)
    _download_reference(settings, paths, index_config, logger)
    _run_core_pipeline(settings, paths, config, assay, logger)

    if settings.run_kb:
        utils.pseudoalign_10x(
            paths=paths,
            fastq_files=paths.multiplexed_files,
            kb_out_dir=paths.kb_all_dir,
            tech=config.technology,
            threads=settings.threads,
            logger=logger,
        )


def load_10x_hashtags(settings: RunSettings, config_file: str, hashtags: str, logger: logging.Logger) -> None:
    '''Get H5AD files for 10X Hashtags for Demultiplexing'''
    config = AnalysisConfig.from_yaml(config_file, hashtags)
    paths = HashtagsPaths.build(settings, config, hashtags)
    paths.ensure_dirs(logger)

    _run_core_pipeline(settings, paths, config, hashtags, logger)

    if settings.run_kb:
        utils.pseudoalign_10x_hashtags(
            paths=paths,
            fastq_files=paths.multiplexed_files,
            kb_out_dir=paths.kb_all_dir,
            tech=config.technology,
            threads=settings.threads,
            trim_length=config.r2_trim_length,
            logger=logger,
        )


def load_parse(settings: RunSettings, config_file: str, assay: str, logger: logging.Logger) -> None:
    '''Get H5AD files for Parse data. Supports SRA and ERA sources, single-species and barnyard references.'''
    config = AnalysisConfig.from_yaml(config_file, assay)
    paths = ParsePaths.build(settings, config, assay)
    paths.ensure_dirs(logger)

    if not parse_config.is_parse_kit(config.technology):
        raise ValueError(
            f"'{config.technology}' is not a valid Parse kit name for assay '{assay}'. "
            "Expected format: '<kit>_v<chem>', e.g. 'WT_v2' or 'WT_mini_v3'."
        )
    x_string = parse_config.generate_parse_configs(
        kit_name=config.technology,
        parse_info_dir=paths.parse_info_dir,
        output_dir=paths.config_dir / config.name / assay,
        wells=config.wells or None,
        logger=logger,
    )
    bc1_location = parse_config.x_string_to_bc1_location(x_string)

    index_config = _load_index_config(settings, config)
    _download_reference(settings, paths, index_config, logger)
    _run_core_pipeline(settings, paths, config, assay, logger)

    filter_parse_fastqs_exist = all(p.is_file() for p in paths.filtered_files)
    if not filter_parse_fastqs_exist or settings.overwrite:
        utils.filter_parse_fastqs(
            paths=paths,
            threads=settings.threads,
            logger=logger,
            bc1_location=bc1_location,
        )

    # Split parse FASTQ into polyT and randO files
    randOpolyT_exist = all(p.is_file() for p in paths.randO_files + paths.polyT_files)
    if not randOpolyT_exist or settings.overwrite:
        utils.extract_rando_polyt(
            paths=paths,
            threads=settings.threads,
            logger=logger
        )

    # Pseudoalign reads separated by barcode type
    if settings.run_kb:
        utils.pseudoalign_parse(
            paths=paths,
            fastq_files=paths.filtered_files,
            kb_out_dir=paths.kb_all_dir,
            tech=x_string,
            threads=settings.threads,
            logger=logger)


def subsample_parse(settings: RunSettings, config_file: str, assay: str, subsample_num: int, logger: logging.Logger) -> None:
    '''Subsample the parse, polyT and randO FASTQ files for comparison.'''
    config = AnalysisConfig.from_yaml(config_file, assay)
    paths = ParsePaths.build(settings, config, assay)
    paths.ensure_dirs(logger)

    if not parse_config.is_parse_kit(config.technology):
        raise ValueError(
            f"'{config.technology}' is not a valid Parse kit name for assay '{assay}'. "
            "Expected format: '<kit>_v<chem>', e.g. 'WT_v2' or 'WT_mini_v3'."
        )
    x_string = parse_config.generate_parse_configs(
        kit_name=config.technology,
        parse_info_dir=paths.parse_info_dir,
        output_dir=paths.config_dir / config.name / assay,
        wells=config.wells or None,
        logger=logger,
    )

    sampled_parse_exist = all(p.is_file() for p in paths.sampled_files)
    if not sampled_parse_exist or settings.overwrite:
        utils.subsample_fastqs(
            fastq_files=paths.filtered_files,
            output_files=paths.sampled_files,
            num_reads=subsample_num,
            threads=settings.threads,
            logger=logger
        )

    sampled_polyT_exist = all(p.is_file() for p in paths.sampled_polyT_files)
    if not sampled_polyT_exist or settings.overwrite:
        utils.subsample_fastqs(
            fastq_files=paths.polyT_files,
            output_files=paths.sampled_polyT_files,
            num_reads=subsample_num,
            threads=settings.threads,
            logger=logger
        )

    sampled_randO_exist = all(p.is_file() for p in paths.sampled_randO_files)
    if not sampled_randO_exist or settings.overwrite:
        utils.subsample_fastqs(
            fastq_files=paths.randO_files,
            output_files=paths.sampled_randO_files,
            num_reads=subsample_num,
            threads=settings.threads,
            logger=logger
        )

    if settings.run_kb:
        utils.pseudoalign_parse(
            paths=paths,
            fastq_files=paths.sampled_files,
            kb_out_dir=paths.kb_sub_dir,
            tech=x_string,
            threads=settings.threads,
            logger=logger,
        )

        utils.pseudoalign_parse(
            paths=paths,
            fastq_files=paths.sampled_polyT_files,
            kb_out_dir=paths.kb_sub_polyT_dir,
            tech=x_string,
            threads=settings.threads,
            logger=logger
        )

        utils.pseudoalign_parse(
            paths=paths,
            fastq_files=paths.sampled_randO_files,
            kb_out_dir=paths.kb_sub_randO_dir,
            tech=x_string,
            threads=settings.threads,
            logger=logger
        )


def subsample_10x(settings: RunSettings, config_file: str, assay: str, subsample_num: int, logger: logging.Logger, tag: str = "") -> None:
    '''Subsample the 10X FASTQ files for comparison.

    Use `tag` to produce multiple subsampled outputs at different depths from the same raw data
    (e.g. tag="mini" and tag="standard" for separate Parse kit comparisons). Each tag gets its
    own FASTA/Sampled_<tag>/ directory and kb_python/sampled_<tag>_out/ directory.
    '''
    config = AnalysisConfig.from_yaml(config_file, assay)
    paths = TenXPaths.build(settings, config, assay)
    paths.ensure_dirs(logger)

    if tag:
        sampled_dir = paths.fasta_dir / f"Sampled_{tag}"
        sampled_files = [sampled_dir / f"{assay}_{i}.fastq.gz" for i in range(2)]
        kb_sub_dir = paths.kb_dir / f"sampled_10x_{tag}_out"
        from .classes import make_dir
        make_dir(sampled_dir, logger)
        make_dir(kb_sub_dir, logger)
    else:
        sampled_files = paths.sampled_files
        kb_sub_dir = paths.kb_sub_dir

    sampled_exist = all(p.is_file() for p in sampled_files)
    if not sampled_exist or settings.overwrite:
        utils.subsample_fastqs(
            fastq_files=paths.multiplexed_files,
            output_files=sampled_files,
            num_reads=subsample_num,
            threads=settings.threads,
            logger=logger
        )

    if settings.run_kb:
        utils.pseudoalign_10x(
            paths=paths,
            fastq_files=sampled_files,
            kb_out_dir=kb_sub_dir,
            tech=config.technology,
            threads=settings.threads,
            logger=logger,
        )


def get_subsample_num(
    settings: RunSettings,
    config_file: str,
    ten_x_assays: list[str],
    parse_assays: list[str],
    logger: logging.Logger,
) -> int:
    '''Count reads in processed FASTQ files and return the minimum across all assays.
    Per-file counts are cached to Configs/<analysis_name>/read_counts.yaml. On subsequent
    calls, only files missing from the cache are recounted; the minimum is always computed
    from the current assay set only. Pass overwrite=True to force a full recount.'''
    config_file = Path(config_file)
    all_assays = ten_x_assays + parse_assays
    analysis_name = AnalysisConfig.from_yaml(config_file, all_assays[0]).name
    cache_file = config_file.parent / analysis_name / "read_counts.yaml"

    cached: dict[str, int] = {}
    if cache_file.is_file():
        with open(cache_file) as f:
            cached = yaml.safe_load(f) or {}

    # Collect the files required for the current assay set
    needed: dict[str, Path] = {}
    for assay in ten_x_assays:
        config = AnalysisConfig.from_yaml(config_file, assay)
        paths = TenXPaths.build(settings, config, assay)
        f = paths.multiplexed_files[0]
        needed[str(f)] = f
    for assay in parse_assays:
        config = AnalysisConfig.from_yaml(config_file, assay)
        paths = ParsePaths.build(settings, config, assay)
        for f in [paths.filtered_files[0], paths.polyT_files[0], paths.randO_files[0]]:
            needed[str(f)] = f

    # Use cached counts where available; recount missing or stale files
    counts: dict[str, int] = {}
    cache_updated = False
    for path_str, path in needed.items():
        if path_str in cached and not settings.overwrite:
            counts[path_str] = cached[path_str]
            logger.info("%s: %d reads (cached)", path_str, counts[path_str])
        else:
            n = utils.count_reads(path)
            counts[path_str] = n
            cached[path_str] = n
            logger.info("%s: %d reads", path_str, n)
            cache_updated = True

    if cache_updated:
        with open(cache_file, "w") as f:
            yaml.dump(cached, f)
        logger.info("Updated read count cache at %s", cache_file)

    minimum = min(counts.values())
    logger.info("Minimum reads across all assays: %d", minimum)
    return minimum
    
def get_genebody_plot(
    settings: RunSettings, 
    config_file: str, 
    tenx_assay: str, parse_assay: str, 
    logger: logging.Logger,
    tag: str = None
) -> None:
    
    '''Generate genebody coverage plot for 10X and Parse data.'''
    star_prefixes = []

    config_10x = AnalysisConfig.from_yaml(config_file, tenx_assay)
    paths_10x = TenXPaths.build(settings, config_10x, tenx_assay)
    paths_10x.ensure_dirs(logger)

    star_prefixes.append(
        utils.run_star_10x(
            paths=paths_10x,
            settings=settings,
            assay=tenx_assay,
            logger=logger,
            tag=tag,
            overwrite=settings.run_kb,
        )
    )

    config_parse = AnalysisConfig.from_yaml(config_file, parse_assay)
    paths_parse = ParsePaths.build(settings, config_parse, parse_assay)
    paths_parse.ensure_dirs(logger)

    star_prefixes.append(
        utils.run_star_parse(
            paths=paths_parse,
            config=config_parse,
            settings=settings,
            assay=parse_assay,
            fastq_files=paths_parse.filtered_files,
            logger=logger,
            tag="parse",
            overwrite=settings.run_kb
        )
    )

    star_prefixes.append(
        utils.run_star_parse(
            paths=paths_parse,
            config=config_parse,
            settings=settings,
            assay=parse_assay,
            fastq_files=paths_parse.polyT_files,
            logger=logger,
            tag="polyT",
            overwrite=settings.run_kb
        )
    )

    star_prefixes.append(
        utils.run_star_parse(
            paths=paths_parse,
            config=config_parse,
            settings=settings,
            assay=parse_assay,
            fastq_files=paths_parse.randO_files,
            logger=logger,
            tag="randO",
            overwrite=settings.run_kb
        )
    )

    index_config = _load_index_config(settings, config_10x)
    bam_files = [prefix + "Aligned.sortedByCoord.out.bam" for prefix in star_prefixes]
    utils.generate_genebody_plot(
        settings=settings,
        gtf_file=paths_10x.gtf_file,
        bed_file=paths_10x.bed_file,
        bam_files=bam_files,
        gene_body_dir=paths_10x.plot_dir,
        logger=logger,
        hk_genes_file=paths_10x.hk_genes_file,
        hrt_atlas_url=index_config.hrt_atlas_url,
        hrt_atlas_url_2=index_config.hrt_atlas_url_2 or None,
        hk_bed_file=paths_10x.hk_bed_file,
        tag=tag,
    )
