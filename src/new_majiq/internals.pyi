from typing import Any, ClassVar, Dict, Iterable, Iterator, List, Optional, overload

import numpy
import numpy.typing as npt

class Contigs:
    def __init__(self, seqids: list) -> None: ...
    def checksum(self) -> int: ...
    def __contains__(self, arg0: str) -> bool: ...
    def __getitem__(self, seqid: str) -> int: ...
    def __len__(self) -> int: ...
    @property
    def seqid(self) -> List[str]: ...

class Events:
    def __init__(
        self,
        introns: GeneIntrons,
        junctions: GeneJunctions,
        ref_exon_idx: npt.ArrayLike,
        event_type: npt.ArrayLike,
        offsets: npt.ArrayLike,
        is_intron: npt.ArrayLike,
        connection_idx: npt.ArrayLike,
    ) -> None: ...
    def connection_denovo(
        self, connection_idx: npt.ArrayLike
    ) -> npt.NDArray[numpy.bool_]: ...
    def connection_end(
        self, connection_idx: npt.ArrayLike
    ) -> npt.NDArray[numpy.int64]: ...
    def connection_gene_idx(
        self, connection_idx: npt.ArrayLike
    ) -> npt.NDArray[numpy.uint64]: ...
    def connection_other_exon_idx(
        self, connection_idx: npt.ArrayLike
    ) -> npt.NDArray[numpy.uint64]: ...
    def connection_start(
        self, connection_idx: npt.ArrayLike
    ) -> npt.NDArray[numpy.int64]: ...
    def __len__(self) -> int: ...
    @property
    def _offsets(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def connection_event_idx(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def connection_idx_end(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def connection_idx_start(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def event_type(self) -> npt.NDArray: ...
    @property
    def idx(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def introns(self) -> GeneIntrons: ...
    @property
    def is_intron(self) -> npt.NDArray[numpy.bool_]: ...
    @property
    def junctions(self) -> GeneJunctions: ...
    @property
    def num_connections(self) -> int: ...
    @property
    def num_events(self) -> int: ...
    @property
    def num_introns(self) -> int: ...
    @property
    def num_junctions(self) -> int: ...
    @property
    def ref_exon_idx(self) -> npt.NDArray[numpy.uint64]: ...

class EventsAlign:
    def __init__(self, left_events: Events, right_events: Events) -> None: ...
    def events_match(self, *args, **kwargs) -> Any: ...
    @property
    def left_event_idx(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def right_event_idx(self) -> npt.NDArray[numpy.uint64]: ...

class EventsCoverage:
    def __init__(
        self,
        events: Events,
        numreads: npt.ArrayLike,
        numbins: npt.ArrayLike,
        bootstraps: npt.ArrayLike,
    ) -> None: ...
    @staticmethod
    def from_sj(
        events: Events,
        sj_junctions: SJJunctionsBins,
        sj_introns: SJIntronsBins,
        num_bootstraps: int,
        pvalue_threshold: float,
    ) -> EventsCoverage: ...
    def __len__(self) -> int: ...
    @property
    def _events(self) -> Events: ...
    @property
    def bootstraps(self) -> npt.NDArray[numpy.float32]: ...
    @property
    def numbins(self) -> npt.NDArray[numpy.float32]: ...
    @property
    def numreads(self) -> npt.NDArray[numpy.float32]: ...

class ExonConnections:
    def __init__(
        self, exons: Exons, introns: GeneIntrons, junctions: GeneJunctions
    ) -> None: ...
    def constitutive(self) -> Events: ...
    def event_description(
        self, exon_idx: npt.ArrayLike, event_type: npt.ArrayLike
    ) -> List[str]: ...
    def event_id(
        self, exon_idx: npt.ArrayLike, event_type: npt.ArrayLike
    ) -> List[str]: ...
    def event_size(
        self, exon_idx: npt.ArrayLike, is_source: npt.ArrayLike
    ) -> npt.NDArray[numpy.uint64]: ...
    def events_for(
        self, exon_idx: npt.ArrayLike, is_source: npt.ArrayLike
    ) -> Events: ...
    def has_intron(
        self, exon_idx: npt.ArrayLike, is_source: npt.ArrayLike
    ) -> npt.NDArray[numpy.bool_]: ...
    def is_LSV(
        self, exon_idx: npt.ArrayLike, is_source: npt.ArrayLike
    ) -> npt.NDArray[numpy.bool_]: ...
    def is_constitutive(
        self, exon_idx: npt.ArrayLike, is_source: npt.ArrayLike
    ) -> npt.NDArray[numpy.bool_]: ...
    def lsvs(self) -> Events: ...
    def passed(
        self, exon_idx: npt.ArrayLike, is_source: npt.ArrayLike
    ) -> npt.NDArray[numpy.bool_]: ...
    def redundant(
        self, exon_idx: npt.ArrayLike, is_source: npt.ArrayLike
    ) -> npt.NDArray[numpy.bool_]: ...
    @property
    def _exons(self) -> Exons: ...
    @property
    def _introns(self) -> GeneIntrons: ...
    @property
    def _junctions(self) -> GeneJunctions: ...
    @property
    def dst_intron_exon_offsets(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def dst_intron_idx(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def dst_junction_exon_offsets(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def dst_junction_idx(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def src_intron_exon_offsets(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def src_intron_idx(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def src_junction_exon_offsets(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def src_junction_idx(self) -> npt.NDArray[numpy.uint64]: ...

class Exons:
    def __init__(
        self,
        genes: Genes,
        gene_idx: npt.ArrayLike,
        start: npt.ArrayLike,
        end: npt.ArrayLike,
        annotated_start: npt.ArrayLike,
        annotated_end: npt.ArrayLike,
    ) -> None: ...
    def checksum(self) -> int: ...
    def index(
        self, gene_idx: npt.ArrayLike, start: npt.ArrayLike, end: npt.ArrayLike
    ) -> object: ...
    def is_denovo(self, exon_idx: npt.ArrayLike) -> npt.NDArray[numpy.bool_]: ...
    def potential_introns(self, make_simplified: bool) -> GeneIntrons: ...
    def __len__(self) -> int: ...
    @property
    def _parent_idx_end(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def _parent_idx_start(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def _parents(self) -> Genes: ...
    @property
    def annotated_end(self) -> npt.NDArray[numpy.int64]: ...
    @property
    def annotated_start(self) -> npt.NDArray[numpy.int64]: ...
    @property
    def end(self) -> npt.NDArray[numpy.int64]: ...
    @property
    def gene_idx(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def start(self) -> npt.NDArray[numpy.int64]: ...

class ExperimentStrandness:
    __members__: ClassVar[dict] = ...  # read-only
    FORWARD: ClassVar[ExperimentStrandness] = ...
    NONE: ClassVar[ExperimentStrandness] = ...
    REVERSE: ClassVar[ExperimentStrandness] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    def __setstate__(self, state: int) -> None: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class ExperimentThresholds:
    def __init__(
        self,
        minreads: int = ...,
        mindenovo: int = ...,
        minpos: int = ...,
        max_pctbins: float = ...,
        junction_acceptance_probability: float = ...,
        intron_acceptance_probability: float = ...,
    ) -> None: ...
    def intron_thresholds_generator(
        self, total_bins: int
    ) -> IntronThresholdsGenerator: ...
    @property
    def intron_acceptance_probability(self) -> float: ...
    @property
    def junction_acceptance_probability(self) -> float: ...
    @property
    def max_pctbins(self) -> float: ...
    @property
    def mindenovo(self) -> int: ...
    @property
    def minpos(self) -> int: ...
    @property
    def minreads(self) -> int: ...

class GFF3FeatureType:
    __members__: ClassVar[dict] = ...  # read-only
    ACCEPT_GENE: ClassVar[GFF3FeatureType] = ...
    ACCEPT_TRANSCRIPT: ClassVar[GFF3FeatureType] = ...
    EXON: ClassVar[GFF3FeatureType] = ...
    HARD_SKIP: ClassVar[GFF3FeatureType] = ...
    REJECT_OTHER: ClassVar[GFF3FeatureType] = ...
    REJECT_SILENT: ClassVar[GFF3FeatureType] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    def __setstate__(self, state: int) -> None: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class GFF3Types:
    def __init__(self) -> None: ...
    def items(self) -> Iterator: ...
    def __bool__(self) -> bool: ...
    def __contains__(self, arg0: str) -> bool: ...
    def __delitem__(self, arg0: str) -> None: ...
    def __getitem__(self, arg0: str) -> GFF3FeatureType: ...
    def __iter__(self) -> Iterator: ...
    def __len__(self) -> int: ...
    def __setitem__(self, arg0: str, arg1: GFF3FeatureType) -> None: ...

class GeneIntrons:
    def __init__(
        self,
        genes: Genes,
        gene_idx: npt.ArrayLike,
        start: npt.ArrayLike,
        end: npt.ArrayLike,
        denovo: npt.ArrayLike,
        passed_build: npt.ArrayLike,
        simplified: npt.ArrayLike,
    ) -> None: ...
    def _pass_all(self) -> None: ...
    def _simplify_all(self) -> None: ...
    def _unsimplify_all(self) -> None: ...
    def build_group(self) -> GroupIntronsGenerator: ...
    def checksum(self) -> int: ...
    def checksum_nodata(self) -> int: ...
    def connect_exons(self, exons: Exons) -> None: ...
    def filter_passed(
        self, keep_annotated: bool = ..., discard_denovo: bool = ...
    ) -> GeneIntrons: ...
    def index(
        self, gene_idx: npt.ArrayLike, start: npt.ArrayLike, end: npt.ArrayLike
    ) -> object: ...
    def update_flags_from(self, donor_introns: GeneIntrons) -> None: ...
    def __len__(self) -> int: ...
    @property
    def _parent_idx_end(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def _parent_idx_start(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def _parents(self) -> Genes: ...
    @property
    def connected_exons(self) -> Optional[Exons]: ...
    @property
    def denovo(self) -> npt.NDArray[numpy.bool_]: ...
    @property
    def end(self) -> npt.NDArray[numpy.int64]: ...
    @property
    def end_exon_idx(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def gene_idx(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def passed_build(self) -> npt.NDArray[numpy.bool_]: ...
    @property
    def simplified(self) -> npt.NDArray[numpy.bool_]: ...
    @property
    def start(self) -> npt.NDArray[numpy.int64]: ...
    @property
    def start_exon_idx(self) -> npt.NDArray[numpy.uint64]: ...

class GeneJunctions:
    def __init__(
        self,
        genes: Genes,
        gene_idx: npt.ArrayLike,
        start: npt.ArrayLike,
        end: npt.ArrayLike,
        denovo: npt.ArrayLike,
        passed_build: npt.ArrayLike,
        simplified: npt.ArrayLike,
    ) -> None: ...
    def _pass_all(self) -> None: ...
    def _simplify_all(self) -> None: ...
    def _unsimplify_all(self) -> None: ...
    def checksum(self) -> int: ...
    def checksum_nodata(self) -> int: ...
    def connect_exons(self, exons: Exons) -> None: ...
    def index(
        self, gene_idx: npt.ArrayLike, start: npt.ArrayLike, end: npt.ArrayLike
    ) -> object: ...
    def __len__(self) -> int: ...
    @property
    def _parent_idx_end(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def _parent_idx_start(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def _parents(self) -> Genes: ...
    @property
    def connected_exons(self) -> Optional[Exons]: ...
    @property
    def denovo(self) -> npt.NDArray[numpy.bool_]: ...
    @property
    def end(self) -> npt.NDArray[numpy.int64]: ...
    @property
    def end_exon_idx(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def gene_idx(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def passed_build(self) -> npt.NDArray[numpy.bool_]: ...
    @property
    def simplified(self) -> npt.NDArray[numpy.bool_]: ...
    @property
    def start(self) -> npt.NDArray[numpy.int64]: ...
    @property
    def start_exon_idx(self) -> npt.NDArray[numpy.uint64]: ...

class GeneStrandness:
    __members__: ClassVar[dict] = ...  # read-only
    __entries: ClassVar[dict] = ...
    ambiguous: ClassVar[GeneStrandness] = ...
    forward: ClassVar[GeneStrandness] = ...
    reverse: ClassVar[GeneStrandness] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    def __setstate__(self, state: int) -> None: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Genes:
    def __init__(
        self,
        contigs: Contigs,
        contig_idx: npt.ArrayLike,
        start: npt.ArrayLike,
        end: npt.ArrayLike,
        strand: npt.ArrayLike,
        gene_id: List[str],
        gene_name: List[str],
    ) -> None: ...
    def checksum(self) -> int: ...
    def __contains__(self, arg0: str) -> bool: ...
    def __getitem__(self, gene_id: str) -> int: ...
    def __len__(self) -> int: ...
    @property
    def _parent_idx_end(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def _parent_idx_start(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def _parents(self) -> Contigs: ...
    @property
    def contig_idx(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def end(self) -> npt.NDArray[numpy.int64]: ...
    @property
    def gene_id(self) -> List[str]: ...
    @property
    def gene_name(self) -> List[str]: ...
    @property
    def start(self) -> npt.NDArray[numpy.int64]: ...
    @property
    def strand(self) -> npt.NDArray: ...

class GroupIntronsGenerator:
    def __init__(self, gene_introns: GeneIntrons) -> None: ...
    def add_experiment(
        self, sj: SJIntronsBins, thresholds: ExperimentThresholds = ...
    ) -> None: ...
    def update_introns(self, min_experiments: float = ...) -> GeneIntrons: ...
    def __len__(self) -> int: ...
    @property
    def _introns(self) -> GeneIntrons: ...
    @property
    def num_experiments(self) -> int: ...
    @property
    def num_passed(self) -> npt.NDArray[numpy.uint64]: ...

class GroupJunctionsGenerator:
    def __init__(self, junctions: GeneJunctions, exons: Exons) -> None: ...
    def add_experiment(
        self,
        sjp: SJJunctionsBins,
        thresholds: ExperimentThresholds = ...,
        add_denovo: bool = ...,
    ) -> None: ...
    def pass_known_inplace(self, min_experiments: float = ...) -> GeneJunctions: ...
    def __len__(self) -> int: ...
    @property
    def num_denovo(self) -> int: ...
    @property
    def num_experiments(self) -> int: ...
    @property
    def num_known(self) -> int: ...

class IntronThresholdsGenerator:
    def __init__(self, *args, **kwargs) -> None: ...
    def __call__(self, intron_lengths: npt.ArrayLike) -> object: ...

class PassedJunctionsGenerator:
    def __init__(self, junctions: GeneJunctions) -> None: ...
    def add_group(
        self, group: GroupJunctionsGenerator, min_experiments: float = ...
    ) -> None: ...
    def get_passed(self, denovo_simplified: bool) -> GeneJunctions: ...
    def __len__(self) -> int: ...
    @property
    def num_denovo(self) -> int: ...
    @property
    def num_known(self) -> int: ...

class SJIntrons:
    def __init__(
        self,
        contigs: Contigs,
        contig_idx: npt.ArrayLike,
        start: npt.ArrayLike,
        end: npt.ArrayLike,
        strand: npt.ArrayLike,
        annotated: npt.ArrayLike,
    ) -> None: ...
    @staticmethod
    def from_exons_and_introns(
        exons: Exons, introns: GeneIntrons, stranded: bool
    ) -> SJIntrons: ...
    def __len__(self) -> int: ...
    @property
    def _parent_idx_end(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def _parent_idx_start(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def _parents(self) -> Contigs: ...
    @property
    def annotated(self) -> Any: ...
    @property
    def contig_idx(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def end(self) -> npt.NDArray[numpy.int64]: ...
    @property
    def start(self) -> npt.NDArray[numpy.int64]: ...
    @property
    def strand(self) -> npt.NDArray: ...

class SJIntronsBins:
    def __init__(
        self,
        sj: SJIntrons,
        bin_reads: npt.ArrayLike,
        bin_idx: npt.ArrayLike,
        _offsets: npt.ArrayLike,
        total_bins: int,
    ) -> None: ...
    @staticmethod
    def from_bam(
        bam_path: str,
        num_bins: int,
        exons: Exons,
        gene_introns: GeneIntrons,
        experiment_strandness: ExperimentStrandness,
        nthreads: int,
    ) -> SJIntronsBins: ...
    def numbins(
        self,
        region_idx: npt.ArrayLike,
        minreads: npt.ArrayLike,
    ) -> object: ...
    def numreads(
        self,
        region_idx: npt.ArrayLike,
        num_stacks: npt.ArrayLike,
    ) -> object: ...
    def numstacks(
        self,
        region_idx: npt.ArrayLike,
        pvalue_threshold: npt.ArrayLike = ...,
    ) -> object: ...
    def __len__(self) -> int: ...
    @property
    def _offsets(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def _regions(self) -> SJIntrons: ...
    @property
    def bin_idx(self) -> npt.NDArray[numpy.int32]: ...
    @property
    def bin_reads(self) -> npt.NDArray[numpy.float32]: ...
    @property
    def total_bins(self) -> int: ...

class SJJunctions:
    def __init__(
        self,
        contigs: Contigs,
        contig_idx: npt.ArrayLike,
        start: npt.ArrayLike,
        end: npt.ArrayLike,
        strand: npt.ArrayLike,
    ) -> None: ...
    def flip_strand(self) -> SJJunctions: ...
    def to_unstranded(self) -> SJJunctions: ...
    def __len__(self) -> int: ...
    @property
    def _contigs(self) -> Contigs: ...
    @property
    def _parent_idx_end(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def _parent_idx_start(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def _parents(self) -> Contigs: ...
    @property
    def contig_idx(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def end(self) -> npt.NDArray[numpy.int64]: ...
    @property
    def start(self) -> npt.NDArray[numpy.int64]: ...
    @property
    def strand(self) -> npt.NDArray: ...

class SJJunctionsBins:
    def __init__(
        self,
        sj: SJJunctions,
        bin_reads: npt.ArrayLike,
        bin_idx: npt.ArrayLike,
        _offsets: npt.ArrayLike,
        total_bins: int,
    ) -> None: ...
    @staticmethod
    def from_bam(
        bam_path: str, experiment_strandness: ExperimentStrandness, nthreads: int
    ) -> SJJunctionsBins: ...
    def numbins(self, region_idx: npt.ArrayLike, minreads: npt.ArrayLike) -> object: ...
    def numreads(
        self,
        region_idx: npt.ArrayLike,
        num_stacks: npt.ArrayLike,
    ) -> object: ...
    def numstacks(
        self,
        region_idx: npt.ArrayLike,
        pvalue_threshold: npt.ArrayLike = ...,
    ) -> object: ...
    def project_reads(self, arg0: SJJunctions, arg1: bool) -> SJJunctionsBins: ...
    def __len__(self) -> int: ...
    @property
    def _offsets(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def _regions(self) -> SJJunctions: ...
    @property
    def bin_idx(self) -> npt.NDArray[numpy.int32]: ...
    @property
    def bin_reads(self) -> npt.NDArray[numpy.int32]: ...
    @property
    def total_bins(self) -> int: ...

class SimplifierGroup:
    def __init__(self, exon_connections: ExonConnections) -> None: ...
    def add_experiment(
        self,
        sg_reads: SpliceGraphReads,
        simplify_min_psi: float = ...,
        simplify_minreads_annotated_junctions: float = ...,
        simplify_minreads_denovo_junctions: float = ...,
        simplify_minreads_introns: float = ...,
    ) -> None: ...
    def update_connections(self, simplifier_min_experiments: float = ...) -> None: ...
    @property
    def _exon_connections(self) -> ExonConnections: ...
    @property
    def introns_passed_dst(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def introns_passed_src(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def junctions_passed_dst(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def junctions_passed_src(self) -> npt.NDArray[numpy.uint64]: ...
    @property
    def num_experiments(self) -> int: ...

class SpliceGraph:
    def __init__(
        self,
        contigs: Contigs,
        genes: Genes,
        exons: Exons,
        junctions: GeneJunctions,
        introns: GeneIntrons,
    ) -> None: ...
    def close_to_annotated_exon(
        self, gene_idx: int, x: int, to_following: bool = ...
    ) -> bool: ...
    @staticmethod
    def from_gff3(
        gff3_path: str, process_ir: bool, gff3_types: Dict[str, GFF3FeatureType]
    ) -> SpliceGraph: ...
    @staticmethod
    def infer_exons(base_exons: Exons, junctions: GeneJunctions) -> Exons: ...
    def make_build_junctions(self) -> PassedJunctionsGenerator: ...
    def make_group_junctions(self) -> GroupJunctionsGenerator: ...
    def sj_introns(self, stranded: bool) -> SJIntrons: ...
    def sj_introns_from_bam(
        self,
        bam_path: str,
        num_bins: int,
        experiment_strandness: ExperimentStrandness = ...,
        nthreads: int = ...,
    ) -> SJIntronsBins: ...
    @property
    def _contigs(self) -> Contigs: ...
    @property
    def _exon_connections(self) -> ExonConnections: ...
    @property
    def _exons(self) -> Exons: ...
    @property
    def _genes(self) -> Genes: ...
    @property
    def _introns(self) -> GeneIntrons: ...
    @property
    def _junctions(self) -> GeneJunctions: ...

class SpliceGraphReads:
    def __init__(
        self,
        introns: GeneIntrons,
        junctions: GeneJunctions,
        introns_reads: npt.ArrayLike,
        junctions_reads: npt.ArrayLike,
    ) -> None: ...
    @staticmethod
    def from_sj(
        introns: GeneIntrons,
        junctions: GeneJunctions,
        sj_introns: SJIntronsBins,
        sj_junctions: SJJunctionsBins,
    ) -> SpliceGraphReads: ...
    @property
    def _introns(self) -> GeneIntrons: ...
    @property
    def _junctions(self) -> GeneJunctions: ...
    @property
    def introns_reads(self) -> npt.NDArray[numpy.float32]: ...
    @property
    def junctions_reads(self) -> npt.NDArray[numpy.float32]: ...

class VectorString:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: VectorString) -> None: ...
    @overload
    def __init__(self, arg0: Iterable) -> None: ...
    def append(self, x: str) -> None: ...
    def clear(self) -> None: ...
    def count(self, x: str) -> int: ...
    @overload
    def extend(self, L: VectorString) -> None: ...
    @overload
    def extend(self, L: Iterable) -> None: ...
    def insert(self, i: int, x: str) -> None: ...
    @overload
    def pop(self) -> str: ...
    @overload
    def pop(self, i: int) -> str: ...
    def remove(self, x: str) -> None: ...
    def __bool__(self) -> bool: ...
    def __contains__(self, x: str) -> bool: ...
    @overload
    def __delitem__(self, arg0: int) -> None: ...
    @overload
    def __delitem__(self, arg0: slice) -> None: ...
    @overload
    def __getitem__(self, s: slice) -> VectorString: ...
    @overload
    def __getitem__(self, arg0: int) -> str: ...
    def __iter__(self) -> Iterator: ...
    def __len__(self) -> int: ...
    @overload
    def __setitem__(self, arg0: int, arg1: str) -> None: ...
    @overload
    def __setitem__(self, arg0: slice, arg1: VectorString) -> None: ...

def _default_gff3_types() -> Dict[str, GFF3FeatureType]: ...
def set_seed(seed: int) -> None: ...
