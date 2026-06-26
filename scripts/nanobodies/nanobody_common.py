"""Shared nanobody plate-building logic.

This module holds everything that the nanobody plate scripts have in common:
loading parents from the lab sheet, matching recommended mutations by sequence,
building parent+mutant variants, laying them onto a plate (row- or column-wise,
96- or 384-well), running the per-variant synthesis pipeline, writing outputs,
and verifying the result. The thin CLI scripts just assemble a RunConfig and
call run_pipeline().
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from prosy.core import io
from prosy.core import layout as layout_mod
from prosy.core import platemap as platemap_mod
from prosy.core.cloning import Flanks
from prosy.core.constraints import ConstraintSet
from prosy.core.optimize import build_fragment
from prosy.core.sequence import (
    PointMutation,
    apply_mutation,
    translate,
    validate_protein,
)

# Default TSV filenames per set (override per task if needed).
DEFAULT_TSV_BY_SET = {
    "ADL1": "recommendations_ADL1.tsv",
    "KRU1": "recommendations_KRU1.tsv",
}


# --------------------------------------------------------------------------- #
# Config & variant model                                                      #
# --------------------------------------------------------------------------- #


@dataclass
class RunConfig:
    """All tunable knobs for one plate-building run."""

    flanks: Flanks
    species: str = "e_coli"
    avoid_enzymes: list[str] = field(default_factory=lambda: ["BsaI"])
    min_length: int = 300
    mutations_per_parent: int = 7
    new_id_start: int = 73
    new_id_prefix: str = "Nb"
    orientation: str = "column"          # "column" or "row"
    plate_size: int = 96                 # 96 or 384
    plate_map: str = "svg"               # "none", "svg", "png" or "both"

    @property
    def group_size(self) -> int:
        return self.mutations_per_parent + 1  # parent + mutants


@dataclass
class Variant:
    new_id: str | None          # new series ID for mutants; own ID for parents
    parent_id: str
    descriptive_id: str         # NbXX (parent) or NbXX_T27F (mutant)
    is_parent: bool
    set_name: str
    parent_tsv_id: str
    mutation: str | None
    score: int | None
    protein: str
    well: str = ""
    dna_coding: str = ""
    dna_final: str = ""

    @property
    def plate_name(self) -> str:
        if self.is_parent:
            return self.parent_id
        return f"{self.new_id} {self.descriptive_id}"


# --------------------------------------------------------------------------- #
# Loading & matching                                                          #
# --------------------------------------------------------------------------- #


def load_parents(
    xlsx_path: Path,
    parent_order: list[str],
    *,
    sheet: str = "Nanobodies",
    id_col: str = "Nanobody ID",
    set_col: str = "Set name",
    seq_col: str = "AA Sequence",
    alias_col: str = "Alias",
) -> dict[str, dict]:
    """Return {NbID: {set, alias, protein}} for the requested parents."""
    rows = io.read_xlsx_sheet(xlsx_path, sheet=sheet)
    wanted = set(parent_order)
    by_id: dict[str, dict] = {}
    for r in rows:
        nid = r.get(id_col)
        if nid in wanted:
            by_id[nid] = {
                "set": (r.get(set_col) or "").strip(),
                "alias": r.get(alias_col),
                "protein": validate_protein(str(r.get(seq_col)).strip()),
            }
    missing = [p for p in parent_order if p not in by_id]
    if missing:
        raise SystemExit(f"Parents not found in sheet {sheet!r}: {missing}")
    return by_id


def load_recommendations(folder: Path, tsv_by_set: dict[str, str]) -> dict:
    """Return {set_name: {full_protein_seq: [mutation rows]}}."""
    out: dict[str, dict[str, list[dict]]] = {}
    for set_name, fname in tsv_by_set.items():
        rows = io.read_delimited(folder / fname, delimiter="\t")
        by_seq: dict[str, list[dict]] = {}
        for r in rows:
            by_seq.setdefault(str(r["sequence"]).strip(), []).append(r)
        out[set_name] = by_seq
    return out


def top_mutations(rows: list[dict], n: int) -> list[dict]:
    """Top-n mutation rows by score (desc), tie-broken by position (asc)."""
    return sorted(rows, key=lambda m: (-int(m["score"]), int(m["position"])))[:n]


# --------------------------------------------------------------------------- #
# Build & lay out variants                                                     #
# --------------------------------------------------------------------------- #


def build_variants(
    parent_order: list[str], parents: dict, recs: dict, cfg: RunConfig
) -> list[Variant]:
    """Parent + mutant variants, grouped by parent in plate-fill order."""
    variants: list[Variant] = []
    next_new_id = cfg.new_id_start

    for parent_id in parent_order:
        meta = parents[parent_id]
        set_name = meta["set"]
        protein = meta["protein"]
        by_seq = recs.get(set_name)
        if by_seq is None:
            raise SystemExit(f"{parent_id}: no recommendations loaded for set {set_name!r}.")
        if protein not in by_seq:
            raise SystemExit(f"{parent_id}: parent sequence not found in {set_name} TSV.")
        rows = by_seq[protein]
        tsv_id = str(rows[0]["seq_id"])
        muts = top_mutations(rows, cfg.mutations_per_parent)
        if len(muts) < cfg.mutations_per_parent:
            raise SystemExit(
                f"{parent_id}: only {len(muts)} mutations available, "
                f"need {cfg.mutations_per_parent}."
            )

        variants.append(Variant(
            new_id=parent_id, parent_id=parent_id, descriptive_id=parent_id,
            is_parent=True, set_name=set_name, parent_tsv_id=tsv_id,
            mutation=None, score=None, protein=protein,
        ))
        for m in muts:
            mut = PointMutation(wt=m["wt"], position=int(m["position"]), mut=m["mut"])
            mutated = apply_mutation(protein, mut, check_wt=True)
            variants.append(Variant(
                new_id=f"{cfg.new_id_prefix}{next_new_id}", parent_id=parent_id,
                descriptive_id=f"{parent_id}_{mut}", is_parent=False,
                set_name=set_name, parent_tsv_id=tsv_id, mutation=str(mut),
                score=int(m["score"]), protein=mutated,
            ))
            next_new_id += 1
    return variants


def assign_wells(variants: list[Variant], parent_order: list[str], cfg: RunConfig) -> None:
    """Map each parent block onto one group of wells via the layout engine."""
    groups = layout_mod.layout_groups(
        n_groups=len(parent_order), group_size=cfg.group_size,
        plate=cfg.plate_size, orientation=cfg.orientation,
    )
    gi = 0
    i = 0
    # variants are emitted parent-block by parent-block, each of group_size.
    for _ in parent_order:
        block = variants[i : i + cfg.group_size]
        wells = groups[gi].wells
        for v, well in zip(block, wells):
            v.well = well
        i += cfg.group_size
        gi += 1


# --------------------------------------------------------------------------- #
# Synthesis pipeline                                                           #
# --------------------------------------------------------------------------- #


def optimize_all(variants: list[Variant], backend, cfg: RunConfig, *, seed: int) -> None:
    constraints = ConstraintSet(avoid_enzymes=cfg.avoid_enzymes)
    pad_constraints = ConstraintSet(
        avoid_enzymes=cfg.avoid_enzymes, max_homopolymer=4,
        forbid_low_complexity=True, gc_bounds=(0.30, 0.70),
    )
    for i, v in enumerate(variants):
        frag = build_fragment(
            v.protein, species=cfg.species, constraints=constraints,
            flanks=cfg.flanks, min_length=cfg.min_length,
            pad_constraints=pad_constraints, backend=backend, seed=seed + i,
        )
        v.dna_coding = frag.coding
        v.dna_final = frag.final


# --------------------------------------------------------------------------- #
# Output & verification                                                        #
# --------------------------------------------------------------------------- #


def write_outputs(variants: list[Variant], out_dir: Path, stamp: str, cfg: RunConfig) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    flanks = cfg.flanks

    plate_records = [
        {"well": v.well, "name": v.plate_name, "sequence": v.dna_final} for v in variants
    ]
    xlsx_path = out_dir / f"{stamp}_nanobody_plate.xlsx"
    io.write_plate_xlsx(plate_records, xlsx_path)
    xls_path = io.convert_to_xls(xlsx_path)

    design_rows = [
        {
            "Name": v.descriptive_id, "Sequence": v.protein,
            "5' nucleotides": flanks.five_prime, "3' nucleotides": flanks.three_prime,
        }
        for v in variants
    ]
    design_path = out_dir / f"{stamp}_design_seqs_nanobodies.csv"
    io.write_csv(design_rows, design_path, io.DESIGN_COLUMNS)

    mapping_rows = [
        {
            "Well": v.well, "NbID": v.new_id, "NbID_parent": v.descriptive_id,
            "Parent": v.parent_id, "Set": v.set_name, "Parent_TSV_ID": v.parent_tsv_id,
            "Mutation": v.mutation or "", "Score": "" if v.score is None else v.score,
            "Type": "parent" if v.is_parent else "mutant", "Protein": v.protein,
            "Protein_len": len(v.protein), "Coding_DNA": v.dna_coding,
            "Final_DNA": v.dna_final, "Final_len": len(v.dna_final),
        }
        for v in variants
    ]
    mapping_path = out_dir / f"{stamp}_nanobody_mapping.csv"
    io.write_csv(mapping_rows, mapping_path, list(mapping_rows[0].keys()))

    result = {"plate_xlsx": xlsx_path, "plate_xls": xls_path,
              "design": design_path, "mapping": mapping_path,
              "map_svg": None, "map_png": None}

    # Plate map (colour-coded by set; parents emphasized).
    if cfg.plate_map != "none":
        cells = {
            v.well: platemap_mod.Cell(
                top=v.new_id if not v.is_parent else v.parent_id,
                bottom=(v.mutation or "") if not v.is_parent else "(parent)",
                category=v.set_name,
                emphasize=v.is_parent,
            )
            for v in variants
        }
        title = f"{stamp} nanobody plate ({cfg.orientation}-major, {cfg.plate_size}-well)"
        if cfg.plate_map in ("svg", "both"):
            result["map_svg"] = platemap_mod.write_svg(
                cells, out_dir / f"{stamp}_nanobody_platemap.svg",
                plate=cfg.plate_size, title=title)
        if cfg.plate_map in ("png", "both"):
            result["map_png"] = platemap_mod.write_png(
                cells, out_dir / f"{stamp}_nanobody_platemap.png",
                plate=cfg.plate_size, title=title)

    return result


def verify(variants: list[Variant], cfg: RunConfig) -> list[str]:
    """Return a list of problems; empty means all checks pass."""
    from prosy.core import enzymes

    problems: list[str] = []
    avoided = enzymes.patterns_for(cfg.avoid_enzymes)
    wells_seen: set[str] = set()
    new_ids: list[str] = []

    for v in variants:
        if translate(v.dna_coding) != v.protein:
            problems.append(f"{v.descriptive_id}: coding DNA does not translate to protein.")
        cd = v.dna_coding.upper()
        for pat in avoided:
            if pat in cd:
                problems.append(f"{v.descriptive_id}: avoided site {pat} in coding region.")
        f5, f3 = cfg.flanks.five_prime.upper(), cfg.flanks.three_prime.upper()
        if f5 not in v.dna_final.upper() or f3 not in v.dna_final.upper():
            problems.append(f"{v.descriptive_id}: flanks not intact in final sequence.")
        if cd not in v.dna_final.upper():
            problems.append(f"{v.descriptive_id}: coding region not preserved in final.")
        if len(v.dna_final) < cfg.min_length:
            problems.append(f"{v.descriptive_id}: final length {len(v.dna_final)} < {cfg.min_length}.")
        if v.well in wells_seen:
            problems.append(f"{v.descriptive_id}: duplicate well {v.well}.")
        wells_seen.add(v.well)
        if not v.is_parent:
            new_ids.append(v.new_id)

    if len(wells_seen) != len(variants):
        problems.append(f"Expected {len(variants)} unique wells, got {len(wells_seen)}.")
    expected = [f"{cfg.new_id_prefix}{cfg.new_id_start + i}" for i in range(len(new_ids))]
    if new_ids != expected:
        problems.append(
            f"New ID series not contiguous from {cfg.new_id_prefix}{cfg.new_id_start}."
        )
    return problems


# --------------------------------------------------------------------------- #
# Orchestration                                                                #
# --------------------------------------------------------------------------- #


def run_pipeline(
    *,
    data_dir: Path,
    xlsx_name: str,
    parent_order: list[str],
    tsv_by_set: dict[str, str],
    stamp: str,
    backend,
    cfg: RunConfig,
    seed: int,
) -> tuple[list[Variant], list[str], dict]:
    """Full run: load -> build -> lay out -> optimize -> verify -> write."""
    xlsx = data_dir / xlsx_name
    parents = load_parents(xlsx, parent_order)
    recs = load_recommendations(data_dir, tsv_by_set)
    variants = build_variants(parent_order, parents, recs, cfg)
    assign_wells(variants, parent_order, cfg)
    optimize_all(variants, backend, cfg, seed=seed)
    problems = verify(variants, cfg)
    out = write_outputs(variants, data_dir, stamp, cfg)
    return variants, problems, out
