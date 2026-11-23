# CS7-8code
# human-gastrulation-cosmx

Analysis code for **"A molecular genetic cartography of human gastrulation"** (Cell, 2026).

The study profiles two 3D-reconstructed human gastrulating embryos with CosMx Whole Transcriptome (WTX, 18,934 genes) spatial single-cell imaging:

| label in paper | Carnegie stage | object name in code | older name in some files |
|---|---|---|---|
| HuE#1 | CS7 | `obj6` | "CS6" |
| HuE#2 | CS8 | `obj7` | "CS7" |

This repository contains **code only**. No count matrices, Seurat/AnnData objects, intermediate tables or figures are included (see `.gitignore`). The scripts are the final versions that produced the figures and tables of the manuscript. Relative to the authors' working copies, only paths, headers and comments were changed. Every module documents this per file in its `PROVENANCE.md`: source, sha256 of the original, and whether the script is verbatim or was reconstructed from the analysis log.

## Repository layout

```
human-gastrulation-cosmx/
├── cosmx/        CosMx processing, annotation/curation records, marker heatmaps, 3D/2D display,
│                 A–P signalling profiles, HOX group score, left–right comparison, QC
├── staging/      mouse-milestone staging (Fig 6A–B), Mittnenzweig validation (SF7),
│                 cross-species primitive-streak dynamics (Fig 6C–E, 6G), LEFTY/FGF detection (Fig 6F)
├── gastruloid/   benchmarking of gastruloid models against the CS7/CS8 embryos (Fig 6H, SF8)
├── aiggot/       AIGG-OT matching of HuE#1→HuE#2, sensitivity, benchmark, Stereo-seq replication (Fig 7A–D, SF9, Table S3)
├── generative/   generative model of the halfway (t = 0.5) embryo and transition analyses (Fig 7D–I, SF10B, Table S4,
│                 Supplementary Materials 3–4)
├── requirements.txt / environment.yml   Python environment
├── install_R_packages.R                 R packages
├── LICENSE (MIT) · CITATION.cff
```

Each module has a `README.md` with a Methods→script→figure table, the expected input layout, the run order and its known gaps. It also has a `PROVENANCE.md` with file-level provenance.

## STAR Methods → code → figure

| STAR Methods section | module / scripts | figure panels / tables |
|---|---|---|
| CosMx data processing, QC | `cosmx/preprocessing/` | SF1E–F, cell numbers in Table S1 |
| Cell-type annotation (RCTD + manual curation) | `cosmx/annotation/` (RCTD for HuE#2, endoderm re-clustering, adjudication, `CURATION_NOTES.md`) | labels used by all figures |
| Marker-gene heatmaps | `cosmx/heatmaps/` | SF2A, SF2B, SF3A |
| 3D reconstruction display | `cosmx/visualisation/build_supp_materials_1_2_and_sections.py` | Supplementary Materials 1–2, 2D section pages |
| Marker positivity, node analyses | `cosmx/visualisation/` | Fig 2C, Fig 2H (cell categories), Fig 2I, SF5F, SF4B-style gene panels |
| HOX group score along the A–P axis | `cosmx/hox/` | Fig 3L, 3M |
| Signalling-pathway profiles along the A–P axis | `cosmx/signalling/` | Fig 4A (specification), 4B, 4C, 4D; SF6A–E |
| Left–right comparison | `cosmx/left_right/` | SF5C, SF5D, SF5E |
| External datasets and pre-processing | `staging/00_inputs/` | inputs of Fig 6 |
| Milestone definition and staging windows | `staging/01_milestones/`, `staging/02_mittnenzweig_validation/` | Fig 6A, 6B; SF7; Table S2 |
| Cross-species primitive-streak dynamics | `staging/03_ps_dynamics/` | Fig 6C, 6D, 6E, 6G |
| LEFTY/FGF detection | `staging/04_lefty_fgf/` | Fig 6F |
| Benchmarking of gastruloid models | `gastruloid/01–05` | Fig 6H; SF8A, SF8B |
| AIGG-OT | `aiggot/lib/`, `aiggot/dataprep/`, `aiggot/matching/` | Fig 7A–D |
| Parameter sensitivity | `aiggot/sensitivity/` | SF9A–D; Table S3 |
| Benchmarking and replication | `aiggot/benchmark/`, `aiggot/stereoseq/`, `aiggot/figures/build_SF9.py` | SF9E–H; Table S3 |
| Generative model of intermediate states | `generative/latent/`, `generative/ai/`, `generative/out_domain/`, `generative/out_gencv/`, `generative/human/` | t = 0.5 cells in Fig 7D; Table S4 |
| Analysis of the reconstructed transition | `generative/analysis/`, `aiggot/matching/07_control_couplings.py` | Fig 7E–I; SF10B |
| Supplementary Materials 3–4 (animations) | `generative/animations/build_supp_material_3_4.py` | Supplementary Materials 3, 4 |

Panels assembled by hand in Adobe Illustrator: final panel arrangement of all figures, and the letters of Figure 6 and SF7. No code is provided for these.

## Data

Raw and processed CosMx data of HuE#1 and HuE#2 will be available from GEO (accession to be added on publication). The external datasets used are:

| dataset | accession | used in |
|---|---|---|
| Mouse gastrulation atlas (Pijuan-Sala et al., 2019) | ArrayExpress E-MTAB-6967 | staging, generative |
| Mouse embryo time course (Mittnenzweig et al., 2021) | GEO GSE169210 | staging (SF7), generative (training embryos) |
| Human CS7 gastrula, Smart-seq2 (Tyser et al., 2021) | ArrayExpress E-MTAB-9388 | staging |
| Human Stereo-seq CS7 / CS8 / CS9 | GSA-Human HRA006197 / HRA005567 / HRA007284 | staging, aiggot (SF9H) |
| Cynomolgus monkey embryos (Zhai et al., 2022) | GEO GSE193007 | staging |
| Gastruloid / embryo-model datasets | GEO GSE208369, GSE194214, GSE144897, GSE271604, GSE262081, GSE232861 | gastruloid |
| In-house 10x scRNA-seq (CS7, CS8) | unpublished | staging, annotation reference |
| Ensembl Compara one-to-one orthologues | Ensembl release 116 | staging |

See Table S1 and the Key Resources Table for details.

### Input layout and environment variables

All scripts read inputs from `$HGC_DATA` and write results to `$HGC_OUT`. Both default to `data/` and `results/`, relative to the working directory. The sub-paths below these roots are the ones used in the original runs and are listed in each module README.

| variable | meaning | default |
|---|---|---|
| `HGC_DATA` | input root | `data` |
| `HGC_OUT` | output root | `results` |
| `HGC_GEN_ROOT` | working tree of the generative model (`generative/`) | `.` |
| `HGC_FONTS` / `HGC_FONT_DIR` | folder with Arial (`cosmx`, `staging`, `gastruloid`) or Liberation Sans (`aiggot`, `generative`) `.ttf` files; optional | — |
| `HGC_TABLE_S2` | path of Table S2 (xlsx) for the staging/gastruloid scripts | `$HGC_DATA/Table_S2_milestones.xlsx` |
| `HGC_SM_INPUT` | input folder of the SM3/SM4 animation script | script folder |
| `STAGE`, `POPS`, `EXPORT`, `PREFIX` | per-run parameters of some `cosmx/signalling` and `cosmx/left_right` scripts (values given in `cosmx/README.md`) | — |

Several workflows pass intermediates from one script to the next through the same relative sub-paths. For a full re-run, set `HGC_OUT` equal to `HGC_DATA`, or copy the intermediates. The module READMEs give the details.

## Installation

```bash
conda env create -f environment.yml      # Python 3.11 environment
conda activate hgc
Rscript install_R_packages.R             # R 4.5 with Seurat 5, spacexr
```

GPU is optional. It is used by scVI, by UCE (33-layer model, installed separately from https://github.com/snap-stanford/UCE) and by the generator training in `generative/`. CellPhoneDB 5.0.1 (database v5.0.0) was used for Figure 5E/5G (see "Known gaps").

## Run order (overview)

1. `cosmx/preprocessing` → `cosmx/heatmaps`, `cosmx/visualisation`, `cosmx/hox`, `cosmx/signalling`, `cosmx/left_right` (independent of each other once the exports exist).
2. `aiggot/dataprep/01–04` (also writes the anatomical frames used by `cosmx/signalling`, `cosmx/left_right` and `generative/analysis`) → `aiggot/matching/05–07` → `aiggot/sensitivity`, `aiggot/benchmark`, `aiggot/stereoseq` → `aiggot/figures/build_SF9.py`.
3. `staging/00_inputs` → `01_milestones` → `02_mittnenzweig_validation` → `03_ps_dynamics` → `04_lefty_fgf`.
4. `gastruloid/01_preprocess` → `02_milestone_scoring` → `03_staging_vs_cosmx` → `04_denovo_annotation` → `05_celltype_comparison`.
5. `generative/latent` → `generative/ai` → `generative/out_domain` → `generative/out_gencv` → `generative/human` → `generative/analysis` → `generative/animations`.

The exact command lists are in each module README.

## Naming used inside the code

Most scripts were written before the final nomenclature was fixed:

| in code | in the paper |
|---|---|
| CS6 / CS7 (old file names), `obj6` / `obj7` | HuE#1 = CS7, HuE#2 = CS8 |
| HuE#2 `EPI`, `aECT`, `pECT` | EPI1, ECT, EPI2 |
| `ADE` | AVE |
| GG-OT, `ggot` | AIGG-OT |
| "Figure 6" (aiggot/generative comments), "Supplementary Figure 11" | Figure 7, Supplementary Figure 10 |
| `t05` | halfway (t = 0.5) state |
| "Heemskerk 2024", "Liu 2023", "Chen 2025" | Jo et al. (GSE271604), peri-gastruloid (GSE232861), GSE262081 |

## Known gaps

The following steps have no script in this repository. Either they were run in vendor software or interactively, or no final script could be recovered. Each module README gives the full list.

- CosMx acquisition, decoding, segmentation and export (AtoMx, vendor software).
- Serial-section 3D alignment with Spateo, and the <100-transcript cell filter applied before the objects used here were built.
- RCTD for HuE#1, and the interactive write-back of the manual curation (rules documented in `cosmx/annotation/CURATION_NOTES.md`).
- Final UMAPs (Fig 1E–F); dorsal/lateral 3D projections with KDE outlines (Fig 1H, 1J, 2L, 5D; SF4B, SF4E, SF5B); tables of Fig 2M, 3F, 3I; Fig 4E. Fig 4A is specified, but its exact five-population run was not saved.
- CellPhoneDB runs for Fig 5E and 5G.
- Mouse-atlas pre-processing (116,312 cells), the ortholog table, the Mittnenzweig metacell annotation table and the pre-audit milestone-call workbooks (`staging`).
- Figure 7 panel builders B/C/D/G/H/I, SF10A, the packaging of the t = 0.5 expression deliverable, and the writer of the bridge bundle (`generative`).

Several scripts were reconstructed from the execution log of the analysis sessions. These are marked "reconstructed" in `PROVENANCE.md` and have not been re-executed end-to-end. All Python files compile (`py_compile`), all R files parse, and all shell scripts pass `bash -n`.

## Software versions (Key Resources Table)

R 4.5.3; Seurat 5.5.1; SeuratObject 5.4.0; harmony 2.0.5; spacexr (RCTD). Python 3.11; scanpy 1.11.5; anndata 0.12.19; numpy 2.4.6; scipy 1.17.1; pandas 2.3.3; statsmodels 0.14.6; matplotlib 3.11.1; scikit-learn 1.9.1; leidenalg 0.12.0; harmonypy 2.0.0; Spateo 1.1.1; CellPhoneDB 5.0.1; OR-Tools 9.15.6755; POT 0.9.6.post1; moscot 0.5.2; ott-jax 0.6.0; jax 0.10.2; PyTorch 2.14.1; scvi-tools 1.4.2; Optuna 5.0.0; UCE (33-layer).

## Licence and citation

Code: MIT licence (see `LICENSE`). If you use this code, please cite the article (see `CITATION.cff`).
