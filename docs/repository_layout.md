# Repository Layout

This document explains how to read the repository when you are new to the project.

## High-Level Structure

- `config/`: runtime configuration files.
- `data/`: small checked-in supporting data such as transition CSV inputs.
- `docs/`: user and maintainer documentation.
- `scripts/`: helper scripts for setup, documentation generation, figure generation, and exploratory workflows.
- `src/code/`: Julia source code for shared library utilities and executable jobs.
- `test/`: automated tests.
- `outputs/`: local generated products from pipeline runs.
- `tmp/`: local scratch space for temporary figures and review artefacts.

## What Is Source Code

These paths are the main source of truth for runtime behavior:

- `src/code/lib/DepolLib.jl`: shared helpers for config loading, file checks, path resolution, and REPL help.
- `src/code/jobs/`: the entrypoints that users execute directly.
- `src/code/instrumental_effect/`: the internal modules used by the instrumental-effect workflow.
- `config/default.toml`: the default local configuration template.

If you want to understand what the software does, start there.

## What Is Generated Locally

These paths are primarily outputs or working directories rather than maintained source:

- `outputs/`: task outputs from local runs.
- `tmp/`: scratch images and temporary review material.
- LaTeX auxiliary files in the repository root such as `paper.aux`, `paper.log`, and `paper.out`.

Those files are useful during analysis, but they are not the canonical implementation.

## How To Read `scripts/`

The `scripts/` directory contains a mix of stable helper scripts and analysis-oriented utilities.

Common categories are:

- environment helpers such as `precompile.jl`
- documentation helpers such as `generate_functions_doc.sh`
- figure-generation scripts used for paper assets
- one-off or exploratory scripts prefixed with `tmp_`

When in doubt, prefer starting from `src/code/jobs/` before diving into `scripts/`, because jobs define the supported workflows more clearly.

## How To Read `docs/`

Recommended reading order:

1. `getting_started.md`
2. `configuration.md`
3. `jobs.md`
4. `instrumental_effect.md`
5. `development.md`

Use `functions.md` only when you need a symbol-level index.

## Recommended Entry Points

Choose your starting point based on what you need:

- To run the pipeline: read `README.md`, then `docs/getting_started.md`.
- To understand outputs and dependencies: read `docs/jobs.md`.
- To work on the instrumental-effect analysis: read `docs/instrumental_effect.md`.
- To modify code safely: read `docs/development.md` and `test/runtests.jl`.

## Working Convention

The repository stays most readable when:

- reusable code goes in `src/code/`
- executable workflows stay in `src/code/jobs/`
- local run products stay in `outputs/`
- temporary experiments stay in `tmp/`
- documentation updates live alongside code changes when behavior changes

This separation is the main convention that keeps the project understandable for external users.
