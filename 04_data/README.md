# Private input contract

The CSV is a structure-only template with invented identifiers. It contains no sample measurement or real input location. Make an edited copy under `local/`, which is ignored by Git. Provide one row per tissue section with a unique `sample_id` and populated `spaceranger_dir`, `group` and `timepoint` fields. Group and timepoint labels may repeat.

Each Space Ranger path points to the directory containing `outs/`; the importer reads filtered HDF5 counts and spatial assets. See `METHODS.md` for required metadata. Do not add real counts, microscopy images, processed objects or the private manifest to the public repository.
