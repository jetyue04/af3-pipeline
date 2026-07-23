# Original AlphaFold 3 input template (verbatim, unmodified)

Pulled directly from the official DeepMind AlphaFold 3 docs (`docs/input.md`,
`google-deepmind/alphafold3` repo) for cross-checking against the filled-in
templates in this folder. Copied as-is — no edits, no placeholder
substitutions.

Note from the source doc: "AlphaFold 3 won't run this input out of the box
as it abbreviates certain fields and the sequences are not biologically
meaningful." The `...` values are the doc's own shorthand for "omitted for
brevity," not something we introduced.

**Cluster caveat (learned the hard way):** this example declares
`"version": 4`, and the docs describe versions up to 4 — but this cluster's
installed AF3 build (`AF3_v3.0.1.sif`) rejects anything above `2` at
runtime. Every actual input in [../inputs](../inputs) uses `"version": 2`
and skips the v4-only `description` field. Kept this example unedited
anyway since it's meant as a reference for the full public schema, not
something to run as-is.

## Top-level structure

```json
{
  "name": "Job name goes here",
  "modelSeeds": [1, 2],  # At least one seed required.
  "sequences": [
    {"protein": {...}},
    {"rna": {...}},
    {"dna": {...}},
    {"ligand": {...}}
  ],
  "bondedAtomPairs": [...],  # Optional.
  "userCCD": "...",  # Optional, mutually exclusive with userCCDPath.
  "userCCDPath": "...",  # Optional, mutually exclusive with userCCD.
  "dialect": "alphafold3",  # Required.
  "version": 4  # Required.
}
```

## Full example (all field types)

```json
{
  "name": "Hello fold",
  "modelSeeds": [10, 42],
  "sequences": [
    {
      "protein": {
        "id": "A",
        "sequence": "PVLSCGEWQL",
        "modifications": [
          {"ptmType": "HY3", "ptmPosition": 1},
          {"ptmType": "P1L", "ptmPosition": 5}
        ],
        "description": "10-residue protein with 2 modifications",
        "unpairedMsa": ...,
        "pairedMsa": ""
      }
    },
    {
      "protein": {
        "id": "B",
        "sequence": "RPACQLW",
        "templates": [
          {
            "mmcif": ...,
            "queryIndices": [0, 1, 2, 4, 5, 6],
            "templateIndices": [0, 1, 2, 3, 4, 8]
          }
        ]
      }
    },
    {
      "dna": {
        "id": "C",
        "sequence": "GACCTCT",
        "modifications": [
          {"modificationType": "6OG", "basePosition": 1},
          {"modificationType": "6MA", "basePosition": 2}
        ]
      }
    },
    {
      "rna": {
        "id": "E",
        "sequence": "AGCU",
        "modifications": [
          {"modificationType": "2MG", "basePosition": 1},
          {"modificationType": "5MC", "basePosition": 4}
        ],
        "unpairedMsa": ...
      }
    },
    {
      "ligand": {
        "id": ["F", "G", "H"],
        "ccdCodes": ["ATP"]
      }
    },
    {
      "ligand": {
        "id": "I",
        "ccdCodes": ["NAG", "FUC"]
      }
    },
    {
      "ligand": {
        "id": "Z",
        "smiles": "CC(=O)OC1C[NH+]2CCC1CC2"
      }
    }
  ],
  "bondedAtomPairs": [
    [["A", 1, "CA"], ["G", 1, "CHA"]],
    [["I", 1, "O6"], ["I", 2, "C1"]]
  ],
  "userCCD": ...,
  "dialect": "alphafold3",
  "version": 4
}
```

Source: https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md

---

**Re: `/opt/apps/community/alphafold3/AF3_Input/alphafold_input.json`** — that
path is on the SLURM cluster (it's what [jobs/af3-test.job](../../jobs/af3-test.job)
copies from), not reachable from this machine. If that's meant to be the
"original template" for cross-checking, `cat` or `scp` it over and I'll drop
the real one in here alongside this public-docs version.
