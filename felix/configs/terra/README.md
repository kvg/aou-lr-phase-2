# Terra-format input exports

Configs downloaded from (or uploadable to) the Terra workflow launcher. They
hold the same values as the repo-format configs one directory up, in Terra's
own JSON dialect: every non-string value is wrapped as a `"${...}"` string, so
arrays, numbers and booleans read as `"${[\"chr22\"]}"`, `"${5}"`, `"${true}"`
rather than as JSON types. Terra writes every declared input, so these files
also carry keys the repo-format configs leave at their WDL defaults.

| File | Workflow | Repo-format counterpart |
|---|---|---|
| `FelixRepeatChr22Limited.inputs.json` | `FelixGenome` | `../felix.chr22_repeat_dosage.inputs.limited.json` |
| `FelixPilotSVLimited.inputs.json` | `FelixPilot` | `../FelixPilotSVLimited.inputs.json` |
| `FelixPilotSVFull.inputs.json` | `FelixPilot` | `../FelixPilotSVFull.inputs.json` |
| `felixpilot.limited.inputs.json` | `FelixPilot` | `../felix.pilot.inputs.limited.json` |
| `felixpilot.full.inputs.json` | `FelixPilot` | `../felix.pilot.inputs.full.json` |

Each pair must be edited together — a key added to one is a key missing from
the other at submit time. `FelixRepeatChr22Limited` and its counterpart both
carry `ru_biallelic_arm`, `compare_ru_arms_script` and `ru_arm_p_threshold`
(the biallelic comparator arm, `felix/REPEAT_DOSAGE.md` §3.1).

Write small floats out in full: `"${0.00000005}"`, not `"${5e-08}"`. Terra
parses the `${...}` body as a WDL expression and rejects exponent notation with
"Error while parsing the expression. Offending symbol is on line 1 at position
0", which reads like a type error but is a lexer failure. This bites
`ru_arm_p_threshold` and `summarize_p_threshold`; both WDL defaults are already
`0.00000005`, so clearing the field in the UI is an equally good fix.

Three things these files do not validate, and Terra will not either:

- `ru_engine` has no validation in `FelixGenome.wdl`. Any value other than
  exactly `felix` or `rust` disables both repeat branches silently, so a
  trailing space turns the run into a no-op with no error.
- The `gs://` script paths must match what `scripts/stage_felix_scripts.sh`
  and `scripts/stage_tractor_scripts.sh` upload. A stale
  `write_admixed_dosage_vcf.py` in the bucket fails the comparator arm's
  writer task on an unrecognized `--coding` argument.
