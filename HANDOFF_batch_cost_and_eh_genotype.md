# Handoff: ExpansionHunter Batch cost model and genotype failures

Written 2026-10-02 from a Claude (Cowork) session. Pick up here in Claude Code, inside this repo.
Scope for now: ExpansionHunter (EH) only. Locityper is parked until it has a real pilot.

All the analysis files are in `scratch/batch_cost_2026-10/` (gitignored). Nothing in tracked files has been changed.

## Goal

Know the per-sample compute and storage cost of the EH and Locityper Batch runs, so we can decide how many samples
we can afford, and see where the money goes. A side investigation found that most EH genotype tasks failed.

## Headline findings

1. **EH costs about 2.3 cents per sample (minicram plus genotype), about $12.4k for the full 535,662-sample cohort on-demand.**
   Per sample: Nearline retrieval 0.92c (about 0.92 GiB read), Nearline read requests 0.40c, minicram VM 0.79c,
   genotype VM 0.18c, other Cloud Storage 0.02c. Nearline reads are 57% of the cost, VMs 42%.
2. **A VM cost model built from Batch task timestamps reproduces billing.** Modeled vs billed per-shard cost for the nine
   labeled EH minicram shards is 1.000 (sum), and modeled E2 core-hours for Sep 29 are 560 vs 556 billed.
   The one fitted parameter is about 85 s of billed VM startup per VM that no task timestamp covers.
3. **EH genotype failed for 1,575 of 1,811 tasks on Sep 29 (236 succeeded).** Exit code 1, median 10 s. The pilot genotype
   (20 tasks) succeeded 20/20. See the next section: the likely cause is a shared work directory on packed VMs.
4. **repo `scripts/vwb_batch/cost.py` overestimates by about 1.6x** (3.6c vs 2.3c actual per sample). It prices each task at
   e2-standard-4 ($0.134/h, about 4.8x the real rate), but Batch actually picked `e2-highcpu-2` and packed 2 tasks per VM.
   It also omits Nearline Class B read requests (about 0.4c per sample). The two errors partly cancel.

## EH genotype failure: likely cause and fix

Evidence:
- `jobs.json` shows every full-scale EH job ran on `e2-highcpu-2`, `taskPack` 2, STANDARD provisioning. The pilot shards (s000)
  ran on `e2-standard-4` with `taskPack` 1.
- Genotype tasks on the same VM start within microseconds of each other. Across 866 VMs with two genotype tasks, both
  succeeded on only 2 VMs (about 15 expected if outcomes were independent at the overall 13% success rate), exactly one
  succeeded on 210, and neither on 654. Outcomes are strongly anti-correlated, which points to the two tasks interfering.
- `expansion_hunter/batch/submit_batch.py` hard-codes `WORKDIR = "/mnt/disks/eh"` (line 18) and uses it, with fixed file
  names, in all three runnables (localize, container, upload) and the container volume mount (`/mnt/disks/eh:/work`, line 318).
  Two tasks on one VM share `tasks.tsv`, `task.env`, `task.json`, `reads.cram`, `reference.fa`, `catalog.json`,
  `worker.log` and `resource_stats.tsv`. The localize and upload runnables source the shared `task.env`, so they can even
  act on the other task's sample.
- Minicram has the same layout but succeeded 1,831 of 1,835, probably because its output files are named by sample ID.
  Not verified.
- Success rate was flat (about 13%) across the launch window and shards, so a ramp or load effect does not explain it.
- Most failures happen before the Python worker starts: only about 121 samples have an `EH.worker.log` but no VCF, against
  1,575 failed tasks. The upload runnable has `alwaysRun: true`, so a task that reached Python would normally leave a log.

Not verified: the exact mechanism. Cloud Logging would settle it, but neither the pet service account nor Kiran's user has
Logs Viewer on `wb-steady-asparagus-6800` (PERMISSION_DENIED for both). Ask the workspace admin for `roles/logging.viewer`.

Proposed fix (cheapest first):
1. Give each task its own work directory, for example `WORKDIR=/mnt/disks/eh/t${BATCH_TASK_INDEX}`, in all three
   runnables and in the container volume (`"${WORKDIR}:/work"`). Keeps the 2-per-VM packing (the cheap config).
2. Or set `taskCountPerNode: 1`. Simple but roughly doubles VM cost (about +1c per sample, about +40% total).
3. Check `locityper/batch/submit_batch.py` for the same pattern before any full Locityper run. Unchecked.

Verify cheaply before scaling: rerun about 20 of the failed samples. Set `PARALLELISM = 20`, run the poll cell first,
set `DISPATCH_GENOTYPE = True`, and add `todo = todo[:20]` before the dispatch loop in
`notebooks/rw/expansion_hunter_04_dispatch.ipynb`. A fixed per-task dir should give about 100% success even at pack 2.

Separate logging bug: in `expansion_hunter/batch/genotype.py`, `main()` catches `Exception`, but input checks raise
`SystemExit(...)`, which is not an `Exception`. Those messages skip `traceback.print_exc()` and never reach `worker.log`
(Kiran's sample 1001792 log held only the resource-monitor line, `wall_sec=0.0`). Catch `BaseException`, or print the message
before raising. ExpansionHunter's own stderr is also not teed into the log (the subprocess writes to fd 2 directly).
Locityper saves a `.host.log`; EH has none.

## Cost model (what the analysis does)

- VM cost per task: (VM wall time covering its tasks + 85 s startup) / 3600 x VM $/h, split across tasks packed on that VM.
  VM $/h is cores x $0.021813 + GiB x $0.002926 + 53 GB boot disk x $0.10/GiB-month (E2 us-central1 on-demand, matching
  the billed SKU rates; `e2-highcpu-2` = 2 vCPU / 2 GiB).
- Nearline: billed Sep 29 totals divided by the 1,815 minicram samples that day: $16.70 retrieval for 1,669.97 GiB
  (0.92 GiB per sample) and $7.31 for 7,305,031 Class B requests. Assumes all Nearline traffic that day was EH minicram.
  Pilot `network_bytes` mean (0.93 GiB) agrees. Source CRAMs are in the AoU controlled bucket (requester pays).
- Billing days are Pacific time. The 20-sample pilot ran 06:17 to 06:40 UTC, so it was billed on Sep 28.
- Labeled billing: Batch VMs carry `batch-job-id` (job name) in billing, so the Console "Group by Labels" gives per-shard
  VM cost. Nearline cost cannot be attributed to a job.
- Unknown: about $6 of Sep 29 sat behind "More results" in the Console SKU table and was not itemized.
  Spot discount in the projection chart (65% off VM cost) is an assumption, not measured.

Sep 29 job IDs: minicram `eh-mc-eh-260929-pilot-s000..s010` (s001 to s010 submitted 07:42 to 07:43 UTC),
genotype `eh-gt-eh-260929-pilot-s000..s010` (s001 to s010 submitted 14:28 to 14:30 UTC). Run ID `eh-260929-pilot`,
ledger `gs://aou-lr-phase2-resources/batchRuns/expansion_hunter/runs/eh-260929-pilot/`.

## Locityper (parked)

The 20-sample pilot (`lt-260929-pilot`) used `locityper/smoke.bed` and `vcf_db.smoke.tar.gz`, so its timings (minicram 43 s on
`e2-standard-4`, genotype 61 s on `e2-highmem-2`) are not representative. Batch used larger machines than the tuned compute
in `run.json` (1 CPU / 2 GiB), worth checking how the Locityper job requests resources. Tasks ran one per VM, so the
85 s startup exceeds the work; packing several samples per task would help. Needs a full-BED pilot and the
`*.data_transfer_stats.tsv` files for Nearline bytes before any cost estimate.

## Next steps (suggested order)

1. Apply fix 1 above to `expansion_hunter/batch/submit_batch.py` (and check Locityper's copy); fix the `SystemExit` logging.
2. Rerun about 20 failed genotype samples at low concurrency and confirm success. Then rerun the remaining about 1,555.
3. Ask for Logs Viewer on the workspace project, and ask the Broad billing admins whether the detailed BigQuery billing
   export is enabled (not retroactive, so Sep 29 may not be in it).
4. Update `scripts/vwb_batch/cost.py` to match actuals: packed `e2-highcpu-2` rates, the 85 s startup, Class B requests
   (about 4,375 per GiB read), and use `batch_cost_report.ipynb` instead of the 20-sample estimate.
5. Recompute full-pipeline cost per finished sample after genotype is fixed (current figure uses only successful tasks).
6. Locityper full-BED pilot, then the same breakdown.

## Files in `scratch/batch_cost_2026-10/`

- `batch_cost_report.ipynb`: estimated and actual Batch cost notebook for a Workbench JupyterLab app (read-only gcloud;
  optional BigQuery billing export or a `batch_job_id,cost` CSV for actuals). Tested offline against the export; the live
  gcloud, repo-helper and BigQuery paths are untested. Could move to `notebooks/rw/`.
- `cost_model.py`: offline model and figure script. `python3 cost_model.py <export_dir> <out_dir>`.
- `gbhistory.sh`: exports Batch jobs/tasks, ledgers and object listings with `gcloud` (read-only).
- `cost_model_out/`: figures `fig1..fig5`, `summary.json`, `per_task.csv`, `per_sample_costs.csv`.
- `batch_export_2026-10-01.tgz`: raw export (jobs.json, tasks_*.json, ledgers, object listings). Do not commit.

Known limits: I could not reliably join failed Batch tasks to sample IDs (the task-index to sample mapping did not line
up with the output listing), so nothing here claims which samples fail. Cost per sample for genotype uses only successful
tasks; failed attempts add about $1.43 of VM time in total.
