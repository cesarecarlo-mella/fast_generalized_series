# Instructions for an agent running `cluster_eval` on the T30 cluster

You are operating a ready-made pipeline on a SLURM cluster. **Do not modify the numerical code**
(`src/`, `data/`, `exact/`). Your job:

1. Build the pipeline.
2. Submit the jobs.
3. Monitor them until the analysis has finished.
4. Fix operational problems as they come up.
5. Report the results.

## What the pipeline computes (context only)

It evaluates all 371 canonical master integrals of the "H" family, weights 1-6, by solving the canonical
differential equation numerically (Bulirsch-Stoer, boost odeint). It does this on 2016 lattice points in
two regions:

- `eucl`: x, y > 0;
- `scat`: scattering region, variables (u, v).

Each point is computed in `double` and in `dd_real` (double-double, the reference). `02_analyze.py`
compares them, checks both against the exact weight-1,2 solution, and writes tables and plots to
`results/`. Every task is independent: one SLURM array task processes one chunk of points.

## Hard facts about this cluster (learned the hard way)

1. **`/scratch` is node-local.** The login node's `/scratch` is not visible on compute nodes.
   Jobs started from `/scratch` fail with exit code **127** and leave **no log file**. The pipeline
   must live in **`$HOME`** (`/home/t30/...`, NFS, shared, no tight quota; the run needs about 2 GB).
2. **Execute bits get lost when files are copied.** Every SLURM task runs `01_run_chunk.sh`
   directly, so run `chmod +x *.sh` first. Otherwise tasks fail at once (exit code **126**,
   "Permission denied").
3. **`NCHUNK` must be identical for `00_setup.sh` and `slurm_submit.sh`.** Setup writes
   `points/<region>_<k>.txt` for k < NCHUNK, and submit launches array indices 0..NCHUNK-1. Always
   `export NCHUNK=...` once in the shell before running either.
4. **`newton4` is a login node, not a compute node.** `--nodelist=newton4` is invalid.
5. The `/space-nobackup` filesystem is over quota. Do not write there.

## Procedure

```bash
# 1. place the pipeline in $HOME
cd ~/cluster_eval                      # if absent: cp -r <source>/cluster_eval ~/cluster_eval
chmod +x *.sh
export NCHUNK=256

# 2. build (login node, ~1 min).  Success = ends with two "# 0.3 0.2 ..." lines and "setup done"
./00_setup.sh

# 3. sanity test of one chunk without the queue (~1 min). Success = "done .../out/eucl_double_0.out (8 points)"
./01_run_chunk.sh eucl double 0

# 4. submit: 4 arrays (eucl/scat x double/dd) of NCHUNK tasks each + 1 analysis job (dependency afterany)
./slurm_submit.sh

# 5. one minute later, verify that tasks really produce output
sleep 60; squeue -u $USER -h -o "%T %j" | sort | uniq -c; ls out | head
```

Step 5 should show `.out.tmp` (running) or `.out` (finished) files in `out/`. If tasks disappear
within seconds and `out/` stays empty, go to **Troubleshooting**.

## Monitoring

```bash
./status.sh       # queue states, points done / total and chunks done per region x precision, failed points, DONE flag
```

Poll every 5-10 minutes. Expected time per point at weight 6 on one core:

- double: 1-3 s;
- dd: 15-45 s.

With NCHUNK=256 a dd task takes about 3-6 minutes. The whole run is about 35-55 core-hours, and
the wall time depends on the queue. The job is finished when `status.sh` prints `DONE` and
`results/summary.md` exists.

## Troubleshooting (decision table)

| symptom | diagnosis command | cause / fix |
|---|---|---|
| tasks vanish, no log in `logs/` | `sacct -u $USER -S today -X --format=JobID%20,State,ExitCode,NodeList \| tail` | exit **127**: folder not on a shared FS, so move it to `$HOME` |
| exit **126** or "Permission denied" in the logs | `ls -l 01_run_chunk.sh` | `chmod +x *.sh`, then resubmit |
| "Illegal instruction" in the logs | `grep -l -i illegal logs/*.log` | rebuild: `ARCHFLAGS="-march=x86-64-v2" ./00_setup.sh` (NCHUNK exported), `rm -rf out`, resubmit |
| sbatch: invalid partition/account | - | `SBATCH_EXTRA="--partition=P --account=A" ./slurm_submit.sh` |
| tasks `TIMEOUT` | `sacct ... --format=State,Elapsed` | resubmit with a longer limit: `TIME_DD=08:00:00 ./slurm_submit.sh` (completed chunks are skipped) |
| `.out.tmp` left behind, or chunks missing after the arrays end | `./status.sh` | resubmit `./slurm_submit.sh`; only incomplete chunks rerun |
| many `FAILED max time` points | `cat out/*.out \| grep FAILED \| head` | raise the per-point limit: `MAX_TIME=7200`; delete those chunk files and resubmit |
| analysis job failed | `tail -n 30 logs/analyze.log` | usually python packages are missing: `pip install --user numpy mpmath matplotlib`, then `python3 02_analyze.py --jobs 8` |
| boost not found or link error at build | build output | `module load boost`, or `BOOST_INC=...`; for QD link errors build QD 2.3.24 and set `QD_DIR=...` |

To resubmit cleanly:

```bash
scancel -u $USER
rm -rf logs/*            # keep out/ to resume, or rm -rf out to restart
./slurm_submit.sh
```

Before running `scancel -u $USER`, check with `squeue -u $USER` that the user has no unrelated jobs
running. If they do, cancel by job name instead: `scancel -u $USER -n ce_eucl_dd` and so on.

**Do not:**

- run the pipeline from `/scratch`;
- change `NCHUNK` without rerunning `00_setup.sh` and deleting `out/`;
- edit `src/`;
- delete `results/ref_*.npz` (the exact-reference cache).

## What to report back

When `results/summary.md` exists, report:

1. **Status counts:** points done and failed per region and precision (from `./status.sh`).
2. **The full contents of `results/summary.md`:**
   - the timing table: median, mean, 90%, max and total core-hours, for double and dd per region;
   - double vs dd: digits per weight (worst point / 1% / median, % of points with ≥ 12, 10 and 8 digits);
   - double and dd vs exact at weights 1 and 2. The dd line should show ≥ 13 digits, which validates
     the reference.
3. **Anything unusual:** failed points (with coordinates), timeouts, warnings in `logs/analyze.log`.
4. **The location of the plots:**
   - `results/maps_eucl.png`, `maps_scat.png`;
   - `exact_*.png`;
   - `per_function_*.png`;
   - `timing.png`.

   They can be copied to the user's laptop with `scp -r newton4:~/cluster_eval/results .`.
