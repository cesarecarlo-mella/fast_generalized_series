# How to run cluster_eval (T30 cluster, newton login node)

Step-by-step instructions for evaluating the H family (371 masters, weights 1-6) on 2016 lattice
points in the Euclidean and scattering regions. Every point is computed in double and in dd_real
(the reference). For background, see `README.md`.

---

## 0. Where to put the folder: `$HOME`, not `/scratch`

On this cluster **`/scratch` is a local disk on each machine.** The `/scratch` of the login node
(`newton4`) is not visible from the compute nodes. A job started from there cannot find the scripts,
and exits with code 127 and no log file.

`$HOME` (`/home/t30/...`, NFS) is shared by every node and has plenty of space. The full run writes
about 2 GB.

```bash
cp -r /path/to/cluster_eval ~/cluster_eval      # or unpack it there
cd ~/cluster_eval
chmod +x *.sh                                   # execute permission is lost when copying from a Mac
```

> Every SLURM task runs `01_run_chunk.sh` directly. Without `chmod +x` all tasks fail immediately
> with "Permission denied".

---

## 1. Build (once, on the login node, about 1 minute)

```bash
cd ~/cluster_eval
export NCHUNK=256          # number of array tasks per (region, precision); must be the SAME for setup and submit
./00_setup.sh
```

The end of the output should look like this:

```
eucl 2016 points, 256 chunks
scat 2016 points, 256 chunks
# 0.3 0.2 ...            <- smoke test, Euclidean
# 0.3 0.2 ...            <- smoke test, scattering
setup done
```

Possible problems:

| message | fix |
|---|---|
| `boost odeint not found` | `module load boost` (or `BOOST_INC=/dir/containing/boost ./00_setup.sh`) |
| link error involving `libqd.a` | build QD 2.3.24 (https://www.davidhbailey.com/dhbsoftware/) and run `QD_DIR=/its/prefix ./00_setup.sh` |
| a job later dies with "Illegal instruction" | the compute nodes have an older CPU: `ARCHFLAGS="-march=x86-64-v2" ./00_setup.sh` |

You only need to recompile when you change machines, `src/`, `WEIGHT` or `NCHUNK`.

---

## 2. Test one chunk by hand (optional, about 1 minute)

```bash
./01_run_chunk.sh eucl double 0
# -> done /home/.../cluster_eval/out/eucl_double_0.out (8 points)
```

---

## 3. Submit

```bash
cd ~/cluster_eval
export NCHUNK=256
./slurm_submit.sh
```

This submits four job arrays of `NCHUNK` tasks each: `eucl double`, `eucl dd`, `scat double` and
`scat dd`. It also submits an analysis job (`ce_analyze`), which waits until all four arrays have
finished.

Optional settings, passed as environment variables:

```bash
SBATCH_EXTRA="--partition=XYZ --account=ABC" ./slurm_submit.sh   # if a partition or account is required
TIME_DD=00:30:00 TIME_DOUBLE=00:10:00 ./slurm_submit.sh          # time limit per task (defaults: 4 h and 1 h)
REGIONS=scat ./slurm_submit.sh                                   # only one region
```

To start more jobs, set `NCHUNK` higher (up to 2016, i.e. one point per task), then rerun
`./00_setup.sh` and `rm -rf out`, then submit again.

Resubmitting is safe: chunks that are already complete are skipped.

---

## 4. Check progress

```bash
./status.sh                 # or keep it running: watch -n 30 ./status.sh
```

`status.sh` prints:

- the jobs in the queue, by state;
- the points done and chunks completed, for each region and precision;
- the number of points that failed;
- `DONE` once the analysis has finished.

The same information with plain commands:

```bash
squeue -u $USER -h -o "%T %j" | sort | uniq -c                      # queue
ls out | head                                                       # *.out.tmp = running, *.out = finished
sacct -u $USER -S today -X --format=State | sort | uniq -c          # COMPLETED / FAILED / TIMEOUT
tail -n 3 logs/eucl_dd_0.log                                        # one task's log
```

Typical times per point (weight 6, one core): **double 1-3 s**, **dd 15-45 s**. With `NCHUNK=256`
a dd task holds 8 points and takes 3-6 minutes. The total time depends mainly on how fast the queue
starts the tasks.

If tasks disappear after a few seconds:

1. Look at `logs/<region>_<prec>_<k>.log`.
2. If there is no log, check `sacct ... ExitCode`:
   - **127** means the folder is not on a shared disk (see step 0);
   - **126** means a missing `chmod +x`.

---

## 5. Results

Everything is written to `results/` when `ce_analyze` finishes.

| file | content |
|---|---|
| `summary.md` | tables: timing per region and precision; digits of double vs dd per weight; double and dd vs exact (w1, w2) |
| `maps_eucl.png`, `maps_scat.png` | double vs dd, one column per weight; rows worst and median component |
| `exact_eucl.png`, `exact_scat.png` | weight 2 vs the exact GPLs, double and dd |
| `per_function_*.png` | per component: median and worst point |
| `timing.png` | histogram of the time per point |
| `digits.npz` | all arrays (per point, weight and component), coordinates and times |

The Euclidean maps are plotted in (x, y), the scattering maps in (u, w).

If the analysis job failed, for example because python lacks `mpmath` or `matplotlib`, run it by
hand:

```bash
pip install --user numpy mpmath matplotlib     # once, if needed
python3 02_analyze.py --jobs 8
```

To copy the results to your laptop:

```bash
scp -r newton4:~/cluster_eval/results .         # run on the laptop
```

---

## 6. Start over

```bash
scancel -u $USER          # cancel everything still queued or running
rm -rf out logs/* results
./slurm_submit.sh
```
