# Setup: the two edgeR environments for Barracoda's `-V` option

Barracoda's edgeR results (p-values / enrichment calls) depend on the **edgeR
version**. `barracoda-2.0.sh` therefore has a `-V` option that selects which
R+edgeR environment runs:

- `-V stable` (default) — the reference R+edgeR that reproduces historical / web-server results.
- `-V latest` — a newer R+edgeR.

`-V` chooses the R interpreter from two paths configured in `barracoda-2.0.sh`
(`R_STABLE`, `R_LATEST`), each overridable by the env vars
`BARRACODA_R_STABLE` / `BARRACODA_R_LATEST`. Only the **R** steps use that
interpreter; bowtie2 / perl / GNU parallel still come from `PATH`, so the
`latest` environment only needs **R + edgeR + the Barracoda R packages**.

The script protects the operator from mis-wiring:
- if the selected R doesn't exist or **can't load edgeR**, the run **exits with a
  clear `[ERROR]`** up front (rather than failing deep inside the analysis);
- if `-V latest` resolves to the **same R as stable**, it prints a `[WARNING]`;
- every job logs `edgeR-version(-V)=… R=… edgeR=<version>` into its own log, so a
  saved result is always traceable to the edgeR that produced it.

This document is a self-contained task you can paste into a Claude Code session
**on the web-server host** to build and wire the two environments.

> ⚠️ Before handing this off, the `-V` feature, the mode-aware `smoke_test.sh`,
> and this doc must be **committed and pushed** to the branch the operator will
> check out. Confirm on that branch that `git log --oneline` shows the
> "edgeR version selection (-V)" commit and `grep -n edger_version barracoda-2.0.sh`
> returns hits.

---

## Paste the following into Claude on the server

> You are setting up two R+edgeR environments for Barracoda-2.0's `-V
> stable|latest` option, on this web-server host, and verifying them with the
> project's smoke test.
>
> **First, get the right code.** Work in a Barracoda-2.0 checkout of the branch
> that has the `-V` feature (branch `pin-web-versions` of
> github.com/SRHgroup/Barracoda-2.0). Verify BEFORE doing anything else:
> `grep -n edger_version barracoda-2.0.sh` must return hits and
> `ls smoke_test.sh tests/smoke/ docs/SETUP-edgeR-environments.md` must succeed.
> If not, `git fetch` and check out that branch (or the commit that adds `-V`).
> Do not proceed until this passes.
>
> Background you must respect:
> - `-V` only switches the **R** interpreter. bowtie2 / perl / GNU parallel come
>   from `PATH` (the existing deployment already has them), so the `latest` env
>   needs only R + edgeR + the R packages below.
> - `stable` must reproduce the **historical** results, so it must be the exact
>   R+edgeR the web server already uses — do **not** rebuild or modify it.
> - Default behaviour must not change: with no `-V`, `stable` is used, and if
>   `R_STABLE` is unset it falls back to `which R`.
>
> **Do NOT**: modify/rebuild the existing stable R, touch `/var/www/webface` or
> change how the web server currently runs, or use sudo/root. Only ADD the
> `latest` env and wire the two paths.
>
> ### Task 1 — identify the STABLE R (the reference)
> The reference is the R the webface actually used for historical runs — not
> necessarily an interactive shell's `which R` (a CGI environment can have a
> different PATH/user). Determine it authoritatively:
> - Grep an archived production run log for the echoed R path — every R step logs
>   `<path>/bin/R --vanilla --slave ...`. That `<path>/bin/R` is the reference.
>   (The production log we have showed `/opt/R-4.3.1/bin/R`; treat that as a
>   hypothesis to confirm, not a given.)
> - If unsure, inspect the webface's own environment/PATH as the **webface user**,
>   not your login shell.
>
> Record its versions and confirm it loads every package the pipeline actually uses:
> ```
> <stableR>/bin/Rscript -e 'cat(R.version.string,"\n");
>   for (p in c("edgeR","limma")) cat(p, as.character(packageVersion(p)), "\n");
>   for (p in c("tidyverse","openxlsx","edgeR","limma","squash","ggplot2","reshape2","scales","data.table"))
>     suppressMessages(stopifnot(require(p, character.only=TRUE)));
>   cat("all required R packages load OK\n")'
> ```
> (The pipeline uses **openxlsx**, not `xlsx`, and does not use `ggseqlogo` — do
> not treat those as required, and do not try to add them to the reference R.)
>
> ### Task 2 — build the LATEST env
> Newest R + latest edgeR + the pipeline's R deps, user-space, no root, via
> conda/mamba (install Miniconda into your home or a shared group dir if conda is
> absent):
> ```
> conda create -y -n barracoda_latest --override-channels -c conda-forge -c bioconda \
>   r-base bioconductor-edger bioconductor-limma \
>   r-openxlsx r-squash r-tidyverse r-ggplot2 r-dplyr r-data.table r-scales r-reshape2
> ```
> (`--override-channels` avoids the Anaconda ToS gate. Do NOT add `r-xlsx` /
> `openjdk` / `r-ggseqlogo` — the pipeline never loads them and `r-xlsx` drags in
> fragile rJava. bioconda's `bioconductor-edger` is built against a specific
> `r-base`; if the solve conflicts, let conda pick the compatible `r-base` rather
> than forcing the newest.)
> Record its R/edgeR/limma versions the same way as Task 1, and **assert it is
> actually newer** than stable — fail loudly if `packageVersion("edgeR")` in the
> latest env is **not greater** than in the stable env (otherwise `-V latest`
> would just reproduce `stable`).
>
> ### Task 3 — wire the two paths
> Set the two R-binary paths so `-V` finds them:
> - `<stableR>/bin/R` for stable, `<latest-env>/bin/R` for latest.
>
> Do it by editing the `R_STABLE=` / `R_LATEST=` lines in `barracoda-2.0.sh`'s
> config block, **or** via the `BARRACODA_R_STABLE` / `BARRACODA_R_LATEST` env
> vars in the environment where the webface invokes Barracoda. Recommended:
> explicitly pin `R_STABLE=<stableR>/bin/R` rather than relying on the `which R`
> fallback (a CGI PATH may resolve a different R and silently change `stable`).
>
> ### Task 4 — verify BOTH modes with the smoke test
> The smoke test runs the full pipeline on a tiny fixture under the selected env
> (needs bowtie2/perl/parallel on `PATH`, which the deployment has). It must be
> able to see the R paths, so **export them in the shell you run the test in**
> (the webface's env vars are not visible here):
> ```
> export BARRACODA_R_STABLE=<stableR>/bin/R
> export BARRACODA_R_LATEST=<latest-env>/bin/R
> bash smoke_test.sh stable    # must print: SMOKE TEST: PASS (stable)
> bash smoke_test.sh latest    # must print: SMOKE TEST: PASS (latest)
> ```
> Note: the smoke test proves each env **executes** the pipeline end-to-end (and
> that the selected R can load edgeR — the script errors otherwise). It does not
> compare numbers; `stable`'s fidelity comes from it being the unmodified
> reference binary. Sanity-check that each run's job log shows the expected edgeR
> version (`grep 'edgeR-version' <store>/.../log_file_*.txt`), and that stable and
> latest report **different** edgeR versions.
>
> ### Task 5 — report back
> The two R paths; the exact R/edgeR/limma versions of both envs; the newer-than-
> stable check result; how you wired them (script edit vs env vars); and both
> smoke-test results. Do not change any deployment default beyond making
> `-V latest` available.

---

## After the server is set up

- Web users keep getting `stable` (reference) results by default — nothing changes.
- To expose the choice on the website, the webface needs a control (e.g. a
  dropdown, default "stable") that passes `-V stable|latest` to `barracoda-2.0.sh`.
  That lives in the web front-end, separate from this repo.
- CLI users can run `barracoda-2.0.sh ... -V latest` directly once the env is wired.
- `test_environment.sh` checks the R on `PATH`, **not** the `-V`-selected R — use
  `smoke_test.sh stable` / `smoke_test.sh latest` to verify each environment.
