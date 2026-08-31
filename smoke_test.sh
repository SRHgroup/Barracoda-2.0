#!/usr/bin/env bash
# Barracoda smoke test - fast pre-commit gate.
#
# Runs the FULL pipeline on a tiny synthetic fixture (2 a-keys x 2 antigens,
# test_data/smoke/) and asserts it completes cleanly end-to-end. This catches
# not just non-zero exits but also the class of bug where barracoda-2.0.sh exits
# 0 while an R step (e.g. summarize) actually crashed - by checking the expected
# outputs exist and no error signatures appear in the R logs.
#
# Requires the barracoda toolchain on PATH (R+edgeR, bowtie2, bowtie2-build,
# perl, GNU parallel) - e.g. `conda activate barracoda`.
#
# Regenerate the fixture (rarely needed) with:
#   Rscript tests/smoke/make_smoke_data.R
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIX="$REPO/tests/smoke"

fail() { echo "SMOKE TEST: FAIL - $*" >&2; exit 1; }

# --- barracoda uses `declare -A`, which needs bash >= 4 ---
BASH_BIN="bash"
if [ "${BASH_VERSINFO:-0}" -lt 4 ]; then
  for b in /opt/homebrew/bin/bash /usr/local/bin/bash; do
    [ -x "$b" ] && BASH_BIN="$b" && break
  done
fi
"$BASH_BIN" -c '[ "${BASH_VERSINFO:-0}" -ge 4 ]' \
  || fail "need bash >= 4 (barracoda uses declare -A); found $("$BASH_BIN" --version | head -1)"

# --- required tools ---
for t in R bowtie2 bowtie2-build perl parallel; do
  command -v "$t" >/dev/null 2>&1 || fail "missing tool on PATH: $t (activate the barracoda env)"
done
[ -f "$FIX/reads.fastq" ] || fail "fixture missing - run: Rscript $FIX/make_smoke_data.R"

# --- run the pipeline into a throwaway store ---
STORE="$(mktemp -d)"
echo "SMOKE TEST: running pipeline (store: $STORE) ..."
"$BASH_BIN" "$REPO/barracoda-2.0.sh" \
  -f "$FIX/reads.fastq" \
  -m "$FIX/sampleID.xlsx" \
  -a "$FIX/annotations.xlsx" \
  -A "$FIX/tags.fasta" \
  -B GAAGTTCCAGCCAGCGTCACAGTTT \
  -C 6 \
  -D "$FIX/epitopeA.fasta" \
  -E GGTCAGCATCATTTCC \
  -F "$FIX/epitopeB.fasta" \
  -G 6 \
  -H GTTATCGGCTCGTTCACACTCGA \
  -s "$STORE" > "$STORE/run.log" 2>&1
rc=$?
[ "$rc" -eq 0 ] || { echo "--- run.log tail ---"; tail -25 "$STORE/run.log"; fail "pipeline exited $rc"; }

OUT="$(echo "$STORE"/store.barracoda_*/output)"
LOGS="$(echo "$STORE"/store.barracoda_*/files_intermediate/logs)"

# --- assertions: the experiment completed and produced its key outputs ---
[ -d "$OUT/experiment_RCC_17" ]                     || fail "experiment_RCC_17 output dir not created"
[ -s "$OUT/experiment_RCC_17/fold_change.xlsx" ]    || fail "fold_change.xlsx missing (summarize did not finish RCC_17)"
[ -s "$OUT/experiment_RCC_17/all-readcounts.xlsx" ] || fail "all-readcounts.xlsx missing"

# --- assertion: no hard-error signatures in any R log ---
if grep -qiE '\[ERROR\]|Execution halted|cannot be repeated|dim\(X\) must have|could not find' "$LOGS"/*.log 2>/dev/null; then
  echo "--- offending log lines ---"
  grep -niE '\[ERROR\]|Execution halted|cannot be repeated|dim\(X\) must have|could not find' "$LOGS"/*.log
  fail "error signature found in R logs"
fi

echo "SMOKE TEST: PASS - experiment_RCC_17 produced fold_change + readcounts, no errors"
rm -rf "$STORE"
