#!/usr/bin/env bash
# Magma-level tests of the verifier (run by hand on Mordell; ~1 minute):
#   pipeline/tests/run_magma_tests.sh
# Runs verify.py on the fixtures into a scratch copy of the data directory, so the census is not touched,
# and checks the outcomes: the positives are accepted with the expected class, the negatives rejected.
set -euo pipefail
cd "$(dirname "$0")/../.."
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
cp -r pipeline "$TMP/pipeline"; mkdir -p "$TMP/data/curves" "$TMP/data/logs" "$TMP/data/rejected" "$TMP/submissions/inbox"
cp pipeline/tests/fixtures/*.json "$TMP/submissions/inbox/"
cd "$TMP"
python3 -B pipeline/verify.py --no-github --no-lmfdb --conductor-timeout 120 --require-claim 2> verify.err || true
cat verify.err
fail=0
check() { if grep -q "$2" verify.err; then echo "PASS: $1"; else echo "FAIL: $1 (expected '$2')"; fail=1; fi; }
check "positive [2,4] simple accepted"        "ACCEPTED as 2.4.a \[certified\] \[2,4\] class=simple"
check "positive [2,2,8] split over Q accepted" "ACCEPTED as 2.2.8.a \[certified\] \[2,2,8\] class=qsplit"
check "singular model rejected"                "negative_singular.json: rejected"
check "wrong claim rejected (--require-claim)" "negative_wrong_claim.json: rejected"
check "genus 1 rejected"                       "negative_genus1.json: rejected"
grep -q "singular" data/rejected/negative_singular.json && echo "PASS: singular reason recorded" || { echo "FAIL: singular reason"; fail=1; }
grep -q "genus 1" data/rejected/negative_genus1.json && echo "PASS: genus reason recorded" || { echo "FAIL: genus reason"; fail=1; }
exit $fail
