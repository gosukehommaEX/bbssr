"""Check that the test files carry the independent reference values.

    python3 tools/reference/check_test_constants.py [path/to/reference-values.txt]

Every value listed in reference-values.txt (written by reference_values.py) must appear
as a number in the named test file, to the relative tolerance given for its key. The
script prints PASS or FAIL for each key and exits with status 1 if any key fails.
"""
import os
import re
import sys

here = os.path.dirname(os.path.abspath(__file__))
pkg = os.path.dirname(os.path.dirname(here))
ref_file = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "reference-values.txt")
number = re.compile(r"(?<![\w.])-?\d+\.\d+(?:[eE][-+]?\d+)?|(?<![\w.])\d+(?![\w.])")

failed = 0
cache = {}
for line in open(ref_file, encoding="utf-8"):
    if line.startswith("#") or not line.strip():
        continue
    test_file, key, rtol, values = line.rstrip("\n").split("\t")
    rtol = float(rtol)
    if test_file not in cache:
        text = open(os.path.join(pkg, "tests", "testthat", test_file), encoding="utf-8").read()
        cache[test_file] = [float(x) for x in number.findall(text)]
    found = cache[test_file]
    missing = [v for v in map(float, values.split())
               if not any(abs(v - x) <= rtol * max(abs(v), 1e-300) for x in found)]
    if missing:
        failed += 1
        print(f"FAIL  {test_file}  {key}  missing: {' '.join('%.15g' % v for v in missing)}")
    else:
        print(f"PASS  {test_file}  {key}")
print(f"{failed} key(s) failed")
sys.exit(1 if failed else 0)
