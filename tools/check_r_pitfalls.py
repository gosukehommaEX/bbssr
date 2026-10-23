"""Detect R constructs that have caused errors in this package before.

    python3 tools/check_r_pitfalls.py [package root]

Two kinds of line are reported:

  pminmax  pmin() or pmax() whose first argument is a number, as in pmin(1, M). The
           result takes its attributes from the first argument, so a matrix in a later
           argument loses its dimensions; the matrix belongs first, as in pmin(M, 1).
  null     an assignment of NULL to an element of a list, as in args$N.max <- NULL,
           which removes the element instead of setting it to NULL. A list passed to
           do.call() then lacks the argument; args['N.max'] <- list(NULL) keeps it.

Every such line must be listed in tools/check_r_pitfalls-reviewed.txt (tab-separated:
kind, file, code of the line without leading and trailing spaces, reason). A line that is
not listed is NEW and makes the script fail, so that each new occurrence is looked at
once. Listed lines that no longer occur are reported as STALE. The report is printed and
written to tools/output/check-r-pitfalls.txt.
"""
import glob
import os
import re
import sys

root = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
patterns = {
    'pminmax': re.compile(r"\bp(?:min|max)\(\s*-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?L?\s*,"),
    'null': re.compile(r"^\s*[\w.]+(?:\$[\w.`]+|\[\[\s*['\"][^'\"]+['\"]\s*\]\])\s*<-\s*NULL\s*(?:#.*)?$"),
}


def code_part(line):
    """The line without a trailing comment, keeping # inside strings."""
    quote = None
    for i, c in enumerate(line):
        if quote:
            if c == '\\':
                continue
            if c == quote:
                quote = None
        elif c in '\'"`':
            quote = c
        elif c == '#':
            return line[:i]
    return line


files = []
for pat in ('R/*.R', 'tests/testthat/*.R', 'inst/**/*.R', 'vignettes/*.Rmd', 'tools/**/*.R'):
    files += glob.glob(os.path.join(root, pat), recursive=True)
found = []
for path in sorted(set(files)):
    rel = os.path.relpath(path, root).replace(os.sep, '/')
    for n, line in enumerate(open(path, encoding='utf-8'), 1):
        code = code_part(line.rstrip('\n'))
        if rel.endswith('.R') and line.lstrip().startswith("#'"):
            continue
        for kind, rx in patterns.items():
            if rx.search(code):
                found.append((kind, rel, n, line.strip()))

reviewed_file = os.path.join(root, 'tools', 'check_r_pitfalls-reviewed.txt')
reviewed = {}
for line in open(reviewed_file, encoding='utf-8'):
    if line.startswith('#') or not line.strip():
        continue
    kind, rel, code, reason = line.rstrip('\n').split('\t')
    reviewed[(kind, rel, code)] = reason
out = []
new = 0
seen = set()
for kind, rel, n, code in found:
    key = (kind, rel, code)
    seen.add(key)
    if key in reviewed:
        out.append(f'REVIEWED  {kind}  {rel}:{n}  {code}')
    else:
        new += 1
        out.append(f'NEW       {kind}  {rel}:{n}  {code}')
stale = [k for k in reviewed if k not in seen]
for kind, rel, code in stale:
    out.append(f'STALE     {kind}  {rel}  {code}')
out.append(f'{len(files)} files, {len(found)} line(s) found, {new} new, {len(stale)} stale')
text = '\n'.join(out) + '\n'
sys.stdout.write(text)
out_dir = os.path.join(root, 'tools', 'output')
if os.path.isdir(out_dir) and len(sys.argv) <= 1:
    open(os.path.join(out_dir, 'check-r-pitfalls.txt'), 'w', encoding='utf-8').write(text)
sys.exit(1 if new else 0)
