"""Check the arguments of every call to a function of the package.

    python3 tools/check_call_args.py [package root]

The functions defined at the top level of R/*.R are collected with their formal
arguments. Every call to one of them is then matched to its formals as R would match it
(exact names, partial names, then positions), in R/*.R, tests/testthat/*.R, inst/*/*.R,
tools/*.R, tools/*/*.R, the @examples sections of the roxygen comments and the R chunks of
vignettes/*.Rmd. A call to do.call() is checked in the same way when its argument list can
be resolved statically: a list() or c() of lists written in the call, modifyList(), a
subset X[setdiff(names(X), ...)], or a variable assigned one of these earlier in the same
file. Lists that cannot be resolved are counted and listed.

ERROR   an argument that matches no formal, more arguments than formals, or a formal
        without a default that the call does not supply; also a replacement given to
        local_mocked_bindings() or with_mocked_bindings() for a function of the package
        that lacks one of its formals and has no ..., so that the package's own calls
        of the function would fail while the mock is in place
WARN    an argument matched by a partial name, or an unnamed element in the argument list
        of do.call(), whose meaning depends on the order of the formals
The script prints every finding and a summary, and exits with status 1 if there is an
ERROR. It reads the files only and needs no R installation.
"""
import glob
import os
import re
import sys

# ---------------------------------------------------------------------------------------
# Tokenizer for R code


TOKEN = re.compile(r"""
    (?P<ws>[ \t\r\f]+)
  | (?P<nl>\n)
  | (?P<comment>\#[^\n]*)
  | (?P<rawstr>[rR]["'](?P<dash>-*)(?P<open>[\(\[\{]))
  | (?P<str>"(?:[^"\\]|\\.)*"|'(?:[^'\\]|\\.)*')
  | (?P<bq>`[^`]*`)
  | (?P<num>(?:0[xX][0-9a-fA-F]+|(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)[Li]?)
  | (?P<ident>(?:[A-Za-z]|\.(?![0-9]))[A-Za-z0-9._]*|\.)
  | (?P<op>%[^%\n]*%|<<-|->>|<-|->|:::|::|==|!=|<=|>=|&&|\|\||\|>|[-+*/^<>=!&|~?:$@\\])
  | (?P<punct>[()\[\]{},;])
""", re.VERBOSE | re.DOTALL)
CLOSE = {'(': ')', '[': ']', '{': '}'}


def tokenize(text):
    """List of (kind, value, line) without whitespace, comments and newlines."""
    out = []
    pos, line = 0, 1
    while pos < len(text):
        m = TOKEN.match(text, pos)
        if m is None:
            raise SyntaxError(f'cannot tokenize at line {line}: {text[pos:pos + 20]!r}')
        kind = m.lastgroup if m.lastgroup not in ('dash', 'open') else None
        if m.group('rawstr'):
            close = CLOSE[m.group('open')] + m.group('dash') + m.group('rawstr')[1]
            end = text.find(close, m.end())
            if end < 0:
                raise SyntaxError(f'unterminated raw string at line {line}')
            value = text[pos:end + len(close)]
            out.append(('str', value, line))
            line += value.count('\n')
            pos = end + len(close)
            continue
        value = m.group(0)
        if m.group('nl'):
            line += 1
        elif m.group('str') or m.group('bq'):
            out.append(('str' if m.group('str') else 'ident', value, line))
            line += value.count('\n')
        elif not (m.group('ws') or m.group('comment')):
            kind = [k for k in ('num', 'ident', 'op', 'punct') if m.group(k)][0]
            out.append((kind, value, line))
        pos = m.end()
    return out


def matching(tokens, i):
    """Index of the bracket closing the one at position i."""
    want = CLOSE[tokens[i][1]]
    depth = 0
    for j in range(i, len(tokens)):
        v = tokens[j][1]
        if tokens[j][0] == 'punct' and v in CLOSE:
            depth += 1
        elif tokens[j][0] == 'punct' and v in (')', ']', '}'):
            depth -= 1
            if depth == 0:
                if v != want:
                    raise SyntaxError(f'mismatched bracket at line {tokens[j][2]}')
                return j
    raise SyntaxError(f'unclosed bracket at line {tokens[i][2]}')


def split_args(tokens, i, j):
    """Arguments between the brackets at i and j, as lists of tokens."""
    args, cur, depth = [], [], 0
    for k in range(i + 1, j):
        kind, v, _ = tokens[k]
        if kind == 'punct' and v in CLOSE:
            depth += 1
        elif kind == 'punct' and v in (')', ']', '}'):
            depth -= 1
        if depth == 0 and kind == 'punct' and v == ',':
            args.append(cur)
            cur = []
        else:
            cur.append(tokens[k])
    if cur or args:
        args.append(cur)
    return args


def unquote(v):
    if v[:1] in ('"', "'", '`'):
        return v[1:-1]
    return v


def arg_name(arg):
    """Name of an argument written as name = value, or None."""
    if len(arg) >= 2 and arg[1] == ('op', '=', arg[1][2]) and arg[0][0] in ('ident', 'str'):
        return unquote(arg[0][1])
    return None


# ---------------------------------------------------------------------------------------
# Function definitions of the package


def definitions(files):
    """Top-level functions name -> list of (formal, has_default)."""
    defs = {}
    for f in files:
        tokens = tokenize(open(f, encoding='utf-8').read())
        depth = 0
        for k, (kind, v, _) in enumerate(tokens):
            if kind == 'punct' and v in CLOSE:
                depth += 1
            elif kind == 'punct' and v in (')', ']', '}'):
                depth -= 1
            elif (depth == 0 and kind == 'ident' and k + 3 < len(tokens)
                  and tokens[k + 1][1] in ('<-', '=') and tokens[k + 2][1] == 'function'
                  and tokens[k + 3][1] == '('):
                j = matching(tokens, k + 3)
                formals = []
                for a in split_args(tokens, k + 3, j):
                    if a:
                        formals.append((unquote(a[0][1]), len(a) > 1))
                defs[unquote(v)] = formals
    return defs


# ---------------------------------------------------------------------------------------
# Matching of a call to the formals, as done by R


def match_call(formals, args):
    """args is a list of (name or None, label). Returns (errors, warnings)."""
    errors, warnings = [], []
    names = [f for f, _ in formals]
    dots = '...' in names
    before = names[:names.index('...')] if dots else names
    matched = {}
    rest = []
    # Exact names
    for name, label in args:
        if name is not None and name in names and name != '...' and name not in matched:
            matched[name] = label
        else:
            rest.append((name, label))
    # Partial names, for the formals before ...
    rest2 = []
    for name, label in rest:
        if name is None:
            rest2.append((name, label))
            continue
        cand = [f for f in before if f.startswith(name) and f not in matched]
        if len(cand) == 1:
            matched[cand[0]] = label
            warnings.append(f'argument {name} matched partially to {cand[0]}')
        elif dots:
            matched.setdefault('...', label)
        else:
            errors.append(f'argument {name} matches no formal'
                          + (' uniquely' if len(cand) > 1 else ''))
    # Positions
    free = [f for f in before if f not in matched]
    for name, label in rest2:
        if free:
            matched[free.pop(0)] = label
        elif dots:
            matched.setdefault('...', label)
        else:
            errors.append('more arguments than formals')
            break
    for f, has_default in formals:
        if f != '...' and not has_default and f not in matched:
            errors.append(f'formal {f} without a default is not supplied')
    return errors, warnings


# ---------------------------------------------------------------------------------------
# Static resolution of the argument list of do.call()


def find_assignment(tokens, name, before):
    """Tokens of the last expression assigned to name before position before."""
    best = None
    for k in range(before):
        if (tokens[k] == ('ident', name, tokens[k][2]) and k + 1 < len(tokens)
                and tokens[k + 1][1] in ('<-', '=') and (k == 0 or tokens[k - 1][1] not in
                                                          ('$', '@', '(', ','))):
            best = k + 2
    if best is None:
        return None
    # The expression is a call: an identifier followed by its brackets
    if best + 1 < len(tokens) and tokens[best][0] == 'ident' and tokens[best + 1][1] == '(':
        return tokens[best:matching(tokens, best + 1) + 1]
    if tokens[best][0] == 'ident':
        return tokens[best:best + 1]
    return None


def resolve(expr, tokens, pos, depth=0):
    """List of (name or None) of the arguments an expression supplies, or None."""
    if depth > 10 or not expr:
        return None
    if len(expr) == 1 and expr[0][0] == 'ident':
        sub = find_assignment(tokens, expr[0][1], pos)
        return resolve(sub, tokens, pos, depth + 1) if sub else None
    if expr[0][0] != 'ident' or len(expr) < 3 or expr[1][1] != '(':
        return None
    j = matching(expr, 1)
    head = expr[0][1]
    if j == len(expr) - 1:
        parts = split_args(expr, 1, j)
        if head == 'list':
            return [arg_name(a) for a in parts if a]
        if head == 'c':
            out = []
            for a in parts:
                r = resolve_any(a, tokens, pos, depth + 1)
                if r is None:
                    r = scalar_element(a)
                if r is None:
                    return None
                out += r
            return out
        if head == 'modifyList' and len(parts) == 2:
            a = resolve(parts[0], tokens, pos, depth + 1)
            b = resolve(parts[1], tokens, pos, depth + 1)
            if a is None or b is None or None in a or None in b:
                return None
            return a + [n for n in b if n not in a]
        return None
    return None


def resolve_subset(expr, tokens, pos):
    """X[setdiff(names(X), drop)] with drop a string or c() of strings."""
    if not (len(expr) >= 4 and expr[0][0] == 'ident' and expr[1][1] == '['
            and matching(expr, 1) == len(expr) - 1):
        return None
    inner = expr[2:-1]
    if not (inner and inner[0][1] == 'setdiff' and inner[1][1] == '('):
        return None
    parts = split_args(inner, 1, matching(inner, 1))
    if len(parts) != 2:
        return None
    first = [t[1] for t in parts[0]]
    if first != ['names', '(', expr[0][1], ')']:
        return None
    drop = [unquote(t[1]) for t in parts[1] if t[0] == 'str']
    base = resolve(expr[:1], tokens, pos)
    if base is None:
        return None
    return [n for n in base if n not in drop]


def resolve_any(expr, tokens, pos, depth=0):
    r = resolve_subset(expr, tokens, pos)
    if r is not None:
        return r
    return resolve(expr, tokens, pos, depth)


SCALARS = {'TRUE', 'FALSE', 'NULL', 'NA', 'NA_integer_', 'NA_real_', 'NA_character_',
           'Inf', 'NaN'}


def scalar_element(arg):
    """[name or None] for an element of c() that is a single literal value, else None."""
    name = arg_name(arg)
    value = arg[2:] if name is not None else arg
    if value and value[0][1] == '-':
        value = value[1:]
    if len(value) == 1 and (value[0][0] in ('str', 'num') or value[0][1] in SCALARS):
        return [name]
    return None


# ---------------------------------------------------------------------------------------
# Sources of R code


def roxygen_examples(text):
    """Code of the @examples sections, with the other lines blanked to keep line numbers."""
    out, inside = [], False
    for line in text.split('\n'):
        m = re.match(r"\s*#'\s?(.*)", line)
        if m and re.match(r'@examples\b', m.group(1).strip()):
            inside = True
            out.append('')
            continue
        if not m or m.group(1).strip().startswith('@'):
            inside = False
        body = m.group(1) if (m and inside) else ''
        # The line opening a \donttest{ or \dontrun{ wrapper. Its closing brace is left in
        # place, which does not affect the calls
        if inside and re.match(r'\s*\\(donttest|dontrun|dontshow)\{\s*$', body):
            body = ''
        out.append(body)
    return '\n'.join(out)


def rmd_chunks(text):
    out, inside = [], False
    for line in text.split('\n'):
        if re.match(r'\s*```\s*\{r[ ,}]', line):
            inside = True
            out.append('')
        elif inside and re.match(r'\s*```\s*$', line):
            inside = False
            out.append('')
        else:
            out.append(line if inside else '')
    return '\n'.join(out)


def check_mocks(tokens, k, label, line, defs, report):
    """Formals of the replacements given to a call of local_mocked_bindings() or
    with_mocked_bindings() at position k, compared with those of the package."""
    j = matching(tokens, k + 1)
    for a in split_args(tokens, k + 1, j):
        name = arg_name(a)
        if name not in defs or len(a) < 4 or a[2][1] != 'function' or a[3][1] != '(':
            continue
        jf = matching(a, 3)
        mock = [unquote(f[0][1]) for f in split_args(a, 3, jf) if f]
        report['mocks'] += 1
        if '...' in mock:
            continue
        missing = [f for f, _ in defs[name] if f != '...' and f not in mock]
        if missing:
            report['error'].append(f'{label}:{line} mock of {name}(): lacks the formal(s) '
                                   + ', '.join(missing))


def check_text(text, label, defs, report):
    tokens = tokenize(text)
    for k, (kind, v, line) in enumerate(tokens):
        if kind != 'ident' or k + 1 >= len(tokens) or tokens[k + 1][1] != '(':
            continue
        prev = tokens[k - 1][1] if k > 0 else ''
        if prev in ('$', '@', 'function'):
            continue
        if unquote(v) in ('local_mocked_bindings', 'with_mocked_bindings'):
            check_mocks(tokens, k, label, line, defs, report)
            continue
        if prev in ('::', ':::') and (k < 2 or tokens[k - 2][1] != 'bbssr'):
            continue
        name = unquote(v)
        j = matching(tokens, k + 1)
        parts = split_args(tokens, k + 1, j)
        if name == 'do.call' and parts:
            target = parts[0]
            if target and target[-1][0] in ('ident', 'str'):
                fname = unquote(target[-1][1])
            else:
                continue
            if fname not in defs:
                continue
            if len(parts) < 2:
                continue
            args = resolve_any(parts[1], tokens, k)
            where = f'{label}:{line} do.call({fname})'
            if args is None:
                report['unresolved'].append(where)
                continue
            report['calls'] += 1
            labelled = [(n, i) for i, n in enumerate(args)]
            if any(n is None for n in args):
                report['warn'].append(f'{where}: unnamed element in the argument list')
            errs, warns = match_call(defs[fname], labelled)
            report['error'] += [f'{where}: {e}' for e in errs]
            report['warn'] += [f'{where}: {w}' for w in warns]
            continue
        if name not in defs:
            continue
        if any(t[1] == '...' for a in parts for t in a):
            report['dots'] += 1
            continue
        report['calls'] += 1
        labelled = [(arg_name(a), i) for i, a in enumerate(parts) if a]
        errs, warns = match_call(defs[name], labelled)
        where = f'{label}:{line} {name}()'
        report['error'] += [f'{where}: {e}' for e in errs]
        report['warn'] += [f'{where}: {w}' for w in warns]


def main(root):
    rfiles = sorted(glob.glob(os.path.join(root, 'R', '*.R')))
    defs = definitions(rfiles)
    report = {'calls': 0, 'dots': 0, 'mocks': 0, 'error': [], 'warn': [], 'unresolved': []}
    sources = rfiles + sorted(glob.glob(os.path.join(root, 'tests', 'testthat', '*.R'))) + \
        sorted(glob.glob(os.path.join(root, 'inst', '*', '*.R'))) + \
        sorted(glob.glob(os.path.join(root, 'tools', '*.R'))) + \
        sorted(glob.glob(os.path.join(root, 'tools', '*', '*.R')))
    for f in sources:
        text = open(f, encoding='utf-8').read()
        rel = os.path.relpath(f, root).replace(os.sep, '/')
        check_text(text, rel, defs, report)
        if f in rfiles:
            check_text(roxygen_examples(text), rel + ' (examples)', defs, report)
    rmd = sorted(glob.glob(os.path.join(root, 'vignettes', '*.Rmd')))
    for f in rmd:
        rel = os.path.relpath(f, root).replace(os.sep, '/')
        check_text(rmd_chunks(open(f, encoding='utf-8').read()), rel, defs, report)
    for e in report['error']:
        print('ERROR  ' + e)
    for w in report['warn']:
        print('WARN   ' + w)
    for u in report['unresolved']:
        print('UNRESOLVED  ' + u)
    print(f"{len(defs)} functions, {len(sources) + len(rmd)} files, {report['calls']} calls "
          f"checked, {report['dots']} calls passing ... skipped, "
          f"{len(report['unresolved'])} do.call() lists unresolved, "
          f"{report['mocks']} mocks checked, "
          f"{len(report['error'])} error(s), {len(report['warn'])} warning(s)")
    return 1 if report['error'] else 0


if __name__ == '__main__':
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))
    sys.exit(main(root))
