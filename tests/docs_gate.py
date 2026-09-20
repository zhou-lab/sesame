#!/usr/bin/env python3
"""The docs gate: check the docs against the binary that ships.

Two modes, because the two questions cost very different things to answer.

  docs_gate.py --consistency     no data, no network, under a second.
      Does the page say things this binary makes true? Versions agree across
      include/sesame.h, conda-recipe/meta.yaml and `sesame version`; every
      subcommand the docs name exists; every flag they quote appears in that
      subcommand's -h; the qc metric count matches SESAME_QC_FIELDS. This runs
      in `make test` and in CI.

  docs_gate.py                   the documented-workflow gate.
      Runs every docs/examples/*.sh against this checkout's binary, in a
      sandbox, as a reader would. Needs a populated $YAME_DATA_HOME and test
      IDATs at $SESAME_TEST_IDATS, so it is `make test-docs` and never CI.

Why both: a doc review in 2026-09 ran the page verbatim and found nine printed
lines that fail as printed -- a deprecated flag, a wrong probe name, a metric
count two short, four plot commands missing a required argument. Every one of
them would have been caught by one of these two modes.
"""
import os, re, subprocess, sys, glob, shutil, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.environ.get("SESAME_BIN") or os.path.join(ROOT, "sesame")
DOCS = [os.path.join(ROOT, p) for p in
        ("README.md", "docs/index.html", "docs/llms.txt")]

fails = []
def check(ok, msg):
    print(("  ok   " if ok else "  FAIL ") + msg)
    if not ok: fails.append(msg)

def run(*args):
    return subprocess.run([BIN, *args], capture_output=True, text=True).stdout + \
           subprocess.run([BIN, *args], capture_output=True, text=True).stderr

## ---------------------------------------------------------------- consistency

def read(path):
    return open(path, encoding="utf-8").read()

def plain(path):
    """The doc as a reader would copy it: HTML stripped, entities restored,
    line continuations joined, so one invocation is one line whatever the
    markup around it."""
    return re.sub(r"\\\n\s*", " ", strip(read(path)))

## A command invocation, not prose. Two rules, and both matter: the line must
## live in a CODE BLOCK (a ``` fence, a <pre>, or an indented llms.txt line),
## because "sesame is entirely offline" is a sentence, not a command; and the
## flags counted are only those before the first pipe or redirect, because
## `sesame attach-probe ... | tabl quantile --wider` ends in tabl's flags.
INVOKE = re.compile(r"^[ \t]*(?:\$ )?sesame ([a-z][a-z-]+)((?: [^\n]*)?)$", re.M)

def code_blocks(path):
    t = read(path)
    if path.endswith(".html"):
        return [strip(b) for b in re.findall(r"<pre[^>]*>(.*?)</pre>", t, re.S)]
    if path.endswith(".md"):
        return re.findall(r"^```[a-z]*\n(.*?)^```", t, re.S | re.M)
    ## llms.txt: commands are the indented lines
    return ["\n".join(l[2:] for l in t.splitlines() if l.startswith("  "))]

def strip(t):
    t = re.sub(r"<[^>]+>", "", t)
    for a, b in (("&amp;", "&"), ("&gt;", ">"), ("&lt;", "<"), ("&quot;", '"'),
                 ("&nbsp;", " "), ("&rarr;", "->"), ("&times;", "x")):
        t = t.replace(a, b)
    return t

def invocations(path):
    out = []
    for block in code_blocks(path):
        block = re.sub(r"\\\n\s*", " ", block)          # join continuations
        for m in INVOKE.finditer(block):
            rest = re.split(r"[|>]", m.group(2))[0]        # drop the next program
            out.append((m.group(1), rest))
    return out

def consistency():
    print("docs gate: consistency")

    ## 1. one version, in three places
    hdr = re.search(r'#define SESAME_VERSION "([^"]+)"', read(f"{ROOT}/include/sesame.h")).group(1)
    rec = re.search(r'set version = "([^"]+)"', read(f"{ROOT}/conda-recipe/meta.yaml")).group(1)
    binv = re.search(r"sesame ([0-9][^\s]*)", run("version"))
    check(hdr == rec, f"sesame.h ({hdr}) and meta.yaml ({rec}) agree on the version")
    check(bool(binv) and binv.group(1) == hdr,
          f"the binary reports {binv.group(1) if binv else '?'}, sesame.h says {hdr}")

    ## 2-4. every invocation printed in the docs names a real subcommand, and
    ## only flags that subcommand accepts. This is the check that would have
    ## caught `dml --platform` (dml has no such flag) before it shipped.
    help_txt = run("help")
    have = set(re.findall(r"^\s{2,}([a-z][a-z-]+)\s", help_txt, re.M))
    check("preprocess" in have, f"`sesame help` lists subcommands ({len(have)} found)")
    accepted = {c: set(re.findall(r"(--[a-z][a-z-]+)", run(c, "-h"))) for c in have}
    for d in DOCS:
        name = os.path.basename(d)
        before, inv = len(fails), invocations(d)
        for cmd, rest in inv:
            if cmd not in have:
                check(False, f"{name}: `sesame {cmd}` is not a subcommand"); continue
            for flag in re.findall(r"(--[a-z][a-z-]+)", rest):
                if flag not in accepted[cmd]:
                    check(False, f"{name}: `{cmd} {flag}` is not in `{cmd} -h`")
        if len(fails) == before:
            check(True, f"{name}: all {len(inv)} printed invocations match the binary")

    ## 5. the qc metric count the docs quote is the one the binary emits
    fields = len(re.findall(r"_\([ID], ",
                 re.search(r"#define SESAME_QC_FIELDS\(_\)(.*?)\n\n",
                           read(f"{ROOT}/include/sesame.h"), re.S).group(1)))
    for d in DOCS:
        for m in re.finditer(r"(\d\d)[- ]metric|holds (\d\d) (?:per-sample )?metrics", plain(d)):
            q = int(next(g for g in m.groups() if g))
            check(q == fields,
                  f"{os.path.basename(d)} quotes {q} qc metrics; the binary has {fields}")

    print(f"docs gate: {len(fails)} failure(s)")
    return 1 if fails else 0

## ------------------------------------------------------------------ workflows

def workflows():
    ex = sorted(glob.glob(os.path.join(ROOT, "docs", "examples", "*.sh")))
    if not ex:
        print("docs gate: no docs/examples/*.sh to run"); return 0
    store = os.environ.get("YAME_DATA_HOME") or os.path.expanduser("~/.local/share/yame")
    idats = os.environ.get("SESAME_TEST_IDATS") or os.path.expanduser("~/repo/InfiniumTestIDATs")
    if not os.path.isdir(store) or not os.path.isdir(idats):
        print(f"docs gate: SKIP (need a store at {store} and IDATs at {idats})"); return 0

    sandbox = os.environ.get("SESAME_DOCS_SANDBOX") or \
              os.path.join(os.environ.get("TMPDIR") or os.path.expanduser("~/tmp"), "sesame-docs")
    work = os.path.join(sandbox, "work")
    shutil.rmtree(work, ignore_errors=True); os.makedirs(work)
    env = dict(os.environ, PATH=ROOT + ":" + os.path.join(ROOT, "YAME") + ":" + os.environ["PATH"],
               YAME_DATA_HOME=store, SESAME_TEST_IDATS=idats)
    for path in ex:
        name = os.path.basename(path)
        norun = next((l[9:].strip() for l in open(path) if l.startswith("## norun:")), None)
        if norun:
            print(f"  SKIP {name}  ({norun})"); continue
        r = subprocess.run(["bash", "-eo", "pipefail", path], cwd=work, env=env,
                           capture_output=True, text=True)
        print(("  ok   " if r.returncode == 0 else "  FAIL ") + name)
        if r.returncode:
            fails.append(name)
            print("\n".join("        " + l for l in (r.stdout + r.stderr).splitlines()[-12:]))
    print(f"docs gate: {len(fails)} failure(s) over {len(ex)} example(s); sandbox {work}")
    return 1 if fails else 0

if __name__ == "__main__":
    sys.exit(consistency() if "--consistency" in sys.argv else workflows())
