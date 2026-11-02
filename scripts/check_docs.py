#!/usr/bin/env python3
"""Check the project's Markdown math and report templates before publication.

This intentionally permits a small set of basic TeX commands, rather than every
command accepted by a local renderer. GitHub's Markdown math syntax is described
at https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/writing-mathematical-expressions
"""
import argparse
from pathlib import Path
import re
import subprocess
import sys

SAFE_COMMANDS = set("""
Lambda alpha begin beta cos end epsilon frac hat infty int lVert lambda ldots
left leq mathrm max qquad rVert right sigma sim sin sqrt sum tau theta times
top varepsilon
""".split())
SAFE_ENVIRONMENTS = {"aligned", "bmatrix", "gathered"}
COMMAND = re.compile(r"(?<!\\)\\([A-Za-z]+)")
R_COMMAND = re.compile(r"(?<!\\)\\\\([A-Za-z]+)")
FENCE = re.compile(r"^ {0,3}(`{3,}|~{3,})(.*)$")


def check_file(path):
    source = path.read_text(encoding="utf-8")
    errors = []
    # R needs two backslashes in a string to emit one TeX command.
    pattern = R_COMMAND if path.suffix.lower() == ".r" else COMMAND
    for match in pattern.finditer(source):
        if match[1] not in SAFE_COMMANDS:
            line = source.count("\n", 0, match.start()) + 1
            errors.append(f"{path}:{line}: unsupported math command: {match[1]}")
    if path.suffix.lower() != ".md":
        return errors

    fence = None
    math_lines = []
    for number, line in enumerate(source.splitlines(), 1):
        marker = FENCE.match(line)
        if marker and fence is None:
            fence = (marker[1][0], len(marker[1]), marker[2].strip(), number)
            continue
        if fence is not None:
            if marker and marker[1][0] == fence[0] and len(marker[1]) >= fence[1] and not marker[2].strip():
                fence = None
            elif fence[2] == "math":
                if "$" in line:
                    errors.append(f"{path}:{number}: use a math fence without dollar delimiters")
                math_lines.append((number, line))
            continue
        # Inline code may contain literal dollar signs.
        prose = re.sub(r"(`+).*?\1", "", line)
        dollars = len(re.findall(r"(?<!\\)\$", prose))
        if dollars % 2:
            errors.append(f"{path}:{number}: unmatched dollar delimiter")
        if dollars:
            math_lines.append((number, prose))
    if fence is not None:
        errors.append(f"{path}:{fence[3]}: unclosed code fence")

    environments = []
    braces = []
    for number, line in math_lines:
        for match in re.finditer(r"\\(begin|end)\{([^}]+)\}", line):
            action, name = match.groups()
            if name not in SAFE_ENVIRONMENTS:
                errors.append(f"{path}:{number}: unsupported math environment: {name}")
            if action == "begin":
                environments.append((name, number))
            elif not environments or environments[-1][0] != name:
                errors.append(f"{path}:{number}: mismatched math environment: {name}")
            else:
                environments.pop()
        for match in re.finditer(r"(?<!\\)[{}]", line):
            if match[0] == "{":
                braces.append(number)
            elif braces:
                braces.pop()
            else:
                errors.append(f"{path}:{number}: unmatched closing brace")
    for name, number in environments:
        errors.append(f"{path}:{number}: unclosed math environment: {name}")
    for number in braces:
        errors.append(f"{path}:{number}: unclosed math brace")
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", type=Path)
    args = parser.parse_args()
    if args.paths:
        paths = args.paths
    else:
        root = Path(__file__).resolve().parents[1]
        tracked = subprocess.check_output(["git", "-C", str(root), "ls-files", "-z"], text=True).split("\0")
        paths = [root / p for p in tracked if p and
                 (Path(p).suffix.lower() == ".md" or Path(p).name == "report.R")]
    if not paths:
        parser.error("No documentation files found")
    errors = [error for path in paths for error in check_file(path)]
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"Documentation math checks passed for {len(paths)} files.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
