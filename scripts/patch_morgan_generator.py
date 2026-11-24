#!/usr/bin/env python3

import ast
from pathlib import Path

files = [
    Path("scripts/scaffold_libinvent_v5_devicefix.py"),
    Path("scripts/cumulative_score.py"),
]

for path in files:
    source = path.read_text()
    tree = ast.parse(source)
    lines = source.splitlines(keepends=True)
    offsets = [0]

    for line in lines:
        offsets.append(offsets[-1] + len(line))

    replacements = []

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue

        if not isinstance(node.func, ast.Attribute):
            continue

        if node.func.attr != "GetMorganFingerprintAsBitVect":
            continue

        args = [ast.unparse(a) for a in node.args]
        kwargs = {
            k.arg: ast.unparse(k.value)
            for k in node.keywords
        }

        supported = {
            "radius", "nBits", "useChirality",
            "useBondTypes", "mol"
        }

        if len(args) > 3 or set(kwargs) - supported:
            print(f"SKIPPED: {path}:{node.lineno}")
            continue

        mol = args[0] if args else kwargs.get("mol")
        radius = args[1] if len(args) > 1 else kwargs.get("radius")
        nbits = args[2] if len(args) > 2 else kwargs.get("nBits", "2048")

        if mol is None or radius is None:
            print(f"SKIPPED: {path}:{node.lineno}")
            continue

        chirality = kwargs.get("useChirality", "False")
        bondtypes = kwargs.get("useBondTypes", "True")

        replacement = (
            "rdFingerprintGenerator.GetMorganGenerator("
            f"radius={radius}, fpSize={nbits}, "
            f"includeChirality={chirality}, "
            f"useBondTypes={bondtypes}"
            f").GetFingerprint({mol})"
        )

        start = offsets[node.lineno - 1] + node.col_offset
        end = offsets[node.end_lineno - 1] + node.end_col_offset

        replacements.append((start, end, replacement))

    if not replacements:
        print(f"{path}: no supported deprecated calls found")
        continue

    for start, end, replacement in sorted(
        replacements, reverse=True
    ):
        source = source[:start] + replacement + source[end:]

    if "from rdkit.Chem import rdFingerprintGenerator" not in source:
        updated_tree = ast.parse(source)

        imports = [
            node for node in updated_tree.body
            if isinstance(node, (ast.Import, ast.ImportFrom))
        ]

        if not imports:
            raise RuntimeError(f"No import section found in {path}")

        source_lines = source.splitlines(keepends=True)
        insert_at = imports[-1].end_lineno

        source_lines.insert(
            insert_at,
            "from rdkit.Chem import rdFingerprintGenerator\n"
        )

        source = "".join(source_lines)

    ast.parse(source)
    path.write_text(source)

    print(f"PATCHED: {path}")
    print(f"Replaced {len(replacements)} deprecated calls")

print("MorganGenerator patch completed.")
