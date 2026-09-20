#!/usr/bin/env python3
"""Generate or verify the parameter reference from DIALibGen's TOPP schema."""
import argparse
from pathlib import Path
import subprocess
import tempfile
import xml.etree.ElementTree as ET

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("binary", type=Path)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
with tempfile.TemporaryDirectory(prefix="dialibgen-parameters-") as directory:
    ini = Path(directory) / "parameters.ini"
    result = subprocess.run([str(args.binary.resolve()), "-write_ini", str(ini)],
                            capture_output=True, text=True)
    if result.returncode:
        raise SystemExit(result.stderr.strip() or f"Schema export failed: {result.returncode}")
    root = ET.parse(ini).getroot()


def cell(value):
    return str(value).replace("|", "\\|").replace("\n", " ")


def rows(node, prefix=""):
    for item in node:
        name = prefix + item.get("name", "")
        if item.tag == "NODE":
            yield from rows(item, name + ":")
        elif item.tag in ("ITEM", "ITEMLIST"):
            default = item.get("value", "")
            if item.tag == "ITEMLIST":
                default = ", ".join(child.get("value", "") for child in item)
            constraints = item.get("restrictions", "")
            kind = item.get("type", "") + (" list" if item.tag == "ITEMLIST" else "")
            yield "| " + " | ".join(map(cell, ("`-" + name + "`", kind,
                                                "`" + default + "`" if default else "(empty)",
                                                constraints, item.get("description", "")))) + " |"


tool = root.find("./NODE")
instance = tool.find("./NODE") if tool is not None else None
if instance is None:
    raise SystemExit("TOPP schema contains no tool instance")
version = tool.find("./ITEM[@name='version']")
label = version.get("value", "unknown") if version is not None else "unknown"
text = ("# Parameter reference\n\nGenerated from DIALibGen " + label +
        " with `scripts/update-parameters.py`. Do not edit this table by hand.\n\n"
        "The schema lists all modes. `-mode generate` is the default; refinement and "
        "training options apply to their respective modes. Generation booleans take "
        "`true`/`false`; list items are separate arguments. See [usage](usage.md).\n\n"
        "| Option | Type | Default | Constraints | Description |\n"
        "|---|---|---|---|---|\n" + "\n".join(rows(instance)) + "\n")
target = Path(__file__).resolve().parents[1] / "docs" / "parameters.md"
if args.check:
    if not target.exists() or target.read_text() != text:
        raise SystemExit("Parameter reference is stale; run scripts/update-parameters.py BINARY")
else:
    target.write_text(text)
print("Parameter reference " + ("matches the executable" if args.check else "updated"))
