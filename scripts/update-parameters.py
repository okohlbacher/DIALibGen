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
    # -write_ini leaves the deprecated names out (an INI must load once they
    # are removed); --helphelp still lists them.
    helphelp = subprocess.run([str(args.binary.resolve()), "--helphelp"], capture_output=True, text=True)
    deprecated = []
    for line in (helphelp.stdout + helphelp.stderr).splitlines():
        marker = "Deprecated (0.12, removed in 0.13): use "
        if line.strip().startswith("-") and marker in line:
            replacement = line[line.index(marker) + len(marker):].split(" (default:")[0]
            deprecated.append("| `" + line.split()[0] + "` | " + replacement.replace("|", "\\|") + " |")

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
            shown = "" if item.get("advanced") == "true" else "yes"
            yield "| " + " | ".join(map(cell, ("`-" + name + "`", kind,
                                                "`" + default + "`" if default else "(empty)",
                                                constraints, shown, item.get("description", "")))) + " |"


tool = root.find("./NODE")
instance = tool.find("./NODE") if tool is not None else None
if instance is None:
    raise SystemExit("TOPP schema contains no tool instance")
version = tool.find("./ITEM[@name='version']")
label = version.get("value", "unknown") if version is not None else "unknown"
text = ("# Parameter reference\n\nGenerated from DIALibGen " + label +
        " with `scripts/update-parameters.py`. Do not edit this table by hand.\n\n"
        "The schema lists all modes. Without `-mode`, a FASTA generates and a library converts; "
        "refinement and training options apply to their respective modes. `--help` shows the "
        "options marked in the column *--help* (`--help -mode <mode>` one mode's), `--helphelp` "
        "all of them. Generation booleans take `true`/`false`; list items are separate arguments. "
        "`-threads` defaults to all physical cores unless the command line or the INI sets it. "
        "See [usage](usage.md).\n\n"
        "| Option | Type | Default | Constraints | --help | Description |\n"
        "|---|---|---|---|---|---|\n" + "\n".join(rows(instance)) + "\n" +
        ("\n## Deprecated options\n\nAccepted with a WARNING in 0.12 and removed in 0.13; `-write_ini` "
         "does not write them.\n\n| Option | Replacement |\n|---|---|\n" + "\n".join(deprecated) + "\n"
         if deprecated else ""))
target = Path(__file__).resolve().parents[1] / "docs" / "parameters.md"
if args.check:
    if not target.exists() or target.read_text() != text:
        raise SystemExit("Parameter reference is stale; run scripts/update-parameters.py BINARY")
else:
    target.write_text(text)
print("Parameter reference " + ("matches the executable" if args.check else "updated"))
