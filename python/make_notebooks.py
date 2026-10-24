"""Build the Jupyter notebooks from the step-by-step scripts.

    python make_notebooks.py            # build ex1_cT.ipynb and ex2_cT_r.ipynb
    python make_notebooks.py --execute  # ... and run them, saving the plots

The .py scripts are the source; edit them, not the notebooks, then re-run this.
"""
import os
import re
import sys
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HERE = os.path.dirname(os.path.abspath(__file__))
BAR = re.compile(r"^# =+\s*$")
TITLES = {"ex1_cT.py": "Example 1: estimating the killing rate $c_T$",
          "ex2_cT_r.py": "Example 2: estimating $c_T$ and the activation rate $r$"}


def split_steps(src):
    """Return (docstring, [(header_lines, code), ...])."""
    m = re.match(r'\s*"""(.*?)"""\s*\n', src, flags=re.S)
    doc, rest = m.group(1), src[m.end():]
    lines = rest.splitlines()
    blocks, i, pre = [], 0, []
    while i < len(lines) and not BAR.match(lines[i]):
        pre.append(lines[i]); i += 1
    while i < len(lines):
        i += 1                                   # opening bar
        header = []
        while i < len(lines) and not BAR.match(lines[i]):
            header.append(re.sub(r"^#\s?", "", lines[i])); i += 1
        i += 1                                   # closing bar
        code = []
        while i < len(lines) and not BAR.match(lines[i]):
            code.append(lines[i]); i += 1
        blocks.append((header, "\n".join(code).strip("\n")))
    return doc, "\n".join(pre).strip("\n"), blocks


def header_to_md(header):
    title, body = header[0].strip(), [h.rstrip() for h in header[1:]]
    para = []
    # indented continuation lines right after the title
    while body and body[0].strip() and body[0].startswith("    "):
        line = body.pop(0).strip()
        if title.endswith((".", ":", ")")):
            para.append(line)
        else:
            title += " " + line
    while body and not body[0].strip():
        body.pop(0)
    md = "### " + title.replace("STEP", "Step")
    if para:
        md += "\n\n" + " ".join(para)
    if any(b.strip() for b in body):
        md += "\n\n```text\n" + "\n".join(body).strip("\n") + "\n```"
    return md


def build(script):
    src = open(os.path.join(HERE, script)).read()
    doc, pre, blocks = split_steps(src)
    doc_lines = [l for l in doc.strip("\n").splitlines() if not set(l.strip()) <= {"="}]
    title = TITLES[script]
    nb = new_notebook()
    nb.cells.append(new_markdown_cell(
        f"# {title}\n\n"
        "Notebook version of `" + script + "` (same code, split by steps).\n\n"
        "```text\n" + "\n".join(doc_lines).strip("\n") + "\n```"))
    nb.cells.append(new_code_cell(pre + "\n%matplotlib inline"))
    for header, code in blocks:
        nb.cells.append(new_markdown_cell(header_to_md(header)))
        if code:
            nb.cells.append(new_code_cell(code))
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    out = os.path.join(HERE, script.replace(".py", ".ipynb"))
    nbformat.write(nb, out)
    return out


if __name__ == "__main__":
    for s in ["ex1_cT.py", "ex2_cT_r.py"]:
        path = build(s)
        if "--execute" in sys.argv:
            from nbconvert.preprocessors import ExecutePreprocessor
            nb = nbformat.read(path, as_version=4)
            ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb, {"metadata": {"path": HERE}})
            nbformat.write(nb, path)
        print("wrote", path)
