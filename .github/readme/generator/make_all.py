"""Regenerate every README figure:  python .github/readme/generator/make_all.py"""
import os
import runpy

HERE = os.path.dirname(os.path.abspath(__file__))
for name in ["banner.py", "roadmap.py", "gif_observer.py", "highlights.py"]:
    runpy.run_path(os.path.join(HERE, name), run_name="__main__")
