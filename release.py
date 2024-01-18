#!/usr/bin/env python3

import io
from os.path import dirname, join
import extern

def get_version(relpath):
  """Read version info from a file without importing it"""
  for line in io.open(join(dirname(__file__), relpath), encoding="utf-8"):
    if "version" in line:
      if '"' in line:
        return line.split('"')[1]
      elif "'" in line:
        return line.split("'")[1]

if __name__ == "__main__":
    print("running release.sh")
    extern.run("./release.sh")

    version = get_version('Cargo.toml')
    print("version is {}".format(version))

    # Replace version in CITATION.cff
    citations_lines = []
    with open("CITATION.cff", "r") as f:
        import re
        r = re.compile(r"( *version: )")
        r2 = re.compile(r"( *date-released: )")
        for line in f:
            if matches := r.match(line):
                line = matches.group(1) + version + "\n"
            elif matches := r2.match(line):
                from datetime import datetime
                line = matches.group(1) + datetime.today().strftime('%Y-%m-%d') + "\n"
            citations_lines.append(line)
    with open("CITATION.cff", "w") as f:
        f.writelines(citations_lines)

    print("Checking if repo is clean ..")
    extern.run('if [[ $(git diff --shortstat 2> /dev/null | tail -n1) != "" ]]; then exit 1; fi')

    extern.run('git tag v{}'.format(version))
    print("Now run 'cargo publish' and then 'git push && git push --tags'".format(version))
    print("Then make a release, adding changelog info, so Zenodo picks it up")
