#!/usr/bin/env python3
# Copyright © 2026 GSI Helmholtzzentrum für Schwerionenforschung gmbh
# <https://gsi.de/>
#
# Everyone is permitted to copy and distribute verbatim copies of this license
# document, but changing it is not allowed.
#
# Preamble
# The GNU Affero General Public License (version 3) is a free, copyleft
# license for software and other kinds of works. The licenses for most
# software and other practical works are designed to take away your freedom
# to share and change the works. By contrast, the GNU Affero General Public
# License is intended to guarantee your freedom to share and change all
# versions of a program--to make sure it remains free software for all its
# users. We, GSI Helmholtzzentrum für Schwerionenforschung, use the GNU
# Affero General Public License for most of our software; it applies also to
# any other work released this way by its authors. You can apply it to your
# programs, too.

"""Doxygen INPUT_FILTER for markdown files.

Doxygen does not understand GitHub-style $...$ math: it requires its own
\\f$...\\f$ (inline) and \\f[...\\f] (block) delimiters. This filter rewrites
the markdown math on the fly so the source files can keep the GitHub syntax.
It also rewrites markdown links to local repo files/folders, which doxygen
cannot resolve (they become unresolvable \\ref commands), into links to the
generated docs page or to the GitLab web UI of this repository.

Usage (called by doxygen, one file path per invocation):
    python3 docs/md_math_filter.py <input-file>
The transformed content is written to stdout.
"""
import hashlib
import os
import re
import subprocess
import sys

path = sys.argv[1]
with open(path, encoding="utf-8") as f:
    src = f.read()


def dir_page_html(name):
    """HTML page name of doxygen's directory-reference page for a
    top-level directory (doxygen names it dir_<md5('name/')>.html)."""
    return "dir_" + hashlib.md5(name.encode()).hexdigest() + ".html"


def _repo_url():
    """Base URL of the GitLab project (no trailing .git/)."""
    base = os.environ.get("CI_PROJECT_URL", "").strip()
    if not base:
        try:
            out = subprocess.check_output(
                ["git", "remote", "-v"], text=True, stderr=subprocess.DEVNULL
            )
            # each line is: <name> <url> (fetch|push)
            url = out.splitlines()[0].split()[1]
        except Exception:
            return ""
        if url.startswith("git@") and ":" in url:
            host, proj = url[4:].split(":", 1)
            base = "https://%s/%s" % (host, proj)
        else:
            base = url
    return base.rstrip("/").removesuffix(".git")


def _repo_branch():
    branch = os.environ.get("CI_COMMIT_REF_NAME", "").strip()
    if not branch:
        try:
            branch = subprocess.check_output(
                ["git", "rev-parse", "--abbrev-ref", "HEAD"],
                text=True,
                stderr=subprocess.DEVNULL,
            ).strip()
        except Exception:
            branch = "main"
    return branch


# repo-relative path (as written in docs/readme.md) -> (url kind, repo path)
REPO_LINKS = {
    "validation/": ("tree", "docs/validation"),
    "../macros/generated": ("tree", "macros/generated"),
    "primaries.csv": ("raw", "docs/primaries.csv"),
}


def _gitlab_link(m):
    target = m.group(2).strip()
    kind, repo_path = REPO_LINKS.get(target, ("", ""))
    if not kind:
        return m.group(0)
    base = _repo_url()
    if not base:
        return m.group(0)
    url = "%s/-/%s/%s/%s" % (base, kind, _repo_branch(), repo_path)
    return '<a href="%s">%s</a>' % (url, m.group(1))


if path.endswith(".md"):
    # block math first (may span lines), then inline math (single line)
    src = re.sub(r"\$\$(.*?)\$\$", r"\\f[\1\\f]", src, flags=re.S)
    src = re.sub(r"(?<!\$)\$([^$\n]+?)\$", r"\\f$\1\\f$", src)
    # markdown links to local files in the docs/ folder cannot be resolved
    # by doxygen (they become unresolvable \ref commands): point them to
    # the generated "docs Directory Reference" page instead
    src = re.sub(
        r"\[([^\]]+)\]\((?:\./)?docs/[^)\s]*\)",
        lambda m: '<a href="%s">%s</a>' % (dir_page_html("docs/"), m.group(1)),
        src,
    )
    # links to other local repo files/folders: point at the GitLab web UI
    src = re.sub(
        r"(?<!!)\[([^\]]+)\]\((validation/|\.\./macros/generated|primaries\.csv)\)",
        _gitlab_link,
        src,
    )

sys.stdout.write(src)
