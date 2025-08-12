# SPDX-FileCopyrightText: Copyright 2020-2025, Contributors to Tyrannosaurus
# SPDX-PackageHomePage: https://github.com/dmyersturnbull/tyrannosaurus
# SPDX-License-Identifier: Apache-2.0
#
# https://just.systems/man/en/
# https://cheatography.com/linux-china/cheat-sheets/justfile/
# CAUTION – Quotes aren't quotes.
# This doesn't do what you'd expect:
# ```
# git +args:
#   git {{args}}
# commit-reformat: (git "-m" "style: reformat code")
# ```
# This results in `git -m style: reformat code` !

set dotenv-load := true
set ignore-comments := true

# Confusingly, Python has a default warning filter, which is **different** from `-W default`.

export PYTHONWARNINGS := env('PYTHONWARNINGS', 'default')

# Enable dev mode (https://docs.python.org/3/library/devmode.html):

export PYTHONDEVMODE := env('PYTHONDEVMODE', '1')

# To complain when `encoding=` is omitted (in addition to Ruff rule):
# export PYTHONWARNDEFAULTENCODING := env('PYTHONWARNDEFAULTENCODING', '1')

git_sha := `git rev-parse --short=16 HEAD`
git_ref := `git rev-parse --abbrev-ref HEAD`
git_rev_date := `git --no-pager log -1 --format=%cd`

###################################################################################################

# List available recipes.
[group('help')]
list:
    @just --list --unsorted --list-prefix "  "

alias help := list

# Print info for a bug report.
[group('help')]
info:
    @echo "just_ver: {{ replace_regex(`just --version`, '^[^ ]+ ([^ ]+)$', '$1') }}"
    @echo "uv_ver: {{ replace_regex(`uv --version`, '^uv ([^ ]+) .+$', '$1') }}"
    @echo "git_ver: {{ replace_regex(`git --version`, '^git version ([^ ]+)$', '$1') }}"
    @echo "os: {{ os() }}"
    -@echo "os_ver: {{ `uname -r` }}"
    @echo "cpu_arch: {{ arch() }}"
    @echo "cpu_count: {{ num_cpus() }}"
    @echo "shell: {{ env('SHELL', '?') }}"
    @echo "lang: {{ env('LANG', '?') }}"
    @echo "repo_dir: {{ file_name(justfile_directory()) }}"
    @echo "project_name: {{ replace_regex(`uv version`, '^([^ ]+) .+$', '$1') }}"
    @echo "project_ver: {{ `uv version --short` }}"
    @echo "git_ref: {{ git_ref }}"
    @echo "git_sha: {{ git_sha }}"
    @echo "git_rev: {{ git_rev_date }}"

###################################################################################################

# Sync the venv and install commit hooks.
[group('setup')]
init:
    -rm .git/hooks/*.sample
    uv sync --all-extras --exact
    uv run --no-sync pre-commit install --install-hooks --overwrite
    uv run --no-sync pre-commit gc
    uv run --no-sync pre-commit run

###################################################################################################

# `just update fix-changes format-changes clean`
[group('project')]
spruce: update fix-changes format-changes clean

# `just update-lock update-hooks`
[group('project')]
update: update-lock update-hooks

alias upgrade := update

# Update the lock file and sync the venv.
[group('project')]
update-lock: sync
    uv run pre-commit gc

alias upgrade-lock := update-lock

# Auto-update commit hooks.
[group('project')]
update-hooks:
    uv run --locked pre-commit autoupdate
    uv run --no-sync pre-commit gc

alias upgrade-hooks := update-hooks

# Lock and sync the venv (uses all extras).
[group('project')]
sync:
    uv sync --all-extras --exact

alias lock := sync

# Prune temp data, including from uv, pre-commit, and git.
[group('project')]
clean: clean-temp clean-misc clean-git

# Delete temp/cache files and directories. Called by `clean`.
[group('project')]
clean-temp:
    # Generated docs:
    -rm -f -r .site/
    # Python temp files:
    -rm -f -r .ruff_cache/
    -rm -f -r .hypothesis/
    -rm -f -r **/__pycache__/
    -rm -f -r **/.pytest_cache/
    -rm -f -r **/cython_debug/
    -rm -f -r **/*.egg-info/
    -rm -f -r **/.node-modules/
    -rm -f .coverage.json
    -rm -f **/*.py[codi]
    # OS-specific files:
    -rm -f **/.directory
    -rm -f **/.DS_Store
    -rm -f **/.localized
    -rm -f **/Thumbs.db

# Delete unused uv and pre-commit cache data.
[group('project')]
clean-misc:
    uv run --no-sync pre-commit gc
    uv cache prune

# Run ultra-safe 'git maintenance' tasks. Called by `clean`.
[group('project')]
[no-exit-message]
clean-git:
    # The repo is small, so use `pack-refs` instead of `incremental-repack`.
    git maintenance run \
      --task=prefetch \
      --task=commit-graph \
      --task=loose-objects \
      --task=pack-refs

# Prune stale tracking refs, reflog items, etc. Pack refs. ❗
[group('project')]
[metadata('advanced')]
prune-git:
    # Needed on macOS.
    @-chflags -R nouchg .git/*
    git remote prune --all
    # Recover objects first to avoid deleting objects simply because they're unreachable.
    git fsck --lost-found --tags
    # `--task=gc` is NOT equivalent to `git gc`:
    # It works safely with other tasks and respects `gc.pruneExpire`.
    # I don't think `--task=pack-refs` is needed here, but it won't hurt.
    # Also, I've assumed the tasks are run in order. (If they're not, no harm done.)
    git maintenance run \
      --task=prefetch \
      --task=commit-graph \
      --task=reflog-expire \
      --task=rerere-gc \
      --task=loose-objects \
      --task=pack-refs \
      --task=gc

# Delete files that are probably temporary. 🧨
[group('project')]
[metadata('advanced')]
@delete-likely-temp:
    # Log files
    -rm -f **/*.log
    -rm -f **/*.log.gz
    -rm -f **/*.log.zst
    # Files and directories ending in `.(bak|junk|temp|tmp|trash)`
    -rm -f -r **/*.bak
    -rm -f -r **/*.junk
    -rm -f -r **/*.temp
    -rm -f -r **/*.tmp
    -rm -f -r **/*.trash
    # Files and directories starting or ending with `~|#|$` or starting with `.~`.
    -rm -f -r **/[~#\$]*
    -rm -f -r **/*[~#\$]
    -rm -f -r **/.~*

# Minify the repo by deleting nearly all recreatable files. 🧨
[group('project')]
[metadata('advanced')]
minify-repo: clean-temp clean-misc prune-git delete-likely-temp _dismantle

# <Internal.>
_dismantle:
    uv run pre-commit clean
    uv run pre-commit uninstall
    -rm -f -r .venv
    -rm -f -r .idea
    -rm -f uv.lock

###################################################################################################

# Format STAGED files (via pre-commit).
[group('format')]
format-changes: _format

alias format := format-changes

# Format ALL files (via pre-commit).
[group('format')]
format-all: (_format "--all-files")

# <Internal.>
_format *args:
    -uv run pre-commit run end-of-file-fixer {{ args }}
    -uv run pre-commit run fix-byte-order-marker {{ args }}
    -uv run pre-commit run trailing-whitespace {{ args }}
    -uv run pre-commit run ruff-format {{ args }}
    -uv run pre-commit run prettier {{ args }}

###################################################################################################

# Run pre-commit auto-fix hooks on STAGED files.
[group('fix')]
[no-exit-message]
fix-changes: _stage_precommit (hook "ruff-check")

alias fix := fix-changes

# Run pre-commit auto-fix hooks on ALL files. ⚠️
[group('fix')]
[no-exit-message]
fix-all: _stage_precommit (hook "ruff-check --all-files")

# Apply Ruff auto-fixes. (Runs on all files by default.) ⚠️
[group('fix')]
[no-exit-message]
fix-ruff *args: (_ruff_fix args)

# Preview Ruff auto-fixes. Runs on all files by default.
[group('fix')]
[no-exit-message]
diff-ruff *args: (_ruff_fix "--diff" args)

# <Internal.> (FYI: Use `--config` so users can still pass `--output-format`.)
[no-exit-message]
_ruff_fix *args: (ruff "check --fix-only --config 'output-format=\"grouped\"'" args)

# <Internal.>
[no-exit-message]
_stage_precommit:
    git add .pre-commit-config.yaml

###################################################################################################

# `just check-simple check-ruff check-ty check-links` (slow)
[group('check')]
check-all: check-core check-ruff check-ty check-schemas check-links

# Check basic rules (via pre-commit).
[group('check')]
check-core:
    # Keep slower hooks lower in the list.
    uv run --no-sync pre-commit run check-filenames
    uv run --no-sync pre-commit run pathvalidate
    uv run --no-sync pre-commit run check-symlinks
    uv run --no-sync pre-commit run check-case-conflict
    uv run --no-sync pre-commit run check-illegal-windows-names
    uv run --no-sync pre-commit run check-shebang-scripts-are-executable
    uv run --no-sync pre-commit run check-github-actions
    uv run --no-sync pre-commit run check-github-workflows
    uv run --no-sync pre-commit run check-compose-spec
    uv run --no-sync pre-commit run forbid-new-submodules
    uv run --no-sync pre-commit run trailing-whitespace
    uv run --no-sync pre-commit run check-added-large-files

# Check JSON schemas (via pre-commit).
[group('check')]
check-schemas:
    uv run --no-sync pre-commit run check-github-actions
    uv run --no-sync pre-commit run check-github-workflows
    uv run --no-sync pre-commit run check-compose-spec
    uv run --no-sync pre-commit run check-metaschema

# Check Ruff rules (without auto-fixing).
[group('check')]
[no-exit-message]
check-ruff *args: (_ruff_check args)

# Show the number of violations per Ruff rule.
[group('check')]
[no-exit-message]
check-ruff-stats *args: (ruff "check --no-fix --statistics" args)

# Check Ruff security (S) rules derived from Bandit.
[group('check')]
[no-exit-message]
check-bandit *args: (_ruff_check "--select S" args)

# Check Astral ty typing rules.
[group('check')]
[no-exit-message]
check-ty *args: (run "ty check" args)

# Detect broken hyperlinks via pre-commit. (slow)
[group('check')]
[no-exit-message]
check-links: (hook "markdown-link-check --hook-stage manual --all-files")

# <Internal.> (FYI: Use `--config` so users can still pass `--output-format`.)
[no-exit-message]
_ruff_check *args: (ruff "check --no-fix --config 'output-format=\"grouped\"'" args)

###################################################################################################

# Run tests not marked ux.
[group('test')]
[no-exit-message]
test-non-ux *args: (pytest "-m 'not ux'" args)

alias test := test-non-ux

# Run tests marked ux for manual interaction or verification.
[group('test')]
[no-exit-message]
test-ux *args: (pytest "-m ux" args)

# Run tests not marked ux, e2e, slow, or net.
[group('test')]
[no-exit-message]
test-fast *args: (pytest "-m 'not (ux or e2e or slow or net)'" args)

# Run doctest tests via PyTest.
[group('test')]
[no-exit-message]
test-doctest *args: (pytest "--doctest-modules src/" args)

alias doctest := test-doctest

# Run tests not marked ux, reporting coverage.
[group('test')]
[no-exit-message]
test-with-cov *args: (pytest "-m 'not ux' --cov" args)

alias test-cov := test-with-cov

# Run PyTest tests, highlighting test durations. ☆
[group('test')]
[no-exit-message]
test-durations *args: (pytest "--durations=0 --durations-min=0" args)

# Run tests marked hypothesis with explain phase enabled. ☆
[group('test')]
[no-exit-message]
test-hypothesis *args: (pytest "-m hypothesis --hypothesis-explain" args)

# Run PyTest tests stepwise, starting at last failure. ☆
[group('test')]
[no-exit-message]
test-stepwise *args: (pytest "--stepwise" args)

# Run PyTest tests, showing minimal output. ☆
[group('test')]
[no-exit-message]
test-quietly *args: (pytest "--capture=no --tb=line --quiet" args)

# Run PyTest tests, showing tracebacks, locals, and INFO. ☆
[group('test')]
[no-exit-message]
test-loudly *args: (pytest "--showlocals --full-trace --log-level=INFO --verbose" args)

# Run PyTest tests with pdb debugger. ☆
[group('test')]
[no-exit-message]
test-with-pdb *args: (pytest "--pdb" args)

###################################################################################################

# Build mkdocs docs from scratch, failing for warnings.
[group('docs')]
[no-exit-message]
build-docs *args: (run "mkdocs build --clean --strict" args)

# Locally serve the mkdocs docs.
[group('docs')]
[no-exit-message]
serve-docs *args: (run "mkdocs serve" args)

###################################################################################################

# `uv`
[group('alias')]
[no-exit-message]
[private]
@uv *args:
    uv {{ args }}

# `pre-commit`
[group('alias')]
[no-exit-message]
[private]
@pre-commit *args: (run "pre-commit" args)

# `uv run --locked`
[group('alias')]
[no-exit-message]
@run +args:
    uv run --locked {{ args }}

# `uv run --locked python`
[group('alias')]
[no-exit-message]
@python *args: (run "python" args)

# `uv run --locked pre-commit run {hook}`
[group('alias')]
[no-exit-message]
@hook name *args: (run "pre-commit run" name args)

# `uv run --locked ruff`
[group('alias')]
[no-exit-message]
@ruff *args: (run "ruff" args)

# `uv run --locked ty`
[group('alias')]
[no-exit-message]
@ty *args: (run "ty" args)

# `uv run --locked pytest`
[group('alias')]
[no-exit-message]
pytest *args: (run "pytest" args)

# `gh pr create --fill-verbose`
[group('alias')]
[no-exit-message]
pr *args="--web":
    gh pr create --fill-verbose {{ args }}
