# Contributing

### Bug reports, feature requests, and questions

To report a bug or request a feature, or if you have any questions,
either submit a [GitHub issue](https://github.com/rcsb/rcsb-chem-search/issues/new)
or contact [RCSB support](https://www.rcsb.org/pages/contactus) by
**emailing <a href="info@rcsb.org">info@rcsb.org</a>**.
Please use the email to report issues with the live service at
[rcsb.org](https://rcsb.org).

### Source code contributions

Contributions are always welcome.
Before writing any code, discuss the proposed changes with the developers by
[creating a GitHub issue](https://github.com/rcsb/rcsb-chem-search/issues/new)

Avoid making changes not related to the proposal so that your contributions are properly listed.
Feel free to open a draft pull request (PR) at any time – even long before the changes are complete.

### Recommended pull request workflow

We recommend the following steps, which use
[uv](https://docs.astral.sh/uv/)
and the
[GitHub CLI](https://cli.github.com/).

<b>Using a fork:</b>

```bash
#gh repo fork https://github.com/rcsb/rcsb-chem-search --default-branch-only --clone
gh repo fork rcsb-chem-search --clone --default-branch-only
cd rcsb-chem-search
uv run pre-commit install
```

<b>Using a branch (RCSB members):</b>

```bash
read -p "Branch name (e.g. fix/RO-1550--dmt or dev-dmt-): " branch
#gh repo fork https://github.com/rcsb/rcsb-chem-search --default-branch-only --clone
gh repo clone rcsb-chem-search -- --single-branch
git switch -c "$branch"
cd rcsb-chem-search
uv run pre-commit install
```

Create a PR by running

```bash
gh pr create --fill --web --draft
```

When you are ready, mark the PR as _ready for review_.
