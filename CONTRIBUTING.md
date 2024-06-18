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

We recommend the following steps, which use the
[GitHub CLI](https://cli.github.com/).

```bash
python -m install --upgrade pip
pip install hatch pre-commit
gh repo fork https://github.com/rcsb/rcsb-chem-search --default-branch-only --clone
cd rcsb-chem-search
pre-commit install
```

Create a PR by running

```bash
gh pr create --fill --web --draft
```

When you are ready, mark the PR as *ready for review*.
