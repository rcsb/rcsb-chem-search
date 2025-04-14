# SPDX-FileCopyrightText: Copyright 2024, Contributors to rcsb-chem-search
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem-search
# SPDX-License-Identifier: BSD-3-Clause
"""
Root of rcsb-chem.

Imports public modules for easy usage:

```python
from rcsbchemsearch import __about__, Tautomer
```
"""

from rcsbchemsearch.core.about import __about__

__uri__ = __about__["urls"]["homepage"]
__title__ = __about__["name"]
__summary__ = __about__["summary"]
__version__ = __about__["version"]
__license__ = __about__["license"]
