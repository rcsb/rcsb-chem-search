# SPDX-FileCopyrightText: Copyright 2024, Contributors to rcsb-chem-search
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem-search
# SPDX-License-Identifier: BSD-3-Clause
"""
Root of rcsb-chem-search.
Imports public modules for easy usage:

```python
from rcsbchemsearch import ProjectMetadata #, ...
```
"""

from rcsbchemsearch._about import about as __about__

__uri__ = __about__.homepage
__title__ = __about__.title
__summary__ = __about__.summary
__version__ = __about__.version
__license__ = __about__.license
