# SPDX-FileCopyrightText: Copyright 2020-2024, Contributors to CICD
# SPDX-PackageHomePage: https://github.com/dmyersturnbull/cicd
# SPDX-License-Identifier: Apache-2.0
#
# SPDX-FileCopyrightText: Copyright 2024, the RCSB PDB and contributors
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem-search
# SPDX-License-Identifier: BSD-3-Clause
#
# Original source code is from:
# CICD (https://github.com/dmyersturnbull/cicd).
# This file includes modifications.
#
"""
Root of rcsb-chem-search.
Imports public modules for easy usage:

```python
from rcsbchemsearch import ProjectMetadata #, ...
```
"""

from rcsbchemsearch._project_metadata import ProjectMetadata as __ProjectMetadata

__metadata__ = __ProjectMetadata
__uri__ = __ProjectMetadata.homepage
__title__ = __ProjectMetadata.title
__summary__ = __ProjectMetadata.summary
__version__ = __ProjectMetadata.version
__license__ = __ProjectMetadata.license
