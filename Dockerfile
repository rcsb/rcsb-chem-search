# SPDX-FileCopyrightText: Copyright 2020-2024, Contributors to Tyrannosaurus
# SPDX-PackageHomePage: https://github.com/dmyersturnbull/tyrannosaurus
# SPDX-License-Identifier: Apache-2.0
#
# SPDX-FileCopyrightText: Copyright 2024, Contributors to rcsb-chem-search
# SPDX-PackageHomePage: https://github.com/rcsb/rcsb-chem-search
# SPDX-License-Identifier: BSD-3-Clause
#
# Adapted from Tyrannosaurus <https://github.com/dmyersturnbull/tyrannosaurus>.
# Also includes elements from:
# https://stackoverflow.com/questions/53835198/integrating-python-poetry-with-docker/54763270#54763270

# :tyranno: FROM python:${.python-version-in-cicd.semver_minor(@)~}
FROM python:3.12

# --------------------------------------
# ------------- Set labels -------------

# See https://github.com/opencontainers/image-spec/blob/master/annotations.md
# :tyranno: LABEL org.opencontainers.image.version="${project.version}"
LABEL org.opencontainers.image.version="0.0.1-alpha0"
# :tyranno: LABEL org.opencontainers.image.vendor="${tool.cicd.data.vendor}"
LABEL org.opencontainers.image.vendor="rcsb"
# :tyranno: LABEL org.opencontainers.image.title="${project.name}"
LABEL org.opencontainers.image.title="rcsbchemsearch"
# :tyranno: LABEL org.opencontainers.image.url="${project.urls.Homepage}"
LABEL org.opencontainers.image.url="https://github.com/rcsb/rcsbchemsearch"
# :tyranno: LABEL org.opencontainers.image.documentation="${project.urls.docs}"
LABEL org.opencontainers.image.documentation="https://github.com/rcsb/py-rcsbchemsearch"

# --------------------------------------
# ---------- Copy and install ----------

# ENV no longer adds a layer in new Docker versions,
# so we don't need to chain these in a single line
ENV PYTHONFAULTHANDLER=1
ENV PYTHONUNBUFFERED=1
ENV PYTHONHASHSEED=random
ENV PIP_DISABLE_PIP_VERSION_CHECK=on
ENV PIP_DEFAULT_TIMEOUT=120

WORKDIR /app
COPY . /app/

RUN pip install .

EXPOSE 80
EXPOSE 443

CMD hypercorn\
  rcsbchemsearch.api:app\
 --bind '[::]:80'\
 --bind '[::]:443'\
 --quic-bind '[::]:443'
