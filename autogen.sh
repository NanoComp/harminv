#!/bin/sh

autoreconf --verbose --install --symlink --force
autoreconf --verbose --install --symlink --force
autoreconf --verbose --install --symlink --force

./configure --enable-maintainer-mode "$@"
