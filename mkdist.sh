#!/bin/sh

cvs update -d

# hackery to build ChangeLog
darcs changes --summary > ChangeLog

sh autogen.sh

make maintainer-clean
./configure --enable-maintainer-mode
make dist
