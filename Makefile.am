bin_PROGRAMS = harminv
lib_LTLIBRARIES = libharminv.la
include_HEADERS = harminv.h
dist_man_MANS = harminv.1
noinst_PROGRAMS = sines

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = harminv.pc

ACLOCAL_AMFLAGS=-I m4

EXTRA_DIST = COPYRIGHT harminv.pc.in m4 sines-test.sh

libharminv_la_SOURCES = harminv.c harminv.h harminv-int.h check.h
libharminv_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

harminv_SOURCES = harminv-main.c copyright.h
harminv_LDADD = libharminv.la

sines_SOURCES = sines.c

BUILT_SOURCES = copyright.h

TESTS = sines-test.sh

# Generate built sources only in maintainer mode
if MAINTAINER_MODE

copyright.h: COPYRIGHT
	(echo '#define COPYRIGHT \'; sed 's,^/\* *,",;s,^ \* *,",;s,$$,\\n" \\,' COPYRIGHT  |egrep -v '^"/|^[^"]'; echo) > $@

maintainer-clean-local:
	rm -f $(BUILT_SOURCES)

endif # MAINTAINER_MODE
