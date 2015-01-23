# A simple GIMP plug-in Makefile by Yeti <yeti@physics.muni.cz>.
# This file is in public domain.
# You may need to change EXTRA_CFLAGS, CFLAGS and/or LDFLAGS, namely when your
# compiler is not gcc.
PACKAGE = voronoi
VERSION = 2.2
GIMPTOOL = gimptool-2.0

GIMP_LDFLAGS = `$(GIMPTOOL) --libs`
GIMP_CFLAGS = `$(GIMPTOOL) --cflags`
EXTRA_CFLAGS = -DLOCALEDIR=\"`$(GIMPTOOL) --prefix`/share/locale\" -DVERSION=\"$(VERSION)\" -Wall -O2 -march=`uname -m`

CFLAGS = $(GIMP_CFLAGS) $(EXTRA_CFLAGS)
LDFLAGS = $(GIMP_LDFLAGS) -lm

DNAME = $(PACKAGE)-$(VERSION)
SOURCES = voronoi.c voronoi2d.c randgrid.c voronoi.h
EXTRA_DIST = NEWS

all: $(PACKAGE)

install: $(PACKAGE)
	$(GIMPTOOL) --install-admin-bin $(PACKAGE)

uninstall: $(PACKAGE)
	$(GIMPTOOL) --uninstall-admin-bin $(PACKAGE)

install-user: $(PACKAGE)
	$(GIMPTOOL) --install-bin $(PACKAGE)

uninstall-user: $(PACKAGE)
	$(GIMPTOOL) --uninstall-bin $(PACKAGE)

clean:
	-rm -f *~ *.o core $(PACKAGE)

dist:
	-rm -rf $(DNAME)
	mkdir $(DNAME) $(DNAME)/libgimp
	chmod a+rx $(DNAME) $(DNAME)/libgimp
	cp $(SOURCES) $(EXTRA_DIST) Makefile README COPYING $(DNAME)
	cp libgimp/stdplugins-intl.h $(DNAME)/libgimp
	chmod a+r $(DNAME)/* $(DNAME)/libgimp/*
	tar cf - $(DNAME) | bzip2 > $(DNAME).tar.bz2
	rm -rf $(DNAME)

$(PACKAGE): voronoi.o voronoi2d.o randgrid.o

voronoi.o: voronoi.c voronoi.h

voronoi2d.o: voronoi2d.c voronoi.h

randgrid.o: randgrid.c voronoi.h

.PHONY: all clean install install-user uninstall uninstall-user dist

