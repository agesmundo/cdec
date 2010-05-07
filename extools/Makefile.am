bin_PROGRAMS = \
  extractor

noinst_PROGRAMS =

extractor_SOURCES = extractor.cc sentence_pair.cc
extractor_LDADD = $(top_srcdir)/decoder/libcdec.a -lz

AM_CPPFLAGS = -W -Wall -Wno-sign-compare $(GTEST_CPPFLAGS) -I$(top_srcdir)/decoder

