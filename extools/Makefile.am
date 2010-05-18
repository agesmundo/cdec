bin_PROGRAMS = \
  extractor \
  mr_stripe_rule_reduce

noinst_PROGRAMS =

mr_stripe_rule_reduce_SOURCES = mr_stripe_rule_reduce.cc extract.cc sentence_pair.cc
mr_stripe_rule_reduce_LDADD = $(top_srcdir)/decoder/libcdec.a -lz
mr_stripe_rule_reduce_LDFLAGS = -static

extractor_SOURCES = sentence_pair.cc extract.cc extractor.cc
extractor_LDADD = $(top_srcdir)/decoder/libcdec.a -lz
extractor_LDFLAGS = -static

AM_CPPFLAGS = -W -Wall -Wno-sign-compare $(GTEST_CPPFLAGS) -I$(top_srcdir)/decoder

