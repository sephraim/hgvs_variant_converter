#!/bin/env python

import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper
import sys

# Initialize the mapper for GRCh37 with splign-based alignments
hdp = hgvs.dataproviders.uta.connect()
evm = hgvs.variantmapper.EasyVariantMapper(hdp)
hp = hgvs.parser.Parser()

# Read variants from file
variants = [line.rstrip('\n') for line in open(sys.argv[1])]

# Print file header
print "input\thgvs_g\thgvs_c\thgvs_p"

for variant in variants:
  var_c = ''
  var_g = ''
  var_p = ''

  try:
    if variant.startswith('NC_') | variant.startswith('NM_'):
      # Parse variant into Python structure
      if variant.startswith('NC_'):
        var_g = hp.parse_hgvs_variant(variant)
      else:
        # Convert C to G
        var_g = evm.c_to_g(hp.parse_hgvs_variant(variant))
      
      # Identify transcripts that overlap this genomic variant
      transcripts = evm.relevant_transcripts(var_g) #=> ['NM_001637.3', 'NM_001177506.1', 'NM_001177507.1']
      
      if not transcripts:
        # No transcripts exist
        print "%s\t\t\t" % variant
      else:
        # Map genomic variant to one of these transcripts
        for t in transcripts:
          var_c = evm.g_to_c(var_g, t)
          var_p = evm.c_to_p(var_c)
          print "%s\t%s\t%s\t%s" %(variant, var_g, str(var_c), str(var_p))
    else:
      # No variant provided
      print "%s\t\t\t" % variant
  except:
    # Something went wrong with this variant
    print "%s\t\t\t" % variant
