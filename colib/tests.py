# coding: utf-8
"""

Testing plan:

1. load genbank file
2. create strain
    2.1. apply breseq/variant call mutations
    2.2. apply feature removal, addition mutations
3. load SBOL part(s)
    3.1. apply part insertion mutations
        3.1.1. apply those relative to features
4. lineage
5. store/load from database
6. store/load components & mutations from JSON (for REST api)

"""
from __future__ import unicode_literals


