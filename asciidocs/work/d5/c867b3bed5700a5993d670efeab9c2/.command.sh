#!/bin/bash -ue
cat chunk_aa chunk_ab | tr '[a-z]' '[A-Z]' > outfile
