add_library(fasta_validator
        fasta_validator.c
        fasta_validator.h
)

set_target_properties(fasta_validator PROPERTIES COMPILE_FLAGS "${MMSEQS_CXX_FLAGS} -w" LINK_FLAGS "${MMSEQS_CXX_FLAGS} -w")

