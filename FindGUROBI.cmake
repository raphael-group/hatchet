FIND_PATH(GUROBI_INCLUDE_DIR
          NAMES "gurobi_c++.h" "gurobi_c.h"
          PATHS /n/fs/ragr-code/general/gurobi751/linux64/include/ /Library/gurobi702/mac64/include/
          DOC "Gurobi include directory")

FIND_LIBRARY(GUROBI_CPP_LIB
             NAMES gurobi_c++ 
             PATHS /n/fs/ragr-code/general/gurobi751/linux64/lib/ /Library/gurobi702/mac64/lib/
             DOC "Gurobi C++ Libraries")

FIND_LIBRARY(GUROBI_LIB
             NAMES gurobi70
             PATHS /n/fs/ragr-code/general/gurobi751/linux64/lib/ /Library/gurobi702/mac64/lib/
             DOC "Gurobi C Libraries")

set(GUROBI_LIBRARIES ${GUROBI_CPP_LIB} ${GUROBI_LIB})

set(GUROBI_FOUND TRUE)

