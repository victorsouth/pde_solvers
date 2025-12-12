python
import sys, os
sys.path.insert(0, os.path.join(os.getcwd(), "eigen-pretty-printers"))
from printers import register_eigen_printers; register_eigen_printers(gdb)
end
