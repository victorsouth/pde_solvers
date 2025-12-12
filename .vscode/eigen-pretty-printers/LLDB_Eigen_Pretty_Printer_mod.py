import lldb
import re
import os
from functools import partial

def __lldb_init_module (debugger, dict):
    import sys
    prefix = sys.modules[__name__].__name__
    # Use the actual module name so LLDB can resolve functions correctly
    # LLDB uses POSIX regex, not Python regex - use .* instead of .*?
    debugger.HandleCommand('type summary add -x "^Eigen::Matrix<.*>$" -F %s.eigen_matrix_print -p -r -w Eigen' % prefix)
    debugger.HandleCommand('type summary add -x "^Eigen::PlainObjectBase<Eigen::Matrix<.*>$" -F %s.eigen_matrix_print -p -r -w Eigen' % prefix)
    debugger.HandleCommand('type summary add -x "^Eigen::Array<.*>$" -F %s.eigen_array_print -p -r -w Eigen' % prefix)
    debugger.HandleCommand('type summary add -x "^Eigen::PlainObjectBase<Eigen::Array<.*>$" -F %s.eigen_array_print -p -r -w Eigen' % prefix)
    debugger.HandleCommand('type summary add -x "^Eigen::Quaternion<.*>$" -F %s.eigen_quaternion_print -p -r -w Eigen' % prefix)
    debugger.HandleCommand('type summary add -x "^Eigen::SparseMatrix<.*>$" -F %s.eigen_sparsematrix_print -p -r -w Eigen' % prefix)
    debugger.HandleCommand('type synthetic add -x "^Eigen::Matrix<.*>$" --python-class %s.EigenMatrixSyntheticProvider -w Eigen' % prefix)
    debugger.HandleCommand('type synthetic add -x "^Eigen::PlainObjectBase<Eigen::Matrix<.*>$" --python-class %s.EigenMatrixSyntheticProvider -w Eigen' % prefix)
    debugger.HandleCommand('type synthetic add -x "^Eigen::Array<.*>$" --python-class %s.EigenArraySyntheticProvider -w Eigen' % prefix)
    debugger.HandleCommand('type synthetic add -x "^Eigen::PlainObjectBase<Eigen::Array<.*>$" --python-class %s.EigenArraySyntheticProvider -w Eigen' % prefix)
    debugger.HandleCommand('type synthetic add -x "^Eigen::Quaternion<.*>$" --python-class %s.EigenQuaternionSyntheticProvider -w Eigen' % prefix)
    debugger.HandleCommand('type synthetic add -x "^Eigen::SparseMatrix<.*>$" --python-class %s.EigenSparseMatrixSyntheticProvider -w Eigen' % prefix)
    debugger.HandleCommand("type category enable Eigen")

def evaluate_expression(valobj, expr):
    return valobj.GetProcess().GetSelectedThread().GetSelectedFrame().EvaluateExpression(expr)

def evaluate_at_index(valobj, index):
    return valobj.GetValueForExpressionPath("["+str(index)+"]")

class Printer:
    def __init__(self, data):
        self.data = data

    def evaluate_real(self, index):
        return "%1.8e" % float(evaluate_at_index(self.data, index).GetValue())

    def evaluate_bool(self, index):
        return "%d" % evaluate_at_index(self.data, index).GetValueAsUnsigned()

    def evaluate_complex_double(self, index):
        val = \
        list(self.data.GetValueForExpressionPath("["+str(index)+"]._M_value").GetValue())
        val[-1] = 'j'
        for n in range(0, len(val)):
            if (val[n] == " "):
                del val[n]
                del val[n+1]
                if val[n+1] == "-":
                    del val[n]
                break
        val = complex("".join(val))

        return '{0:1.5e} {1} {2:1.5e}i'.format(val.real,\
                                               '+-'[val.imag < 0],\
                                               abs(val.imag))

    def evaluate_complex_int(self, index):
        # val = self.data.GetValueForExpressionPath("["+str(index)+"]")
        val = evaluate_at_index(self.data, index)
        real = val.GetValueForExpressionPath("._M_real").GetValueAsSigned()
        imag = val.GetValueForExpressionPath("._M_imag").GetValueAsSigned()
        val = real + imag * 1j

        return '{0:1.5e} {1} {2:1.5e}i'.format(val.real,\
                                               '+-'[val.imag < 0],\
                                               abs(val.imag))

    def evaluate_complex_bool(self, index):
        # val = self.data.GetValueForExpressionPath("["+str(index)+"]")
        val = evaluate_at_index(self.data, index)
        real = val.GetValueForExpressionPath("._M_real").GetValueAsUnsigned()
        imag = val.GetValueForExpressionPath("._M_imag").GetValueAsUnsigned()

        return '{0:d} + {1:d}i'.format(real, imag)

class Matrix(Printer):
    def __init__(self, variety, val):
        try:
            valtype = val.GetType()
            type_name = valtype.GetName()
            
            # Handle PlainObjectBase - extract the template argument
            if "PlainObjectBase" in type_name:
                # Extract the Matrix/Array type from PlainObjectBase<Eigen::Matrix<...>>
                match = re.search(r'PlainObjectBase<(Eigen::(?:Matrix|Array)<[^>]+>)', type_name)
                if match:
                    type_str = match.group(1)
                else:
                    # Fallback: try to get from base class
                    if valtype.GetNumberOfDirectBaseClasses() > 0:
                        type_str = valtype.GetDirectBaseClassAtIndex(0).GetName()
                    else:
                        type_str = type_name
            else:
                # Try to get from base class first
                if valtype.GetNumberOfDirectBaseClasses() > 0:
                    type_str = valtype.GetDirectBaseClassAtIndex(0).GetName()
                else:
                    type_str = type_name

            begin = "Eigen::"+variety+"<"
            complex_scalar = "std::complex<"
            bool_key = "bool"

            if (type_str.find(complex_scalar) >= 0):
                regex = re.compile(begin + complex_scalar + ".*?>,.*?>")
                is_complex = True

                if (type_str.find(begin + complex_scalar + bool_key) >= 0):
                    is_bool = True
                else:
                    is_bool = False

            else:
                regex = re.compile(begin+".*?>")
                is_complex = False

                if (type_str.find(begin + bool_key) >= 0):
                    is_bool = True
                else:
                    is_bool = False

            self.variety = regex.findall(type_str)[0]
            m = self.variety[len(begin):-1]
            template_params = m.split(",")
            # template_params = [x.replace(">", "") for x in template_params]
            template_params = [x.replace(" ", "") for x in template_params]

            self.rows = int(template_params[1])
            if self.rows == -1:
                self.rows = val.GetValueForExpressionPath(".m_storage.m_rows").GetValueAsSigned()

            self.cols = int(template_params[2])
            if self.cols == -1:
                self.cols = val.GetValueForExpressionPath(".m_storage.m_cols").GetValueAsSigned()

            # Fix for vectors: if rows is 0 but cols is 1 (or vice versa), it's a vector
            if self.rows == 0 and self.cols == 1:
                # Column vector - try to get size from allocatedSize
                alloc = val.GetValueForExpressionPath(".m_storage.m_allocatedSize")
                if alloc.IsValid():
                    alloc_size = alloc.GetValueAsSigned()
                    if alloc_size > 0:
                        self.rows = alloc_size
            elif self.rows == 1 and self.cols == 0:
                # Row vector - try to get size from allocatedSize
                alloc = val.GetValueForExpressionPath(".m_storage.m_allocatedSize")
                if alloc.IsValid():
                    alloc_size = alloc.GetValueAsSigned()
                    if alloc_size > 0:
                        self.cols = alloc_size
            elif self.rows == 0 and self.cols == 0:
                # Both are 0 - try allocatedSize
                alloc = val.GetValueForExpressionPath(".m_storage.m_allocatedSize")
                if alloc.IsValid():
                    alloc_size = alloc.GetValueAsSigned()
                    if alloc_size > 0:
                        # Assume it's a column vector
                        self.rows = alloc_size
                        self.cols = 1

            self.options = 0
            if len(template_params) > 3:
                self.options = int(template_params[3])

            self.rowMajor = (int(self.options) & 0x1)
            self.innerType = template_params[0]

            if int(template_params[1]) == -1 or int(template_params[2]) == -1:
                data = val.GetValueForExpressionPath(".m_storage.m_data")
            else:
                data = val.GetValueForExpressionPath(".m_storage.m_data.array")

            Printer.__init__(self, data)
            if is_complex:
                val = self.data.GetValueForExpressionPath("[0]")
                if val.GetValueForExpressionPath("._M_value").IsValid():
                    self.get = partial(Printer.evaluate_complex_double, self)
                elif val.GetValueForExpressionPath("._M_real").IsValid():
                    if is_bool:
                        self.get = partial(Printer.evaluate_complex_bool, self)
                    else:
                        self.get = partial(Printer.evaluate_complex_int, self)
                else:
                    self.variety = -1
            else:
                if is_bool:
                    self.get = partial(Printer.evaluate_bool, self)
                else:
                    self.get = partial(Printer.evaluate_real, self)

        except:
            self.variety = -1

    def to_string(self):
        padding = 1
        for i in range(0, self.rows * self.cols):
            padding = max(padding, len(str(self.get(i))))

        output = "rows: %d, cols: %d\n[" % (self.rows, self.cols)
        if (self.options):
            for i in range(0, self.rows):
                if i!=0:
                    output += " "

                for j in range(0, self.cols):
                    val = self.get(i*self.cols+j)
                    if j!=0:
                        output += val.rjust(padding+2, ' ')
                    else:
                        output += val.rjust(padding+1, ' ')

                if i!=self.rows-1:
                    output += ";\n"
        else:
            for i in range(0, self.rows):
                if i!=0:
                    output += " "

                for j in range(0, self.cols):
                    val = self.get(i+j*self.rows)
                    if j!=0:
                        output += val.rjust(padding+2, ' ')
                    else:
                        output += val.rjust(padding+1, ' ')

                if i!=self.rows-1:
                    output += ";\n"

        output+=" ]\n"

        return output

class SparseMatrix(Printer):
    def __init__(self, val):
        try:
            valtype = val.GetType()
            type_str = valtype.GetDirectBaseClassAtIndex(0).GetName()

            begin = "Eigen::SparseMatrix<"
            complex_scalar = "std::complex<"
            bool_key = "bool"

            if (type_str.find(complex_scalar)>=0):
                regex = re.compile(begin + complex_scalar + ".*?>,.*?>")
                is_complex = True

                if (type_str.find(begin + complex_scalar + bool_key) >= 0):
                    is_bool = True
                else:
                    is_bool = False
            else:
                regex = re.compile(begin+".*?>")
                is_complex = False

                if (type_str.find(begin + bool_key) >= 0):
                    is_bool = True
                else:
                    is_bool = False

            self.variety = regex.findall(type_str)[0]

            m = self.variety[len(begin):-1]
            template_params = m.split(",")
            template_params = [x.replace(" ", "") for x in template_params]

            self.options = 0

            if len(template_params) > 1:
                self.options = int(template_params[1])

            self.rowMajor = (int(self.options) & 0x1)
            self.innerType = template_params[0]

            if self.rowMajor:
                self.rows = \
                val.GetValueForExpressionPath(".m_outerSize").GetValueAsSigned()
                self.cols = \
                val.GetValueForExpressionPath(".m_innerSize").GetValueAsSigned()
            else:
                self.rows = \
                val.GetValueForExpressionPath(".m_innerSize").GetValueAsSigned()
                self.cols = \
                val.GetValueForExpressionPath(".m_outerSize").GetValueAsSigned()

            self.outerStarts = val.GetValueForExpressionPath(".m_outerIndex")

            self.innerNNZs = []
            innerNNZs = val.GetValueForExpressionPath(".m_innerNonZeros")

            if (innerNNZs.GetValueAsSigned()):
                for k in range(0, self.rows if self.rowMajor else self.cols):
                    self.innerNNZs.append(evaluate_at_index(innerNNZs,
                                                            k).GetValueAsSigned())
            else:
                for k in range(0, self.rows if self.rowMajor else self.cols):
                    self.innerNNZs.append(evaluate_at_index(self.outerStarts,
                                                            k+1).GetValueAsSigned()
                                         - evaluate_at_index(self.outerStarts,
                                                            k).GetValueAsSigned())

            self.nonzeros = sum(self.innerNNZs)
            self.size = \
            val.GetValueForExpressionPath(".m_data.m_size").GetValueAsSigned()
            self.indices = val.GetValueForExpressionPath(".m_data.m_indices")

            data = val.GetValueForExpressionPath(".m_data.m_values")
            Printer.__init__(self, data)

            if is_complex:
                val = self.data.GetValueForExpressionPath("[0]")
                if val.GetValueForExpressionPath("._M_value").IsValid():
                    self.get = partial(Printer.evaluate_complex_double, self)
                elif val.GetValueForExpressionPath("._M_real").IsValid():
                    if is_bool:
                        self.get = partial(Printer.evaluate_complex_bool, self)
                    else:
                        self.get = partial(Printer.evaluate_complex_int, self)
                else:
                    self.variety = -1
            else:
                if is_bool:
                    self.get = partial(Printer.evaluate_bool, self)
                else:
                    self.get = partial(Printer.evaluate_real, self)
        except:
            self.variety = -1

    def to_string(self):
        padding = 1

        for i in range(0, self.size):
            padding = max(padding, len(str(self.get(i))))

        output = "rows: %d, cols: %d, nonzeros: %d\n{ " % (self.rows,
                                                           self.cols,
                                                           self.nonzeros)

        if (self.rowMajor):
            for i in range(0, self.rows):
                index = evaluate_at_index(self.outerStarts, \
                                          i).GetValueAsSigned()
                size = self.innerNNZs[i]

                for count in range(0, size):
                    j = evaluate_at_index(self.indices, \
                                          index+count).GetValueAsSigned()
                    output += "[%d, %d] =" % (i, j)
                    output += self.get(index+count).rjust(padding+1, ' ') + ", "
        else:
            rows = list()
            cols = list()
            vals = list()
            for j in range(0, self.cols):
                index = evaluate_at_index(self.outerStarts, \
                                          j).GetValueAsSigned()
                size = self.innerNNZs[j]

                for count in range(0, size):
                    i = evaluate_at_index(self.indices, \
                                          index+count).GetValueAsSigned()

                    val = self.get(index+count)
                    rows.append(i)
                    cols.append(j)
                    vals.append(val)

            indices = sorted(range(len(rows)), key=lambda k: rows[k])

            for index in indices:
                output += "[%d, %d] =" % (rows[index], cols[index])
                output += vals[index].rjust(padding+1, ' ') + ", "

        output += "\b\b }\n" if self.nonzeros else " }\n"

        return output


class Quaternion(Printer):
    def __init__(self, val):
        try:
            valtype = val.GetType()
            type_str = valtype.GetDirectBaseClassAtIndex(0).GetName()
            begin = "Eigen::Quaternion<"
            regex = re.compile(begin + ".*?>")

            self.variety = regex.findall(type_str)[0]

            m = self.variety[len(begin):-1]
            template_params = m.split(",")
            template_params = [x.replace(" ", "") for x in template_params]

            self.options = 0

            if len(template_params) > 1:
                self.options = int(template_params[1])

            self.DontAlign = (int(self.options) & 0x2)
            self.innerType = template_params[0]

            data = val.GetValueForExpressionPath(".m_coeffs.m_storage.m_data.array")
            Printer.__init__(self, data)
            self.variety = "Quaternion"
            self.get = partial(Printer.evaluate_real, self)
        except:
            self.variety = -1

    def to_string(self):
        padding = 1
        for i in range(0, 4):
            padding = max(padding, len(str(self.get(i))))

        output ="{ [x] = " + self.get(0) + ", " + \
                  "[y] = " + self.get(1) + ", " + \
                  "[z] = " + self.get(2) + ", " + \
                  "[w] = " + self.get(3) + " }\n"

        return output


class EigenMatrixSyntheticProvider:
    """Synthetic child provider for Eigen::Matrix to enable element-by-element inspection"""
    
    def __init__(self, valobj, internal_dict):
        self.valobj = valobj
        self.matrix = Matrix("Matrix", valobj)
        self.update()
    
    def num_children(self):
        if self.matrix.variety == -1:
            return 0
        return self.matrix.rows * self.matrix.cols
    
    def get_child_index(self, name):
        try:
            if name.startswith('[') and name.endswith(']'):
                indices = name[1:-1].split(',')
                if len(indices) == 1:
                    index = int(indices[0])
                    return index
                elif len(indices) == 2:
                    row = int(indices[0])
                    col = int(indices[1])
                    if self.matrix.rowMajor:
                        return row * self.matrix.cols + col
                    else:
                        return col * self.matrix.rows + row
        except:
            pass
        return None
    
    def get_child_at_index(self, index):
        if self.matrix.variety == -1 or index >= self.num_children():
            return None
        
        if self.matrix.rowMajor:
            row = index // self.matrix.cols
            col = index % self.matrix.cols
            linear_index = index
        else:
            row = index % self.matrix.rows
            col = index // self.matrix.rows
            linear_index = col * self.matrix.rows + row
        
        if self.matrix.cols == 1:
            name = '[%d]' % row
        elif self.matrix.rows == 1:
            name = '[%d]' % col
        else:
            name = '[%d,%d]' % (row, col)
        
        try:
            template_params = self.matrix.variety.split(',')
            rows_param = template_params[1].strip() if len(template_params) > 1 else "-1"
            cols_param = template_params[2].strip() if len(template_params) > 2 else "-1"
            
            is_fixed = (rows_param != "-1" and cols_param != "-1")
            
            if is_fixed:
                data_path = ".m_storage.m_data.array[%d]" % linear_index
            else:
                data_path = ".m_storage.m_data[%d]" % linear_index
            
            child = self.valobj.GetValueForExpressionPath(data_path)
            if child.IsValid():
                return child.CreateChildAtOffset(name, 0, child.GetType())
        except Exception as e:
            pass
        
        return None
    
    def update(self):
        pass


class EigenArraySyntheticProvider:
    """Synthetic child provider for Eigen::Array to enable element-by-element inspection"""
    
    def __init__(self, valobj, internal_dict):
        self.valobj = valobj
        self.array = Matrix("Array", valobj)
        self.update()
    
    def num_children(self):
        if self.array.variety == -1:
            return 0
        return self.array.rows * self.array.cols
    
    def get_child_index(self, name):
        try:
            if name.startswith('[') and name.endswith(']'):
                indices = name[1:-1].split(',')
                if len(indices) == 1:
                    index = int(indices[0])
                    return index
                elif len(indices) == 2:
                    row = int(indices[0])
                    col = int(indices[1])
                    if self.array.rowMajor:
                        return row * self.array.cols + col
                    else:
                        return col * self.array.rows + row
        except:
            pass
        return None
    
    def get_child_at_index(self, index):
        if self.array.variety == -1 or index >= self.num_children():
            return None
        
        if self.array.rowMajor:
            row = index // self.array.cols
            col = index % self.array.cols
            linear_index = index
        else:
            row = index % self.array.rows
            col = index // self.array.rows
            linear_index = col * self.array.rows + row
        
        if self.array.cols == 1:
            name = '[%d]' % row
        elif self.array.rows == 1:
            name = '[%d]' % col
        else:
            name = '[%d,%d]' % (row, col)
        
        try:
            template_params = self.array.variety.split(',')
            rows_param = template_params[1].strip() if len(template_params) > 1 else "-1"
            cols_param = template_params[2].strip() if len(template_params) > 2 else "-1"
            
            is_fixed = (rows_param != "-1" and cols_param != "-1")
            
            if is_fixed:
                data_path = ".m_storage.m_data.array[%d]" % linear_index
            else:
                data_path = ".m_storage.m_data[%d]" % linear_index
            
            child = self.valobj.GetValueForExpressionPath(data_path)
            if child.IsValid():
                return child.CreateChildAtOffset(name, 0, child.GetType())
        except Exception as e:
            pass
        
        return None
    
    def update(self):
        pass


class EigenQuaternionSyntheticProvider:
    """Synthetic child provider for Eigen::Quaternion to enable element-by-element inspection"""
    
    def __init__(self, valobj, internal_dict):
        self.valobj = valobj
        self.quaternion = Quaternion(valobj)
        self.update()
    
    def num_children(self):
        if self.quaternion.variety == -1:
            return 0
        return 4
    
    def get_child_index(self, name):
        element_names = ['x', 'y', 'z', 'w']
        if name.startswith('[') and name.endswith(']'):
            elem_name = name[1:-1]
            if elem_name in element_names:
                return element_names.index(elem_name)
        return None
    
    def get_child_at_index(self, index):
        if self.quaternion.variety == -1 or index >= 4:
            return None
        
        element_names = ['x', 'y', 'z', 'w']
        name = '[%s]' % element_names[index]
        
        try:
            data_path = ".m_coeffs.m_storage.m_data.array[%d]" % index
            child = self.valobj.GetValueForExpressionPath(data_path)
            if child.IsValid():
                return child.CreateChildAtOffset(name, 0, child.GetType())
        except:
            pass
        
        return None
    
    def update(self):
        pass


class EigenSparseMatrixSyntheticProvider:
    """Synthetic child provider for Eigen::SparseMatrix to enable element-by-element inspection"""
    
    def __init__(self, valobj, internal_dict):
        self.valobj = valobj
        self.sparsematrix = SparseMatrix(valobj)
        self.update()
    
    def num_children(self):
        if self.sparsematrix.variety == -1:
            return 0
        return self.sparsematrix.nonzeros
    
    def get_child_index(self, name):
        try:
            if name.startswith('[') and name.endswith(']'):
                indices = name[1:-1].split(',')
                if len(indices) == 2:
                    row = int(indices[0])
                    col = int(indices[1])
                    return None
        except:
            pass
        return None
    
    def get_child_at_index(self, index):
        if self.sparsematrix.variety == -1 or index >= self.num_children():
            return None
        
        try:
            if self.sparsematrix.rowMajor:
                row_idx = 0
                col_idx = 0
                current_idx = 0
                for i in range(self.sparsematrix.rows):
                    start = evaluate_at_index(self.sparsematrix.outerStarts, i).GetValueAsSigned()
                    size = self.sparsematrix.innerNNZs[i]
                    if current_idx + size > index:
                        row_idx = i
                        col_idx = evaluate_at_index(self.sparsematrix.indices, start + (index - current_idx)).GetValueAsSigned()
                        break
                    current_idx += size
            else:
                row_idx = 0
                col_idx = 0
                current_idx = 0
                for j in range(self.sparsematrix.cols):
                    start = evaluate_at_index(self.sparsematrix.outerStarts, j).GetValueAsSigned()
                    size = self.sparsematrix.innerNNZs[j]
                    if current_idx + size > index:
                        col_idx = j
                        row_idx = evaluate_at_index(self.sparsematrix.indices, start + (index - current_idx)).GetValueAsSigned()
                        break
                    current_idx += size
            
            name = '[%d,%d]' % (row_idx, col_idx)
            data_path = ".m_data.m_values[%d]" % index
            child = self.valobj.GetValueForExpressionPath(data_path)
            if child.IsValid():
                return child.CreateChildAtOffset(name, 0, child.GetType())
        except:
            pass
        
        return None
    
    def update(self):
        pass


def eigen_matrix_print (valobj,internal_dict):

    # if valobj.GetType().IsReferenceType():
        # val = valobj.Dereference()
    # else:
        # val = valobj

    matrix = Matrix("Matrix", valobj)

    if (matrix.variety == -1):
        return ""
    else:
        return matrix.to_string()

def eigen_array_print (valobj,internal_dict):

    array = Matrix("Array", valobj)

    if (array.variety == -1):
        return ""
    else:
        return array.to_string()

def eigen_quaternion_print (valobj,internal_dict):

    quaternion = Quaternion(valobj)

    if (quaternion.variety == -1):
        return ""
    else:
        return quaternion.to_string()

def eigen_sparsematrix_print (valobj,internal_dict):

    sparsematrix = SparseMatrix(valobj)

    if (sparsematrix.variety == -1):
        return ""
    else:
        return sparsematrix.to_string()
