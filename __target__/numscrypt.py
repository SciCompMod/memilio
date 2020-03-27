# For performance reasons, real arrays or scalars and complex arrays can only be mixed in a limited way
# In general real arrays in natural order are fastest
# Real arrays in non-natural order are slower
# Complex arrays are slowest

from org.transcrypt.stubs.browser import __pragma__
    
__pragma__ ('skip')
Int32Array = Float32Array = Float64Array = Array = 0
__pragma__ ('noskip')
    
import itertools

ns_ctors = {
    'int32': Int32Array,
    'float32': Float32Array,
    'float64': Float64Array,
}

def ns_complex (dtype):
    return dtype in ('complex64', 'complex128')

def ns_buffertype (dtype):
    return (
            'float32'
        if dtype == 'complex64' else
            'float64'
        if dtype == 'complex128' else
            dtype
    )
    
def ns_complextype (dtype):
    return (
            'complex64'
        if dtype == 'float32' else
            'complex128'
        if dtype == 'float64' else
            None
    )
    
def ns_createbuf (imag, dtype, size):
    # The buffer will truly be created if it is the real part of a matrix or if that matrix is complex
    # So if the matrix is real and the buffer isn't the real part, it will not be created
    return (
            __new__ (ns_ctors [ns_buffertype (dtype)] (size))
        if not imag or ns_complex (dtype) else
            None
    )
    
class ndarray:
    def __init__ (
        self,
        shape,
        dtype,
        
        # Any fully constructed real array instance will have realbuf != None, imagbuf == None
        # Any fully constructed complex array instance will have realbuf != None, imagbuf != None
        realbuf = None,
        imagbuf = None
    ):
        self.dtype = dtype
        self.ns_complex = ns_complex (dtype)
        
        self.realbuf = realbuf
        if self.ns_complex:
            self.imagbuf = imagbuf
            
        self.setshape (shape)
        
    def setshape (self, shape):
        self.shape = shape
        self.ndim = shape.length
        self.ns_nrows = shape [0]
        
        if self.ndim == 1:
            self.size = self.ns_nrows
        else:
            self.ns_ncols = shape [1]
            self.size = self.ns_nrows * self.ns_ncols
        
    def astype (self, dtype):   # Do not use to convert between real and complex arrays
        result = empty (self.shape, dtype)
        
        result.realbuf.set (self.realbuf)
        if self.ns_complex:
            result.imagbuf.set (self.imagbuf)
            
        return result
                        
    def tolist (self):
        if self.ns_complex:
            flat = [complex (real, imag) for real, imag in zip (list (self.realbuf), list (self.imagbuf))]
        else:
            flat = self.realbuf
            
        if self.ndim == 1:
            return list (flat)
        else:
            return [[flat [self.ns_ncols * irow + icol] for icol in range (self.ns_ncols)] for irow in range (self.ns_nrows)]
                    
    def __repr__ (self):
        return 'array({})'.format (repr (self.tolist ()))

    def __str__ (self):
        if self.ndim == 1:
            return str (self.tolist ())
        else:
            return '[\n\t{}\n]\n'.format ('\n\t'.join ([str (row) for row in self.tolist ()]))
        
    def reshape (self, shape):
        if self.ndim == 1:
            return array (self, self.dtype),
        else:
            result = array (self, self.dtype)
            result.setshape (self.ns_ncols, self.ns_nrows)
            return result
            
    def transpose (self):
        if self.ndim == 1:
            result = array (self, dtype)
        else:
            result = empty ((self.ns_ncols, self.ns_nrows), self.dtype)
                
            itarget = 0
            if self.ns_complex:
                for icol in range (self.ns_ncols):
                    isource = icol
                    for irow in range (self.ns_nrows):
                        isource = self.ns_ncols * irow + icol
                        result.imagbuf [itarget] = self.imagbuf [isource]
                        result.realbuf [__postinc__ (itarget)] = self.realbuf [isource]
                        isource += self.ns_ncols
            else:
                for icol in range (self.ns_ncols):
                    isource = icol
                    for irow in range (self.ns_nrows):
                        result.realbuf [__postinc__ (itarget)] = self.realbuf [isource]
                        isource += self.ns_ncols
                    
        return result
                
    def __getitem__ (self, key):
        if self.ndim == 1:          
            if type (key) == tuple:
                # Slice of single dim array

                if key [1] == None:
                    key [1] = self.size
                elif key [1] < 0:
                    key [1] += self.size
                
                result = empty ([(key [1] - key [0]) / key [2]], self.dtype)
                
                itarget = 0
                if self.ns_complex:
                    for isource in range (*self.shape):
                        result.realbuf [itarget] = self.realbuf [isource]
                        result.imagbuf [__postinc__ (itarget)] = self.imagbuf [isource]
                else:
                    for isource in range (*self.shape):
                        result.realbuf [__postinc__ (itarget)] = self.realbuf [isource]
                
                return result
            else:
                # Element of single dim array
            
                if self.ns_complex:
                    return complex (self.realbuf [key], self.imagbuf [key])
                else:
                    return self.realbuf [key]
        else:
            rowkey = key [0]
            colkey = key [1]
            
            rowistup = type (rowkey) == tuple
            colistup = type (colkey) == tuple
            
            if rowistup:
                if rowkey [1] == None:
                    rowkey [1] = self.ns_nrows
                elif rowkey [1] < 0:
                    rowkey [1] += self.ns_nrows
                    
            if colistup:
                if colkey [1] == None:
                    colkey [1] = self.ns_ncols
                elif colkey [1] < 0:
                    colkey [1] += self.ns_ncols
            
            if rowistup or colistup:
                # Slice of multidim array
            
                if not rowistup:
                    result = empty (((colkey [1] - colkey [0]) / colkey [2], ), self.dtype)
                    
                    itarget = 0
                    if self.ns_complex:
                        for isourcecol in range (*colkey):
                            isource = self.ns_ncols * rowkey + isourcecol
                            result.realbuf [itarget] = self.realbuf [isource]
                            result.imagbuf [__postinc__ (itarget)] = self.imagbuf [isource]
                    else:
                        for isourcecol in range (*colkey):
                            result.realbuf [__postinc__ (itarget)] = self.realbuf [self.ns_ncols * rowkey + isourcecol]
                elif not colistup:
                    result = empty (((rowkey [1] - rowkey [0]) / rowkey [2], ), self.dtype)
                    
                    itarget = 0
                    if self.ns_complex:
                        for isourcerow in range (*rowkey):
                            isource = self.ns_ncols * isourcerow + colkey
                            result.realbuf [itarget] = self.realbuf [isource]
                            result.imagbuf [__postinc__ (itarget)] = self.imagbuf [isource]
                    else:
                        for isourcerow in range (*rowkey):
                            result.realbuf [__postinc__ (itarget)] = self.realbuf [self.ns_ncols * isourcerow + colkey]
                else:
                    result = empty ((
                        (key[0][1] - key[0][0]) / key [0][2],
                        (key[1][1] - key[1][0]) / key [1][2]
                    ), self.dtype)
                                
                    itarget = 0
                    if self.ns_complex:
                        for isourcerow in range (*rowkey):
                            for isourcecol in range (*colkey):
                                isource = self.ns_ncols * isourcerow + isourcecol
                                result.realbuf [itarget] = self.realbuf [isource]
                                result.imagbuf [__postinc__ (itarget)] = self.imagbuf [isource]                         

                    else:
                        for isourcerow in range (*rowkey):
                            for isourcecol in range (*colkey):
                                result.realbuf [__postinc__ (itarget)] = self.realbuf [self.ns_ncols * isourcerow + isourcecol]
                            
                return result
            else:
                # Element of multi dim array
            
                if self.ns_complex:
                    isource = self.ns_ncols * key [0] + key [1]
                    return complex (self.realbuf [isource], self.imagbuf [isource])
                else:
                    return self.realbuf [self.ns_ncols * key [0] + key [1]]
                    
    def __setitem__ (self, key, value):
        if self.ndim == 1:      
            if type (key) == tuple:
                # Slice of single dim array
                
                if key [1] == None:
                    key [1] = self.size
                elif key [1] < 0:
                    key [1] += self.size
                
                isource = 0
                if self.ns_complex:
                    for itarget in range (*self.shape):
                        self.realbuf [itarget] = value.realbuf [isource]
                        self.imagbuf [itarget] = value.imagbuf [__postinc__ (isource)]
                else:
                    for itarget in range (*self.shape):
                        self.realbuf [itarget] = value.realbuf [__postinc__ (isource)]
                
                return result
            else:
                # Element of single dim array
                
                if self.ns_complex:
                    if __typeof__ (value) == 'number':
                        self.realbuf [key] = value
                        self.imagbuf [key] = 0
                    else:                   
                        self.realbuf [key] = value.real
                        self.imagbuf [key] = value.imag
                else:
                    self.realbuf [key] = value
        else:
            rowkey = key [0]
            colkey = key [1]
            
            rowistup = type (rowkey) == tuple
            colistup = type (colkey) == tuple
            
            if rowistup:
                if rowkey [1] == None:
                    rowkey [1] = self.ns_nrows
                elif rowkey [1] < 0:
                    rowkey [1] += self.ns_nrows
                    
            if colistup:
                if colkey [1] == None:
                    colkey [1] = self.ns_ncols
                elif colkey [1] < 0:
                    colkey [1] += self.ns_ncols
                    
            if rowistup or colistup:
                # Slice of multi dim array
                
                if not rowistup:            
                    isource = 0
                    if self.ns_complex:
                        for itargetcol in range (*colkey):
                            itarget = self.ns_ncols * rowkey + itargetcol
                            self.realbuf [itarget] = value.realbuf [isource]
                            self.imagbuf [itarget] = value.imagbuf [__postinc__ (isource)]
                    else:
                        for itargetcol in range (*colkey):
                            result.realbuf [self.ns_ncols * rowkey + itargetcol] = self.realbuf [__postinc__ (isource)]
                elif not colistup:
                    isource = 0
                    if self.ns_complex:
                        for itargetrow in range (*rowkey):
                            itarget = self.ns_ncols * itargetrow + colkey
                            self.realbuf [itarget] = value.realbuf [isource]
                            self.imagbuf [itarget] = value.imagbuf [__postinc__ (isource)]
                    else:
                        for isourcerow in range (*rowkey):
                            self.realbuf [self.ns_ncols * isourcerow + colkey] = value [__postinc__ (isource)]
                else:           
                    isource = 0
                    if self.ns_complex:
                        for itargetrow in range (*rowkey):
                            for itargetcol in range (*colkey):
                                itarget = self.ns_ncols * itargetrow + itargetcol
                                self.realbuf [itarget] = value.realbuf [isource]
                                self.imagbuf [itarget] = value.imagbuf [__postinc__ (isource)]                          
                    else:
                        for isourcerow in range (*rowkey):
                            for isourcecol in range (*colkey):
                                self.realbuf [self.ns_ncols * itargetrow + itargetcol] = value.realbuf [__postinc__ (isource)]                      
            else:
                # Element of multi dim array
                
                if self.ns_complex:
                    itarget = self.ns_ncols * key [0] + key [1]
                    
                    if __typeof__ (value) == 'number':
                        self.realbuf [itarget] = value
                        self.imagbuf [itarget] = 0
                    else:
                        self.realbuf [itarget] = value.real
                        self.imagbuf [itarget] = value.imag
                    
                else:
                    self.realbuf [self.ns_ncols * key [0] + key [1]] = value
              
    def real (self):    # Returns a view, so you can assign to self via it
        return ndarray (self.shape, ns_buffertype (self.dtype), self.realbuf)
    
    def imag (self):    # Returns a view, so you can assign self via it
        return ndarray (self.shape, ns_buffertype (self.dtype), self.imagbuf)
        
    def __conj__ (self):
        if self.ns_complex:
            result = empty (self.shape, self.dtype)
            result.realbuf.set (self.realbuf)
            result.imagbuf = ns_createbuf (True, self.dtype, self.size)
            for i in range (self.size):
                result.imagbuf [i] = -self.imagbuf [i]
            return result
        else:
            return copy (self)
        
    def conjugate (self):
        return self.__conj__ ()
        
    def __neg__ (self):
        result = empty (self.shape, self.dtype)
        if self.ns_complex:
            for i in range (self.size):
                result.realbuf [i] = -self.realbuf [i]
                result.imagbuf [i] = -self.imagbuf [i]
        else:
            for i in range (self.size):
                result.realbuf [i] = -self.realbuf [i]
            
        return result   
            
    def __ns_inv__ (self):
        result = empty (self.shape, self.dtype)
            
        if self.ns_complex:
            for i in range (self.size):
                real = self.realbuf [i]
                imag = self.imagbuf [i]
                denom = real * real + imag * imag
                
                result.realbuf [i] = real / denom
                result.imagbuf [i] = -imag / denom
        else:
            for i in range (self.size):
                result.realbuf [i] = 1 / self.realbuf [i]
            
        return result   
            
    def __add__ (self, other):
        result = empty (self.shape, self.dtype)
        
        if type (other) == ndarray:
            if self.ns_complex:
                for i in range (self.size):
                    result.realbuf [i] = self.realbuf [i] + other.realbuf [i]
                    result.imagbuf [i] = self.imagbuf [i] + other.imagbuf [i]
            else:
                for i in range (self.size):
                    result.realbuf [i] = self.realbuf [i] + other.realbuf [i]
        else:
            if self.ns_complex:
                for i in range (self.size):
                    result.realbuf [i] = self.realbuf [i] + other.real
                    result.imagbuf [i] = self.imagbuf [i] + other.imag
            else:
                for i in range (self.size):
                    result.realbuf [i] = self.realbuf [i] + other
                    
        return result
        
    def __radd__ (self, scalar):    # scalar + array -> array.__radd__ (scalar)
        return self.__add__ (scalar)
        
    def __sub__ (self, other):
        result = empty (self.shape, self.dtype)
        
        if type (other) == ndarray:
            if self.ns_complex:
                for i in range (self.size):
                    result.realbuf [i] = self.realbuf [i] - other.realbuf [i]
                    result.imagbuf [i] = self.imagbuf [i] - other.imagbuf [i]
            else:
                for i in range (self.size):
                    result.realbuf [i] = self.realbuf [i] - other.realbuf [i]
        else:
            if self.ns_complex:
                for i in range (self.size):
                    result.realbuf [i] = self.realbuf [i] - other.real
                    result.imagbuf [i] = self.imagbuf [i] - other.imag
            else:
                for i in range (self.size):
                    result.realbuf [i] = self.realbuf [i] - other
                    
        return result
        
    def __rsub__ (self, scalar):    # scalar - array -> array.__rsub__ (scalar)
        return self.__neg__ () .__add__ (scalar)
        
    def __mul__ (self, other):
        result = empty (self.shape, self.dtype)
        
        if type (other) == ndarray:
            if self.ns_complex:
                for i in range (self.size):
                    result.realbuf [i] = self.realbuf [i] * other.realbuf [i] - self.imagbuf [i] * other.imagbuf [i]
                    result.imagbuf [i] = self.realbuf [i] * other.imagbuf [i] + self.imagbuf [i] * other.realbuf [i]
            else:
                for i in range (self.size):
                    result.realbuf [i] = self.realbuf [i] * other.realbuf [i]
        else:
            if self.ns_complex:
                for i in range (self.size):
                    result.realbuf [i] = self.realbuf [i] * other.real - self.imagbuf [i] * other.imag
                    result.imagbuf [i] = self.realbuf [i] * other.imag + self.imagbuf [i] * other.real
            else:
                for i in range (self.size):
                    result.realbuf [i] = self.realbuf [i] * other
                    
        return result
        
    def __rmul__ (self, scalar):    # scalar * array -> array.__rmul__ (scalar)
        return self.__mul__ (scalar)
        
    def __div__ (self, other):
        result = empty (self.shape, self.dtype)
        
        if type (other) == ndarray:
            if self.ns_complex:
                for i in range (self.size):
                    real = other.realbuf [i]
                    imag = other.imagbuf [i]
                    denom = real * real + imag * imag
                
                    result.realbuf [i] = (self.realbuf [i] * real + self.imagbuf [i] * imag) / denom
                    result.imagbuf [i] = (self.imagbuf [i] * real - self.realbuf [i] * imag) / denom
            else:
                for i in range (self.size):
                    result.realbuf [i] = self.realbuf [i] / other.realbuf [i]
        else:
            if self.ns_complex:
                real = other.real
                imag = other.imag
                denom = real * real + imag * imag
                
                for i in range (self.size):
                    result.realbuf [i] = (self.realbuf [i] * real + self.imagbuf [i] * imag) / denom
                    result.imagbuf [i] = (self.imagbuf [i] * real - self.realbuf [i] * imag) / denom
            else:
                for i in range (self.size):
                    result.realbuf [i] = self.realbuf [i] / other
                    
        return result
        
    def __rdiv__ (self, scalar):    # scalar / array -> array.__rdiv__ (scalar)
        return self.__ns_inv__ () .__mul__ (scalar)
        
    def __matmul__ (self, other):
        result = empty ((self.ns_nrows, other.ns_ncols), self.dtype)
        
        if self.ns_complex:
            iresult = 0
            for irow in range (self.ns_nrows):
                for icol in range (other.ns_ncols):
                    result.realbuf [iresult] = 0
                    result.imagbuf [iresult] = 0
                    iself = self.ns_ncols * irow
                    iother = icol
                    for iterm in range (self.ns_ncols):
                        result.realbuf [iresult] += self.realbuf [iself] * other.realbuf [iother] - self.imagbuf [iself] * other.imagbuf [iother]
                        result.imagbuf [iresult] += self.realbuf [iself] * other.imagbuf [iother] + self.imagbuf [__postinc__ (iself)] * other.realbuf [iother]
                        iother += other.ns_ncols
                    iresult += 1
        else:
            iresult = 0
            for irow in range (self.ns_nrows):
                for icol in range (other.ns_ncols):
                    result.realbuf [iresult] = 0
                    iself = self.ns_ncols * irow
                    iother = icol
                    for iterm in range (self.ns_ncols):
                        result.realbuf [iresult] += self.realbuf [__postinc__ (iself)] * other.realbuf [iother]
                        iother += other.ns_ncols
                    iresult += 1
            
        return result
        
def empty (shape, dtype = 'float64'):
    result = ndarray (
        shape,
        dtype
    )
    result.realbuf = ns_createbuf (False, dtype, result.size)
    result.imagbuf = ns_createbuf (True, dtype, result.size)
    return result
    
def array (obj, dtype = 'float64'): 
    if Array.isArray (obj):
        if len (obj):
            if Array.isArray (obj [0]):
                result = empty ((obj.length, obj [0] .length), dtype)
                iresult = 0
                if result.ns_complex:
                    for irow in range (result.ns_nrows):
                        for icol in range (result.ns_ncols):
                            element = complex (obj [irow][icol])
                            result.realbuf [iresult] = element.real
                            result.imagbuf [__postinc__ (iresult)] = element.imag
                else:
                    for irow in range (result.ns_nrows):
                        for icol in range (result.ns_ncols):
                            result.realbuf [__postinc__ (iresult)] = obj [irow][icol]
            else:
                result = empty ((obj.length, ), dtype)
                if result.ns_complex:
                    for i in range (result.size):
                        element = complex (obj [i])
                        result.realbuf [i] = element.real 
                        result.imagbuf [i] = element.imag
                else:
                    for i in range (result.size):
                        result.realbuf [i] = obj [i]
        else:
            result = empty ((0, ), dtype)
    else:   # Assume obj is an ndarray        
        result = empty (obj.shape, obj.dtype)
        result.realbuf.set (obj.realbuf)        
        if obj.ns_complex:
            result.imagbuf.set (obj.imagbuf)       
    return result
        
def copy (obj):
    return array (obj, obj.dtype)
    
def hsplit (ary, nparts):
    result = [empty ((ary.ns_nrows, ary.ns_ncols / nparts), ary.dtype) for ipart in range (nparts)]
    
    isource = 0
    if ary.ns_complex:
        for irow in range (ary.ns_nrows):
            for part in result:
                itarget = part.ns_ncols * irow
                for icol in range (part.ns_ncols):
                    part.realbuf [itarget] = ary.realbuf [isource]
                    part.imagbuf [__postinc__ (itarget)] = ary.imagbuf [__postinc__ (isource)]
    else:
        for irow in range (ary.ns_nrows):
            for part in result:
                itarget = part.ns_ncols * irow
                for icol in range (part.ns_ncols):
                    part.realbuf [__postinc__ (itarget)] = ary.realbuf [__postinc__ (isource)]
                    
    return result
    
def vsplit (ary, nparts):
    result = [empty ((ary.ns_nrows / nparts, ary.ns_ncols), array.dtype) for ipart in range (nparts)]

    isource = 0
    if ary.ns_complex:
        for part in result:
            for itarget in range (part.size):
                part.realbuf [itarget] = ary.realbuf [isource]
                part.imagbuf [itarget] = ary.imagbuf [__postinc__ (isource)]
    else:
        for part in result:
            for itarget in range (part.size):
                part.realbuf [itarget] = ary.realbuf [__postinc__ (isource)]
                    
    return result
    
def hstack (tup):
    ncols = 0
    for part in tup:
        ncols += part.ns_ncols
                
    result = empty ((tup [0] .ns_nrows, ncols), tup [0] .dtype)
    
    itarget = 0
    if result.ns_complex:
        for irow in range (result.ns_nrows):
            for part in tup:
                isource = part.ns_ncols * irow
                for icol in range (part.ns_ncols):
                    result.realbuf [itarget] = part.realbuf [isource]
                    result.imagbuf [__postinc__ (itarget)] = part.imagbuf [__postinc__ (isource)]
    else:
        for irow in range (result.ns_nrows):
            for part in tup:
                isource = part.ns_ncols * irow
                for icol in range (part.ns_ncols):
                    result.realbuf [__postinc__ (itarget)] = part.realbuf [__postinc__ (isource)]
                
    return result
    
def vstack (tup):
    nrows = 0
    for part in tup:
        nrows += part.ns_nrows
        
    result = empty ((nrows, tup [0].ns_ncols), tup [0] .dtype)
    
    itarget = 0
    if result.ns_complex:
        for part in tup:
            for isource in range (part.size):
                result.realbuf [itarget] = part.realbuf [isource]
                result.imagbuf [__postinc__ (itarget)] = part.imagbuf [isource]
    else:
        for part in tup:
            for isource in range (part.size):
                result.realbuf [__postinc__ (itarget)] = part.realbuf [isource]
                
    return result
            
def round (a, decimals = 0):    # Truncation rather than bankers rounding, for speed
    result = empty (a.shape, a.dtype)

    if a.ns_complex:
        for i in range (a.size):
            result.realbuf [i] = a.realbuf [i] .toFixed (decimals)
            result.imagbuf [i] = a.imagbuf [i] .toFixed (decimals)
    else:
        for i in range (a.size):
            result.realbuf [i] = a.realbuf [i] .toFixed (decimals)
        
    return result
        
def zeros (shape, dtype = 'float64'):
    result = empty (shape, dtype)

    if result.ns_complex:
        for i in range (result.size):
            result.realbuf [i] = 0
            result.imagbuf [i] = 0
    else:
        for i in range (result.size):
            result.realbuf [i] = 0
            
    return result
    
def ones (shape, dtype = 'float64'):
    result = empty (shape, dtype)
    
    if result.ns_complex:
        for i in range (result.size):
            result.realbuf [i] = 1
            result.imagbuf [i] = 0
    else:
        for i in range (result.size):
            result.realbuf [i] = 1
            
    return result
    
def identity (n, dtype = 'float64'):
    result = zeros ((n, n), dtype)
    
    i = 0
    shift = n + 1
    for j in range (n):
        result.realbuf [i] = 1
        i += shift
    
    return result
    
def conjugate (x):
    return x.__conj__ ()
        
