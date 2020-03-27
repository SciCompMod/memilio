// Transcrypt'ed from Python, 2020-03-27 23:22:05
var itertools = {};
import {AssertionError, AttributeError, BaseException, DeprecationWarning, Exception, IndexError, IterableError, KeyError, NotImplementedError, RuntimeWarning, StopIteration, UserWarning, ValueError, Warning, __JsIterator__, __PyIterator__, __Terminal__, __add__, __and__, __call__, __class__, __conj__, __envir__, __eq__, __floordiv__, __ge__, __get__, __getcm__, __getitem__, __getslice__, __getsm__, __gt__, __i__, __iadd__, __iand__, __idiv__, __ijsmod__, __ilshift__, __imatmul__, __imod__, __imul__, __in__, __init__, __ior__, __ipow__, __irshift__, __isub__, __ixor__, __jsUsePyNext__, __jsmod__, __k__, __kwargtrans__, __le__, __lshift__, __lt__, __matmul__, __mergefields__, __mergekwargtrans__, __mod__, __mul__, __ne__, __neg__, __nest__, __or__, __pow__, __pragma__, __proxy__, __pyUseJsNext__, __rshift__, __setitem__, __setproperty__, __setslice__, __sort__, __specialattrib__, __sub__, __super__, __t__, __terminal__, __truediv__, __withblock__, __xor__, abs, all, any, assert, bool, bytearray, bytes, callable, chr, complex, deepcopy, delattr, dict, dir, divmod, enumerate, filter, float, getattr, hasattr, input, int, isinstance, issubclass, len, list, map, max, min, object, ord, pow, print, property, py_TypeError, py_iter, py_metatype, py_next, py_reversed, py_typeof, range, repr, set, setattr, sorted, str, sum, tuple, zip} from './org.transcrypt.__runtime__.js';
import * as __module_itertools__ from './itertools.js';
__nest__ (itertools, '', __module_itertools__);
var __name__ = 'numscrypt';
export var ns_ctors = dict ({'int32': Int32Array, 'float32': Float32Array, 'float64': Float64Array});
export var ns_complex = function (dtype) {
	return __in__ (dtype, tuple (['complex64', 'complex128']));
};
export var ns_buffertype = function (dtype) {
	return (dtype == 'complex64' ? 'float32' : (dtype == 'complex128' ? 'float64' : dtype));
};
export var ns_complextype = function (dtype) {
	return (dtype == 'float32' ? 'complex64' : (dtype == 'float64' ? 'complex128' : null));
};
export var ns_createbuf = function (imag, dtype, size) {
	return (!(imag) || ns_complex (dtype) ? new ns_ctors [ns_buffertype (dtype)] (size) : null);
};
export var ndarray =  __class__ ('ndarray', [object], {
	__module__: __name__,
	get __init__ () {return __get__ (this, function (self, shape, dtype, realbuf, imagbuf) {
		if (typeof realbuf == 'undefined' || (realbuf != null && realbuf.hasOwnProperty ("__kwargtrans__"))) {;
			var realbuf = null;
		};
		if (typeof imagbuf == 'undefined' || (imagbuf != null && imagbuf.hasOwnProperty ("__kwargtrans__"))) {;
			var imagbuf = null;
		};
		self.dtype = dtype;
		self.ns_complex = ns_complex (dtype);
		self.realbuf = realbuf;
		if (self.ns_complex) {
			self.imagbuf = imagbuf;
		}
		self.setshape (shape);
	});},
	get setshape () {return __get__ (this, function (self, shape) {
		self.shape = shape;
		self.ndim = shape.length;
		self.ns_nrows = shape [0];
		if (self.ndim == 1) {
			self.size = self.ns_nrows;
		}
		else {
			self.ns_ncols = shape [1];
			self.size = self.ns_nrows * self.ns_ncols;
		}
	});},
	get astype () {return __get__ (this, function (self, dtype) {
		var result = empty (self.shape, dtype);
		result.realbuf.set (self.realbuf);
		if (self.ns_complex) {
			result.imagbuf.set (self.imagbuf);
		}
		return result;
	});},
	get tolist () {return __get__ (this, function (self) {
		if (self.ns_complex) {
			var flat = (function () {
				var __accu0__ = [];
				for (var [real, imag] of zip (list (self.realbuf), list (self.imagbuf))) {
					__accu0__.append (complex (real, imag));
				}
				return __accu0__;
			}) ();
		}
		else {
			var flat = self.realbuf;
		}
		if (self.ndim == 1) {
			return list (flat);
		}
		else {
			return (function () {
				var __accu0__ = [];
				for (var irow = 0; irow < self.ns_nrows; irow++) {
					__accu0__.append ((function () {
						var __accu1__ = [];
						for (var icol = 0; icol < self.ns_ncols; icol++) {
							__accu1__.append (flat [self.ns_ncols * irow + icol]);
						}
						return __accu1__;
					}) ());
				}
				return __accu0__;
			}) ();
		}
	});},
	get __repr__ () {return __get__ (this, function (self) {
		return 'array({})'.format (repr (self.tolist ()));
	});},
	get __str__ () {return __get__ (this, function (self) {
		if (self.ndim == 1) {
			return str (self.tolist ());
		}
		else {
			return '[\n\t{}\n]\n'.format ('\n\t'.join ((function () {
				var __accu0__ = [];
				for (var row of self.tolist ()) {
					__accu0__.append (str (row));
				}
				return __accu0__;
			}) ()));
		}
	});},
	get reshape () {return __get__ (this, function (self, shape) {
		if (self.ndim == 1) {
			return tuple ([array (self, self.dtype)]);
		}
		else {
			var result = array (self, self.dtype);
			result.setshape (self.ns_ncols, self.ns_nrows);
			return result;
		}
	});},
	get transpose () {return __get__ (this, function (self) {
		if (self.ndim == 1) {
			var result = array (self, dtype);
		}
		else {
			var result = empty (tuple ([self.ns_ncols, self.ns_nrows]), self.dtype);
			var itarget = 0;
			if (self.ns_complex) {
				for (var icol = 0; icol < self.ns_ncols; icol++) {
					var isource = icol;
					for (var irow = 0; irow < self.ns_nrows; irow++) {
						var isource = self.ns_ncols * irow + icol;
						result.imagbuf [itarget] = self.imagbuf [isource];
						result.realbuf [itarget++] = self.realbuf [isource];
						isource += self.ns_ncols;
					}
				}
			}
			else {
				for (var icol = 0; icol < self.ns_ncols; icol++) {
					var isource = icol;
					for (var irow = 0; irow < self.ns_nrows; irow++) {
						result.realbuf [itarget++] = self.realbuf [isource];
						isource += self.ns_ncols;
					}
				}
			}
		}
		return result;
	});},
	get __getitem__ () {return __get__ (this, function (self, key) {
		if (self.ndim == 1) {
			if (py_typeof (key) == tuple) {
				if (key [1] == null) {
					key [1] = self.size;
				}
				else if (key [1] < 0) {
					key [1] += self.size;
				}
				var result = empty ([(key [1] - key [0]) / key [2]], self.dtype);
				var itarget = 0;
				if (self.ns_complex) {
					for (var isource of range (...self.shape)) {
						result.realbuf [itarget] = self.realbuf [isource];
						result.imagbuf [itarget++] = self.imagbuf [isource];
					}
				}
				else {
					for (var isource of range (...self.shape)) {
						result.realbuf [itarget++] = self.realbuf [isource];
					}
				}
				return result;
			}
			else if (self.ns_complex) {
				return complex (self.realbuf [key], self.imagbuf [key]);
			}
			else {
				return self.realbuf [key];
			}
		}
		else {
			var rowkey = key [0];
			var colkey = key [1];
			var rowistup = py_typeof (rowkey) == tuple;
			var colistup = py_typeof (colkey) == tuple;
			if (rowistup) {
				if (rowkey [1] == null) {
					rowkey [1] = self.ns_nrows;
				}
				else if (rowkey [1] < 0) {
					rowkey [1] += self.ns_nrows;
				}
			}
			if (colistup) {
				if (colkey [1] == null) {
					colkey [1] = self.ns_ncols;
				}
				else if (colkey [1] < 0) {
					colkey [1] += self.ns_ncols;
				}
			}
			if (rowistup || colistup) {
				if (!(rowistup)) {
					var result = empty (tuple ([(colkey [1] - colkey [0]) / colkey [2]]), self.dtype);
					var itarget = 0;
					if (self.ns_complex) {
						for (var isourcecol of range (...colkey)) {
							var isource = self.ns_ncols * rowkey + isourcecol;
							result.realbuf [itarget] = self.realbuf [isource];
							result.imagbuf [itarget++] = self.imagbuf [isource];
						}
					}
					else {
						for (var isourcecol of range (...colkey)) {
							result.realbuf [itarget++] = self.realbuf [self.ns_ncols * rowkey + isourcecol];
						}
					}
				}
				else if (!(colistup)) {
					var result = empty (tuple ([(rowkey [1] - rowkey [0]) / rowkey [2]]), self.dtype);
					var itarget = 0;
					if (self.ns_complex) {
						for (var isourcerow of range (...rowkey)) {
							var isource = self.ns_ncols * isourcerow + colkey;
							result.realbuf [itarget] = self.realbuf [isource];
							result.imagbuf [itarget++] = self.imagbuf [isource];
						}
					}
					else {
						for (var isourcerow of range (...rowkey)) {
							result.realbuf [itarget++] = self.realbuf [self.ns_ncols * isourcerow + colkey];
						}
					}
				}
				else {
					var result = empty (tuple ([(key [0] [1] - key [0] [0]) / key [0] [2], (key [1] [1] - key [1] [0]) / key [1] [2]]), self.dtype);
					var itarget = 0;
					if (self.ns_complex) {
						for (var isourcerow of range (...rowkey)) {
							for (var isourcecol of range (...colkey)) {
								var isource = self.ns_ncols * isourcerow + isourcecol;
								result.realbuf [itarget] = self.realbuf [isource];
								result.imagbuf [itarget++] = self.imagbuf [isource];
							}
						}
					}
					else {
						for (var isourcerow of range (...rowkey)) {
							for (var isourcecol of range (...colkey)) {
								result.realbuf [itarget++] = self.realbuf [self.ns_ncols * isourcerow + isourcecol];
							}
						}
					}
				}
				return result;
			}
			else if (self.ns_complex) {
				var isource = self.ns_ncols * key [0] + key [1];
				return complex (self.realbuf [isource], self.imagbuf [isource]);
			}
			else {
				return self.realbuf [self.ns_ncols * key [0] + key [1]];
			}
		}
	});},
	get __setitem__ () {return __get__ (this, function (self, key, value) {
		if (self.ndim == 1) {
			if (py_typeof (key) == tuple) {
				if (key [1] == null) {
					key [1] = self.size;
				}
				else if (key [1] < 0) {
					key [1] += self.size;
				}
				var isource = 0;
				if (self.ns_complex) {
					for (var itarget of range (...self.shape)) {
						self.realbuf [itarget] = value.realbuf [isource];
						self.imagbuf [itarget] = value.imagbuf [isource++];
					}
				}
				else {
					for (var itarget of range (...self.shape)) {
						self.realbuf [itarget] = value.realbuf [isource++];
					}
				}
				return result;
			}
			else if (self.ns_complex) {
				if (typeof value == 'number') {
					self.realbuf [key] = value;
					self.imagbuf [key] = 0;
				}
				else {
					self.realbuf [key] = value.real;
					self.imagbuf [key] = value.imag;
				}
			}
			else {
				self.realbuf [key] = value;
			}
		}
		else {
			var rowkey = key [0];
			var colkey = key [1];
			var rowistup = py_typeof (rowkey) == tuple;
			var colistup = py_typeof (colkey) == tuple;
			if (rowistup) {
				if (rowkey [1] == null) {
					rowkey [1] = self.ns_nrows;
				}
				else if (rowkey [1] < 0) {
					rowkey [1] += self.ns_nrows;
				}
			}
			if (colistup) {
				if (colkey [1] == null) {
					colkey [1] = self.ns_ncols;
				}
				else if (colkey [1] < 0) {
					colkey [1] += self.ns_ncols;
				}
			}
			if (rowistup || colistup) {
				if (!(rowistup)) {
					var isource = 0;
					if (self.ns_complex) {
						for (var itargetcol of range (...colkey)) {
							var itarget = self.ns_ncols * rowkey + itargetcol;
							self.realbuf [itarget] = value.realbuf [isource];
							self.imagbuf [itarget] = value.imagbuf [isource++];
						}
					}
					else {
						for (var itargetcol of range (...colkey)) {
							result.realbuf [self.ns_ncols * rowkey + itargetcol] = self.realbuf [isource++];
						}
					}
				}
				else if (!(colistup)) {
					var isource = 0;
					if (self.ns_complex) {
						for (var itargetrow of range (...rowkey)) {
							var itarget = self.ns_ncols * itargetrow + colkey;
							self.realbuf [itarget] = value.realbuf [isource];
							self.imagbuf [itarget] = value.imagbuf [isource++];
						}
					}
					else {
						for (var isourcerow of range (...rowkey)) {
							self.realbuf [self.ns_ncols * isourcerow + colkey] = value [isource++];
						}
					}
				}
				else {
					var isource = 0;
					if (self.ns_complex) {
						for (var itargetrow of range (...rowkey)) {
							for (var itargetcol of range (...colkey)) {
								var itarget = self.ns_ncols * itargetrow + itargetcol;
								self.realbuf [itarget] = value.realbuf [isource];
								self.imagbuf [itarget] = value.imagbuf [isource++];
							}
						}
					}
					else {
						for (var isourcerow of range (...rowkey)) {
							for (var isourcecol of range (...colkey)) {
								self.realbuf [self.ns_ncols * itargetrow + itargetcol] = value.realbuf [isource++];
							}
						}
					}
				}
			}
			else if (self.ns_complex) {
				var itarget = self.ns_ncols * key [0] + key [1];
				if (typeof value == 'number') {
					self.realbuf [itarget] = value;
					self.imagbuf [itarget] = 0;
				}
				else {
					self.realbuf [itarget] = value.real;
					self.imagbuf [itarget] = value.imag;
				}
			}
			else {
				self.realbuf [self.ns_ncols * key [0] + key [1]] = value;
			}
		}
	});},
	get real () {return __get__ (this, function (self) {
		return ndarray (self.shape, ns_buffertype (self.dtype), self.realbuf);
	});},
	get imag () {return __get__ (this, function (self) {
		return ndarray (self.shape, ns_buffertype (self.dtype), self.imagbuf);
	});},
	get __conj__ () {return __get__ (this, function (self) {
		if (self.ns_complex) {
			var result = empty (self.shape, self.dtype);
			result.realbuf.set (self.realbuf);
			result.imagbuf = ns_createbuf (true, self.dtype, self.size);
			for (var i = 0; i < self.size; i++) {
				result.imagbuf [i] = -(self.imagbuf [i]);
			}
			return result;
		}
		else {
			return copy (self);
		}
	});},
	get conjugate () {return __get__ (this, function (self) {
		return self.__conj__ ();
	});},
	get __neg__ () {return __get__ (this, function (self) {
		var result = empty (self.shape, self.dtype);
		if (self.ns_complex) {
			for (var i = 0; i < self.size; i++) {
				result.realbuf [i] = -(self.realbuf [i]);
				result.imagbuf [i] = -(self.imagbuf [i]);
			}
		}
		else {
			for (var i = 0; i < self.size; i++) {
				result.realbuf [i] = -(self.realbuf [i]);
			}
		}
		return result;
	});},
	get __ns_inv__ () {return __get__ (this, function (self) {
		var result = empty (self.shape, self.dtype);
		if (self.ns_complex) {
			for (var i = 0; i < self.size; i++) {
				var real = self.realbuf [i];
				var imag = self.imagbuf [i];
				var denom = real * real + imag * imag;
				result.realbuf [i] = real / denom;
				result.imagbuf [i] = -(imag) / denom;
			}
		}
		else {
			for (var i = 0; i < self.size; i++) {
				result.realbuf [i] = 1 / self.realbuf [i];
			}
		}
		return result;
	});},
	get __add__ () {return __get__ (this, function (self, other) {
		var result = empty (self.shape, self.dtype);
		if (py_typeof (other) == ndarray) {
			if (self.ns_complex) {
				for (var i = 0; i < self.size; i++) {
					result.realbuf [i] = self.realbuf [i] + other.realbuf [i];
					result.imagbuf [i] = self.imagbuf [i] + other.imagbuf [i];
				}
			}
			else {
				for (var i = 0; i < self.size; i++) {
					result.realbuf [i] = self.realbuf [i] + other.realbuf [i];
				}
			}
		}
		else if (self.ns_complex) {
			for (var i = 0; i < self.size; i++) {
				result.realbuf [i] = self.realbuf [i] + other.real;
				result.imagbuf [i] = self.imagbuf [i] + other.imag;
			}
		}
		else {
			for (var i = 0; i < self.size; i++) {
				result.realbuf [i] = self.realbuf [i] + other;
			}
		}
		return result;
	});},
	get __radd__ () {return __get__ (this, function (self, scalar) {
		return self.__add__ (scalar);
	});},
	get __sub__ () {return __get__ (this, function (self, other) {
		var result = empty (self.shape, self.dtype);
		if (py_typeof (other) == ndarray) {
			if (self.ns_complex) {
				for (var i = 0; i < self.size; i++) {
					result.realbuf [i] = self.realbuf [i] - other.realbuf [i];
					result.imagbuf [i] = self.imagbuf [i] - other.imagbuf [i];
				}
			}
			else {
				for (var i = 0; i < self.size; i++) {
					result.realbuf [i] = self.realbuf [i] - other.realbuf [i];
				}
			}
		}
		else if (self.ns_complex) {
			for (var i = 0; i < self.size; i++) {
				result.realbuf [i] = self.realbuf [i] - other.real;
				result.imagbuf [i] = self.imagbuf [i] - other.imag;
			}
		}
		else {
			for (var i = 0; i < self.size; i++) {
				result.realbuf [i] = self.realbuf [i] - other;
			}
		}
		return result;
	});},
	get __rsub__ () {return __get__ (this, function (self, scalar) {
		return self.__neg__ ().__add__ (scalar);
	});},
	get __mul__ () {return __get__ (this, function (self, other) {
		var result = empty (self.shape, self.dtype);
		if (py_typeof (other) == ndarray) {
			if (self.ns_complex) {
				for (var i = 0; i < self.size; i++) {
					result.realbuf [i] = self.realbuf [i] * other.realbuf [i] - self.imagbuf [i] * other.imagbuf [i];
					result.imagbuf [i] = self.realbuf [i] * other.imagbuf [i] + self.imagbuf [i] * other.realbuf [i];
				}
			}
			else {
				for (var i = 0; i < self.size; i++) {
					result.realbuf [i] = self.realbuf [i] * other.realbuf [i];
				}
			}
		}
		else if (self.ns_complex) {
			for (var i = 0; i < self.size; i++) {
				result.realbuf [i] = self.realbuf [i] * other.real - self.imagbuf [i] * other.imag;
				result.imagbuf [i] = self.realbuf [i] * other.imag + self.imagbuf [i] * other.real;
			}
		}
		else {
			for (var i = 0; i < self.size; i++) {
				result.realbuf [i] = self.realbuf [i] * other;
			}
		}
		return result;
	});},
	get __rmul__ () {return __get__ (this, function (self, scalar) {
		return self.__mul__ (scalar);
	});},
	get __div__ () {return __get__ (this, function (self, other) {
		var result = empty (self.shape, self.dtype);
		if (py_typeof (other) == ndarray) {
			if (self.ns_complex) {
				for (var i = 0; i < self.size; i++) {
					var real = other.realbuf [i];
					var imag = other.imagbuf [i];
					var denom = real * real + imag * imag;
					result.realbuf [i] = (self.realbuf [i] * real + self.imagbuf [i] * imag) / denom;
					result.imagbuf [i] = (self.imagbuf [i] * real - self.realbuf [i] * imag) / denom;
				}
			}
			else {
				for (var i = 0; i < self.size; i++) {
					result.realbuf [i] = self.realbuf [i] / other.realbuf [i];
				}
			}
		}
		else if (self.ns_complex) {
			var real = other.real;
			var imag = other.imag;
			var denom = real * real + imag * imag;
			for (var i = 0; i < self.size; i++) {
				result.realbuf [i] = (self.realbuf [i] * real + self.imagbuf [i] * imag) / denom;
				result.imagbuf [i] = (self.imagbuf [i] * real - self.realbuf [i] * imag) / denom;
			}
		}
		else {
			for (var i = 0; i < self.size; i++) {
				result.realbuf [i] = self.realbuf [i] / other;
			}
		}
		return result;
	});},
	get __rdiv__ () {return __get__ (this, function (self, scalar) {
		return self.__ns_inv__ ().__mul__ (scalar);
	});},
	get __matmul__ () {return __get__ (this, function (self, other) {
		var result = empty (tuple ([self.ns_nrows, other.ns_ncols]), self.dtype);
		if (self.ns_complex) {
			var iresult = 0;
			for (var irow = 0; irow < self.ns_nrows; irow++) {
				for (var icol = 0; icol < other.ns_ncols; icol++) {
					result.realbuf [iresult] = 0;
					result.imagbuf [iresult] = 0;
					var iself = self.ns_ncols * irow;
					var iother = icol;
					for (var iterm = 0; iterm < self.ns_ncols; iterm++) {
						result.realbuf [iresult] += self.realbuf [iself] * other.realbuf [iother] - self.imagbuf [iself] * other.imagbuf [iother];
						result.imagbuf [iresult] += self.realbuf [iself] * other.imagbuf [iother] + self.imagbuf [iself++] * other.realbuf [iother];
						iother += other.ns_ncols;
					}
					iresult++;
				}
			}
		}
		else {
			var iresult = 0;
			for (var irow = 0; irow < self.ns_nrows; irow++) {
				for (var icol = 0; icol < other.ns_ncols; icol++) {
					result.realbuf [iresult] = 0;
					var iself = self.ns_ncols * irow;
					var iother = icol;
					for (var iterm = 0; iterm < self.ns_ncols; iterm++) {
						result.realbuf [iresult] += self.realbuf [iself++] * other.realbuf [iother];
						iother += other.ns_ncols;
					}
					iresult++;
				}
			}
		}
		return result;
	});}
});
export var empty = function (shape, dtype) {
	if (typeof dtype == 'undefined' || (dtype != null && dtype.hasOwnProperty ("__kwargtrans__"))) {;
		var dtype = 'float64';
	};
	var result = ndarray (shape, dtype);
	result.realbuf = ns_createbuf (false, dtype, result.size);
	result.imagbuf = ns_createbuf (true, dtype, result.size);
	return result;
};
export var array = function (obj, dtype) {
	if (typeof dtype == 'undefined' || (dtype != null && dtype.hasOwnProperty ("__kwargtrans__"))) {;
		var dtype = 'float64';
	};
	if (Array.isArray (obj)) {
		if (len (obj)) {
			if (Array.isArray (obj [0])) {
				var result = empty (tuple ([obj.length, obj [0].length]), dtype);
				var iresult = 0;
				if (result.ns_complex) {
					for (var irow = 0; irow < result.ns_nrows; irow++) {
						for (var icol = 0; icol < result.ns_ncols; icol++) {
							var element = complex (obj [irow] [icol]);
							result.realbuf [iresult] = element.real;
							result.imagbuf [iresult++] = element.imag;
						}
					}
				}
				else {
					for (var irow = 0; irow < result.ns_nrows; irow++) {
						for (var icol = 0; icol < result.ns_ncols; icol++) {
							result.realbuf [iresult++] = obj [irow] [icol];
						}
					}
				}
			}
			else {
				var result = empty (tuple ([obj.length]), dtype);
				if (result.ns_complex) {
					for (var i = 0; i < result.size; i++) {
						var element = complex (obj [i]);
						result.realbuf [i] = element.real;
						result.imagbuf [i] = element.imag;
					}
				}
				else {
					for (var i = 0; i < result.size; i++) {
						result.realbuf [i] = obj [i];
					}
				}
			}
		}
		else {
			var result = empty (tuple ([0]), dtype);
		}
	}
	else {
		var result = empty (obj.shape, obj.dtype);
		result.realbuf.set (obj.realbuf);
		if (obj.ns_complex) {
			result.imagbuf.set (obj.imagbuf);
		}
	}
	return result;
};
export var copy = function (obj) {
	return array (obj, obj.dtype);
};
export var hsplit = function (ary, nparts) {
	var result = (function () {
		var __accu0__ = [];
		for (var ipart = 0; ipart < nparts; ipart++) {
			__accu0__.append (empty (tuple ([ary.ns_nrows, ary.ns_ncols / nparts]), ary.dtype));
		}
		return __accu0__;
	}) ();
	var isource = 0;
	if (ary.ns_complex) {
		for (var irow = 0; irow < ary.ns_nrows; irow++) {
			for (var part of result) {
				var itarget = part.ns_ncols * irow;
				for (var icol = 0; icol < part.ns_ncols; icol++) {
					part.realbuf [itarget] = ary.realbuf [isource];
					part.imagbuf [itarget++] = ary.imagbuf [isource++];
				}
			}
		}
	}
	else {
		for (var irow = 0; irow < ary.ns_nrows; irow++) {
			for (var part of result) {
				var itarget = part.ns_ncols * irow;
				for (var icol = 0; icol < part.ns_ncols; icol++) {
					part.realbuf [itarget++] = ary.realbuf [isource++];
				}
			}
		}
	}
	return result;
};
export var vsplit = function (ary, nparts) {
	var result = (function () {
		var __accu0__ = [];
		for (var ipart = 0; ipart < nparts; ipart++) {
			__accu0__.append (empty (tuple ([ary.ns_nrows / nparts, ary.ns_ncols]), array.dtype));
		}
		return __accu0__;
	}) ();
	var isource = 0;
	if (ary.ns_complex) {
		for (var part of result) {
			for (var itarget = 0; itarget < part.size; itarget++) {
				part.realbuf [itarget] = ary.realbuf [isource];
				part.imagbuf [itarget] = ary.imagbuf [isource++];
			}
		}
	}
	else {
		for (var part of result) {
			for (var itarget = 0; itarget < part.size; itarget++) {
				part.realbuf [itarget] = ary.realbuf [isource++];
			}
		}
	}
	return result;
};
export var hstack = function (tup) {
	var ncols = 0;
	for (var part of tup) {
		ncols += part.ns_ncols;
	}
	var result = empty (tuple ([tup [0].ns_nrows, ncols]), tup [0].dtype);
	var itarget = 0;
	if (result.ns_complex) {
		for (var irow = 0; irow < result.ns_nrows; irow++) {
			for (var part of tup) {
				var isource = part.ns_ncols * irow;
				for (var icol = 0; icol < part.ns_ncols; icol++) {
					result.realbuf [itarget] = part.realbuf [isource];
					result.imagbuf [itarget++] = part.imagbuf [isource++];
				}
			}
		}
	}
	else {
		for (var irow = 0; irow < result.ns_nrows; irow++) {
			for (var part of tup) {
				var isource = part.ns_ncols * irow;
				for (var icol = 0; icol < part.ns_ncols; icol++) {
					result.realbuf [itarget++] = part.realbuf [isource++];
				}
			}
		}
	}
	return result;
};
export var vstack = function (tup) {
	var nrows = 0;
	for (var part of tup) {
		nrows += part.ns_nrows;
	}
	var result = empty (tuple ([nrows, tup [0].ns_ncols]), tup [0].dtype);
	var itarget = 0;
	if (result.ns_complex) {
		for (var part of tup) {
			for (var isource = 0; isource < part.size; isource++) {
				result.realbuf [itarget] = part.realbuf [isource];
				result.imagbuf [itarget++] = part.imagbuf [isource];
			}
		}
	}
	else {
		for (var part of tup) {
			for (var isource = 0; isource < part.size; isource++) {
				result.realbuf [itarget++] = part.realbuf [isource];
			}
		}
	}
	return result;
};
export var round = function (a, decimals) {
	if (typeof decimals == 'undefined' || (decimals != null && decimals.hasOwnProperty ("__kwargtrans__"))) {;
		var decimals = 0;
	};
	var result = empty (a.shape, a.dtype);
	if (a.ns_complex) {
		for (var i = 0; i < a.size; i++) {
			result.realbuf [i] = a.realbuf [i].toFixed (decimals);
			result.imagbuf [i] = a.imagbuf [i].toFixed (decimals);
		}
	}
	else {
		for (var i = 0; i < a.size; i++) {
			result.realbuf [i] = a.realbuf [i].toFixed (decimals);
		}
	}
	return result;
};
export var zeros = function (shape, dtype) {
	if (typeof dtype == 'undefined' || (dtype != null && dtype.hasOwnProperty ("__kwargtrans__"))) {;
		var dtype = 'float64';
	};
	var result = empty (shape, dtype);
	if (result.ns_complex) {
		for (var i = 0; i < result.size; i++) {
			result.realbuf [i] = 0;
			result.imagbuf [i] = 0;
		}
	}
	else {
		for (var i = 0; i < result.size; i++) {
			result.realbuf [i] = 0;
		}
	}
	return result;
};
export var ones = function (shape, dtype) {
	if (typeof dtype == 'undefined' || (dtype != null && dtype.hasOwnProperty ("__kwargtrans__"))) {;
		var dtype = 'float64';
	};
	var result = empty (shape, dtype);
	if (result.ns_complex) {
		for (var i = 0; i < result.size; i++) {
			result.realbuf [i] = 1;
			result.imagbuf [i] = 0;
		}
	}
	else {
		for (var i = 0; i < result.size; i++) {
			result.realbuf [i] = 1;
		}
	}
	return result;
};
export var identity = function (n, dtype) {
	if (typeof dtype == 'undefined' || (dtype != null && dtype.hasOwnProperty ("__kwargtrans__"))) {;
		var dtype = 'float64';
	};
	var result = zeros (tuple ([n, n]), dtype);
	var i = 0;
	var shift = n + 1;
	for (var j = 0; j < n; j++) {
		result.realbuf [i] = 1;
		i += shift;
	}
	return result;
};
export var conjugate = function (x) {
	return x.__conj__ ();
};

//# sourceMappingURL=numscrypt.map