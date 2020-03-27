// Transcrypt'ed from Python, 2020-03-27 23:22:06
import {AssertionError, AttributeError, BaseException, DeprecationWarning, Exception, IndexError, IterableError, KeyError, NotImplementedError, RuntimeWarning, StopIteration, UserWarning, ValueError, Warning, __JsIterator__, __PyIterator__, __Terminal__, __add__, __and__, __call__, __class__, __conj__, __envir__, __eq__, __floordiv__, __ge__, __get__, __getcm__, __getitem__, __getslice__, __getsm__, __gt__, __i__, __iadd__, __iand__, __idiv__, __ijsmod__, __ilshift__, __imatmul__, __imod__, __imul__, __in__, __init__, __ior__, __ipow__, __irshift__, __isub__, __ixor__, __jsUsePyNext__, __jsmod__, __k__, __kwargtrans__, __le__, __lshift__, __lt__, __matmul__, __mergefields__, __mergekwargtrans__, __mod__, __mul__, __ne__, __neg__, __nest__, __or__, __pow__, __pragma__, __proxy__, __pyUseJsNext__, __rshift__, __setitem__, __setproperty__, __setslice__, __sort__, __specialattrib__, __sub__, __super__, __t__, __terminal__, __truediv__, __withblock__, __xor__, abs, all, any, assert, bool, bytearray, bytes, callable, chr, complex, copy, deepcopy, delattr, dict, dir, divmod, enumerate, filter, float, getattr, hasattr, input, int, isinstance, issubclass, len, list, map, max, min, object, ord, pow, print, property, py_TypeError, py_iter, py_metatype, py_next, py_reversed, py_typeof, range, repr, round, set, setattr, sorted, str, sum, tuple, zip} from './org.transcrypt.__runtime__.js';
import * as np from './numscrypt.js';
var __name__ = 'integrators';
export var transpile = true;
if (transpile) {
	if (__envir__.executor_name == __envir__.transpiler_name) {
	}
}
else {
}
export var Integrator =  __class__ ('Integrator', [object], {
	__module__: __name__,
	get __init__ () {return __get__ (this, function (self, dfun, xzero, timerange) {
		var __left0__ = timerange;
		self.timestart = __left0__ [0];
		self.timeend = __left0__ [1];
		self.time = self.timestart;
		var xtest = dfun (self.time, xzero);
		if (!(isinstance (xtest, np.ndarray))) {
			var xtest = np.array (xtest);
			var array_dfun = function (t, x) {
				var x = dfun (t, x);
				var xarray = np.array (x);
				return xarray;
			};
			self.dfun = array_dfun;
		}
		else {
			self.dfun = dfun;
		}
		if (len (xtest.shape) != 1) {
			var __except0__ = Exception ('dfun: {} output is not one dimensional'.format (dfun));
			__except0__.__cause__ = null;
			throw __except0__;
		}
		if (!(isinstance (xzero, np.ndarray))) {
			var xzero = np.array (xzero);
		}
		self.x = xzero;
		self.stepcounter = 0;
	});},
	get __iter__ () {return __get__ (this, function (self) {
		return self;
	});},
	[Symbol.iterator] () {return this.__iter__ ()}
});
export var ConstantTimestep =  __class__ ('ConstantTimestep', [Integrator], {
	__module__: __name__,
	get __init__ () {return __get__ (this, function (self, dfun, xzero, timerange, timestep) {
		__super__ (ConstantTimestep, '__init__') (self, dfun, xzero, timerange);
		self.timestep = timestep;
		self.direction = np.sign (timestep);
		self.steps = np.ceil ((self.timeend - self.timestart) / timestep);
		self.status = 'initialized';
	});}
});
export var Euler =  __class__ ('Euler', [ConstantTimestep], {
	__module__: __name__,
	get __next__ () {return __get__ (this, function (self) {
		if (self.stepcounter < self.steps) {
			if (self.status == 'initialized') {
				self.status = 'running';
				return tuple ([self.time, self.x]);
			}
			else {
				self.stepcounter++;
				var dx = self.dfun (self.time, self.x);
				var __left0__ = tuple ([self.timestart + self.stepcounter * self.timestep, self.x + self.timestep * dx]);
				self.time = __left0__ [0];
				self.x = __left0__ [1];
				return tuple ([self.time, self.x]);
			}
		}
		else {
			self.status = 'finished';
			var __except0__ = StopIteration;
			__except0__.__cause__ = null;
			throw __except0__;
		}
	});},
	next: __jsUsePyNext__
});
export var euler = function (dfun, xzero, timerange, timestep) {
	var __left0__ = zip (...list (Euler (dfun, xzero, timerange, timestep)));
	var t_column = __left0__ [0];
	var X = __left0__ [1];
	var t_column = np.array (t_column);
	var X_columns = np.vstack (X).T;
	return tuple ([t_column, X_columns]);
};

//# sourceMappingURL=integrators.map