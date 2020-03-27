// Transcrypt'ed from Python, 2020-03-27 23:22:05
import {AssertionError, AttributeError, BaseException, DeprecationWarning, Exception, IndexError, IterableError, KeyError, NotImplementedError, RuntimeWarning, StopIteration, UserWarning, ValueError, Warning, __JsIterator__, __PyIterator__, __Terminal__, __add__, __and__, __call__, __class__, __conj__, __envir__, __eq__, __floordiv__, __ge__, __get__, __getcm__, __getitem__, __getslice__, __getsm__, __gt__, __i__, __iadd__, __iand__, __idiv__, __ijsmod__, __ilshift__, __imatmul__, __imod__, __imul__, __in__, __init__, __ior__, __ipow__, __irshift__, __isub__, __ixor__, __jsUsePyNext__, __jsmod__, __k__, __kwargtrans__, __le__, __lshift__, __lt__, __matmul__, __mergefields__, __mergekwargtrans__, __mod__, __mul__, __ne__, __neg__, __nest__, __or__, __pow__, __pragma__, __proxy__, __pyUseJsNext__, __rshift__, __setitem__, __setproperty__, __setslice__, __sort__, __specialattrib__, __sub__, __super__, __t__, __terminal__, __truediv__, __withblock__, __xor__, abs, all, any, assert, bool, bytearray, bytes, callable, chr, complex, copy, deepcopy, delattr, dict, dir, divmod, enumerate, filter, float, getattr, hasattr, input, int, isinstance, issubclass, len, list, map, max, min, object, ord, pow, print, property, py_TypeError, py_iter, py_metatype, py_next, py_reversed, py_typeof, range, repr, round, set, setattr, sorted, str, sum, tuple, zip} from './org.transcrypt.__runtime__.js';
import * as np from './numscrypt.js';
import {euler} from './integrators.js';
var __name__ = '__main__';
export var transpile = true;
if (transpile) {
	if (__envir__.executor_name == __envir__.transpiler_name) {
	}
}
else {
}
export var beta = 1.75;
export var deriv = function (y, t, N, alpha, beta, gamma) {
	var __left0__ = y;
	var S = __left0__ [0];
	var E = __left0__ [1];
	var I = __left0__ [2];
	var R = __left0__ [3];
	var dSdt = ((-(beta) * S) * I) / N;
	var dEdt = ((beta * S) * I) / N - alpha * E;
	var dIdt = alpha * E - gamma * I;
	var dRdt = gamma * I;
	return tuple ([dSdt, dEdt, dIdt, dRdt]);
};
export var deriv_social = function (y, t, N, alpha, beta, gamma, rho) {
	var __left0__ = y;
	var S = __left0__ [0];
	var E = __left0__ [1];
	var I = __left0__ [2];
	var R = __left0__ [3];
	var dSdt = (((-(rho) * beta) * S) * I) / N;
	var dEdt = (((rho * beta) * S) * I) / N - alpha * E;
	var dIdt = alpha * E - gamma * I;
	var dRdt = gamma * I;
	return tuple ([dSdt, dEdt, dIdt, dRdt]);
};
export var my_odeint = function (func, y0, t, args) {
	if (typeof args == 'undefined' || (args != null && args.hasOwnProperty ("__kwargtrans__"))) {;
		var args = tuple ([]);
	};
	var t_min = t [0];
	var t_max = t [-(1)];
	var dt = t [1] - t [0];
	var fun_helper = function (t, y) {
		return func (y, t, ...args);
	};
	var __left0__ = euler (fun_helper, y0, tuple ([t_min, t_max]), dt);
	var tp = __left0__ [0];
	var vals = __left0__ [1];
	return vals.__getitem__ ([tuple ([0, null, 1]), tuple ([0, t.shape [0], 1])]).T;
};
export var seir_solve = function (N, S0, E0, I0, R0, t_max, dt, alpha, gamma, rho, t_min) {
	if (typeof t_min == 'undefined' || (t_min != null && t_min.hasOwnProperty ("__kwargtrans__"))) {;
		var t_min = 0;
	};
	var y0 = tuple ([S0, E0, I0, R0]);
	var t = np.linspace (t_min, t_max, int ((t_max - t_min) / dt));
	var ret = my_odeint (deriv_social, y0, t, __kwargtrans__ ({args: tuple ([N, alpha, beta, gamma, rho])}));
	var __left0__ = ret.T;
	var S = __left0__ [0];
	var E = __left0__ [1];
	var I = __left0__ [2];
	var R = __left0__ [3];
	return tuple ([S, E, I, R]);
};
export var EpiParms =  __class__ ('EpiParms', [object], {
	__module__: __name__,
	N: 10000,
	t_max: 200,
	dt: 0.1,
	t_incubation: 5.2,
	t_infectious: 2,
	rho: 0.5
});
var __left0__ = tuple ([1, 0, 0]);
EpiParms.E0 = __left0__ [0];
EpiParms.I0 = __left0__ [1];
EpiParms.R0 = __left0__ [2];
export var epi = function (parameters) {
	var N = parameters.N;
	var E0 = parameters.E0;
	var I0 = parameters.I0;
	var R0 = parameters.R0;
	var t_max = parameters.t_max;
	var dt = parameters.dt;
	var rho = parameters.rho;
	var S0 = ((parameters.N - parameters.I0) - parameters.R0) - parameters.E0;
	var alpha = 1.0 / parameters.t_incubation;
	var gamma = 1.0 / parameters.t_infectious;
	var __left0__ = seir_solve (N, S0, E0, I0, R0, t_max, dt, alpha, gamma, rho, __kwargtrans__ ({t_min: 0}));
	var SA = __left0__ [0];
	var EA = __left0__ [1];
	var IA = __left0__ [2];
	var RA = __left0__ [3];
	var day_party = 100;
	var N_party = 50;
	var rho_party = 1.0;
	var length_party = 0.5;
	var __left0__ = seir_solve (N, S0, E0, I0, R0, day_party, dt, alpha, gamma, rho, __kwargtrans__ ({t_min: 0}));
	var S_till_party = __left0__ [0];
	var E_till_party = __left0__ [1];
	var I_till_party = __left0__ [2];
	var R_till_party = __left0__ [3];
	var party_factor = N_party / N;
	var IP0 = I_till_party [-(1)] * party_factor;
	var EP0 = E_till_party [-(1)] * party_factor;
	var RP0 = R_till_party [-(1)] * party_factor;
	var SP0 = ((N_party - IP0) - RP0) - EP0;
	var __left0__ = seir_solve (N_party, SP0, EP0, IP0, RP0, day_party + length_party, dt, alpha, gamma, rho_party, __kwargtrans__ ({t_min: day_party}));
	var S_party = __left0__ [0];
	var E_party = __left0__ [1];
	var I_party = __left0__ [2];
	var R_party = __left0__ [3];
	var rest_factor = (N - N_party) / N;
	var I0 = I_till_party [-(1)] * rest_factor;
	var E0 = E_till_party [-(1)] * rest_factor;
	var R0 = R_till_party [-(1)] * rest_factor;
	var N_rest = N - N_party;
	var S0 = ((N_rest - I0) - R0) - E0;
	var __left0__ = seir_solve (N_rest, S0, E0, I0, R0, day_party + length_party, dt, alpha, gamma, rho, __kwargtrans__ ({t_min: day_party}));
	var SR = __left0__ [0];
	var ER = __left0__ [1];
	var IR = __left0__ [2];
	var RR = __left0__ [3];
	var E0 = ER [-(1)] + E_party [-(1)];
	var I0 = IR [-(1)] + I_party [-(1)];
	var R0 = RR [-(1)] + R_party [-(1)];
	var S0 = ((N - I0) - R0) - E0;
	var __left0__ = seir_solve (N, S0, E0, I0, R0, t_max, dt, alpha, gamma, rho, __kwargtrans__ ({t_min: day_party + length_party}));
	var SL = __left0__ [0];
	var EL = __left0__ [1];
	var IL = __left0__ [2];
	var RL = __left0__ [3];
	var t = np.linspace (0, t_max, int (t_max / dt));
	var S_final = np.concatenate (tuple ([S_till_party, S_party + SR, SL]), __kwargtrans__ ({axis: 0}));
	var E_final = np.concatenate (tuple ([E_till_party, E_party + ER, EL]), __kwargtrans__ ({axis: 0}));
	var I_final = np.concatenate (tuple ([I_till_party, I_party + IR, IL]), __kwargtrans__ ({axis: 0}));
	var R_final = np.concatenate (tuple ([R_till_party, R_party + RR, RL]), __kwargtrans__ ({axis: 0}));
	var E_party_max = np.amax (E_final);
	var I_party_max = np.amax (I_final);
	var R_party_max = np.amax (R_final);
	var EA_max = np.amax (EA);
	var IA_max = np.amax (IA);
	var RA_max = np.amax (RA);
	var r1 = R_final [-(1)] - RA [-(1)];
	print ('Population:', N, 'Day of party:', day_party, 'Party people:', N_party, 'Death rate assumption: 2%');
	print ('At the worst day in our scenario, ', int (E_party_max + I_party_max), ' people are infected.');
	print ('A the end of the crisis', int (R_party_max), 'people have been infected,');
	print ('without the party it would have been', round (r1, 3), 'less.');
	print ('You think that this does not make a difference?');
	var factor = 83783942 / N;
	print ('In all of Germany this makes a difference of', int (r1 * factor), 'people who are infected.');
	print (int (((r1 * factor) * 2) / 100), 'people will have died due to the parties.');
	return tuple ([S_final, E_final, I_final, R_final, t]);
};

//# sourceMappingURL=epi.map