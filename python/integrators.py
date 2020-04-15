transpile = True

if transpile:
    from org.transcrypt.stubs.browser import *
    from org.transcrypt.stubs.browser import __main__, __envir__, __pragma__

# Imports for Transcrypt, skipped runtime by CPython
    if __envir__.executor_name == __envir__.transpiler_name:
        import numscrypt as np

else:

    # Imports for CPython, skipped compile time by Transcrypt
    #__pragma__ ('skip')
    import numpy as np
    from scipy.integrate import odeint
    #__pragma__ ('noskip')

class Integrator:
    def __init__(self, dfun, xzero, timerange):
        '''Initializes dfun, xzero, timestart, timeend,
        sets time at start and steps to zero.'''
        assert len(timerange) == 2
        self.timestart, self.timeend = timerange
        self.time = self.timestart
        # check output of dfun, wrap w/ np.array if necessary
        xtest = dfun(self.time, xzero)
        if not isinstance(xtest, np.ndarray):
            xtest = np.array(xtest)
            def array_dfun(t, x):
                x = dfun(t, x)
                xarray = np.array(x)
                return xarray
            self.dfun = array_dfun
        else:
            self.dfun = dfun
        if len(xtest.shape) != 1:
            raise Exception(f'dfun: {dfun} output is not one dimensional')
        if not isinstance(xzero, np.ndarray):
            xzero = np.array(xzero)
        assert len(xzero.shape) == 1, 'xzero must be one dimensional'
        self.x = xzero
        self.stepcounter = 0

    def __iter__(self):
        return self

class ConstantTimestep(Integrator):
    '''The __init__ function of this class sets instance variables for
    integrators with a constant timestep.'''
    def __init__(self, dfun, xzero, timerange, timestep):
        super().__init__(dfun, xzero, timerange)
        assert ((self.timeend - self.timestart) / timestep) > 0, (
            "'timerange' and 'timestep' not consistant. "
            "Check signs and order.")
        self.timestep = timestep
        self.direction = np.sign(timestep)
        self.steps = np.ceil((self.timeend - self.timestart) / timestep)
        self.status = 'initialized'


class Euler(ConstantTimestep):
    '''Euler method integration. This class implements a generator.

        :param dfun:
            derivative function of the system.
            The differential system arranged as a series of first-order
            equations: :math:`\dot{X} = \mathrm{dfun}(t, x)`.
            Returns :math:`\dot{X}` should be a single dimensional array
            or list.
        :param xzero:
            the initial condition of the system
        :param timerange:
            the start and end times as (starttime, endtime) tuple/list/array.
        :param timestep:
            the timestep
        :returns: t, x for each iteration. t is a number. x is an array.
    '''

    def __next__(self):
        if self.stepcounter < self.steps:
            if self.status == 'initialized':
                self.status = 'running'
                return self.time, self.x
            else:
                self.stepcounter += 1
                dx = self.dfun(self.time, self.x)
                self.time, self.x = (
                    self.timestart + (self.stepcounter * self.timestep),
                    self.x + (self.timestep * dx))
                return self.time, self.x
        else:
            self.status = 'finished'
            raise StopIteration


def euler(dfun, xzero, timerange, timestep):
    '''Euler method integration. This function wraps the Euler class.

        :param All: All parameters are identical to the Euler class above.
        :returns: t, x as arrays.
    '''
    t_column, X = zip(*list(Euler(dfun, xzero, timerange, timestep)))
    t_column = np.array(t_column)
    X_columns = np.vstack(X).T
    return t_column, X_columns