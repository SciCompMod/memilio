#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Rene Schmieding
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################
import sys
import time
import threading
from os import get_terminal_size, name as os_name
from warnings import warn

class ProgressIndicator:
    """! Print an animation to show that something is happening.

    Animations are rendered in a new thread, which is set up as deamon so that
    it stops on main thread exit.
    The methods Dots, Spinner and Percentage provide some default animations.
    Start the animation by either using the `start`/`stop` functions or a
    'with as' block, e.g
    ```
    import memilio.ProgressIndicator
    with ProgressIndicator.Spinner():
        do_something()
    ```
    or
    ```
    with ProgressIndicator.Percentage() as indicator:
        for i in range(n) :
            do_something()
            indicator.set_progress((i+1)/n)
    ```
    """

    _first_init = True

    def __init__(self, animator, delay):
        """! Create an ProgressIndicator.

        @param animator generator expression. The expression must loop, i.e.
            it should never return StopIteration. `next(animator)` must be a
            string of length < os.get_terminal_size().columns
        @param delay positive real number. Sets delay in seconds between
            drawing animation frames.
        """
        assert(delay > 0)
        self._animator = animator
        self._delay = delay
        self._thread = None
        self._enabled = False
        if ProgressIndicator._first_init:
            ProgressIndicator._first_init = False
            ProgressIndicator._console_setup()

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exception, value, trace):
        self.stop()

    @staticmethod
    def _console_setup():
        if os_name == 'nt': # os name can be nt, posix, or java
            # Windows uses nt, which does not support carriage returns by
            # default. the following Windows specific module should fix this.
            try:
                from ctypes import windll
                k = windll.kernel32
                # GetStdHandle(-11) is the standart output device, the flag 7
                # is equal to
                #    ENABLE_PROCESSED_OUTPUT | ENABLE_WRAP_AT_EOL_OUTPUT
                #        | ENABLE_VIRTUAL_TERMINAL_PROCESSING
                # the last flag enables several control sequences, like \r
                k.SetConsoleMode(k.GetStdHandle(-11), 7)
            except (ImportError):
                msg = "Failed to set console mode for 'nt' system (e.g."\
                      " Windows). ProgressIndicator(s) may be displayed"\
                      " incorrectly." 
                warn(msg, category=RuntimeWarning, stacklevel=4)

    def _render(self):
        """! Regularly update the animation. Do not call manually!"""
        while self._enabled:
            self.show()
            # wait before writing next frame
            time.sleep(self._delay)

    def show(self):
        """! Print and advance the animation. """
        # write animation and proceed to next frame
        sys.stdout.write(" {}\r".format(next(self._animator)))
        sys.stdout.flush()

    def start(self):
        """! Start the animation in a new thread.

        Must call stop() afterwards.
        """
        if not self._enabled:
            self._enabled = True
            # start new threat to render the animator in the background
            self._thread = threading.Thread(target=self._render)
            self._thread.setDaemon(True) # stops thread on main thread exit
            self._thread.start()

    def stop(self):
        """! Stop the animation and join the thread. """
        if self._enabled:
            self._enabled = False
            if self._thread and self._thread.is_alive():
                self._thread.join()
            sys.stdout.write("\033[K") # clear line

class Spinner(ProgressIndicator):
    """! Subclass of ProgressIndicator with a predefined animation. """
    def __init__(self, delay=0.1, message=""):
        """! initializes a ProgressIndicator with a rotating line animation.

        This method spawns a new thread to print the animation.
        Start the animation by either using the `start`/`stop` functions or a
        'with' block.

        @param delay [Default: 0.1] positive real number. Sets delay in
            seconds between drawing animation frames.
        @param message [Default: ""] string. Text shown before the indicator
            (consider appeding a space as separator).
        """
        # prevent spamming output with messages longer than a single line
        if get_terminal_size().columns < len(message) + 2:
            sys.stdout.write(message)
            message = ""

        def _spin():
            while True: # loop animation
                for s in "|/-\\": # iterate animation frames
                    yield "{}{}".format(message, s) # return single frame
        super().__init__(_spin(), delay)

class Dots(ProgressIndicator):
    """! Subclass of ProgressIndicator with a predefined animation. """
    def __init__(self, delay=1, message="", num_dots=3, dot=".", blank=" "):
        """! initializes ProgressIndicator with a 'dot, dot, dot' animation.

        This method spawns a new thread to print the animation.
        Start the animation by either using the `start`/`stop` functions or a
        'with' block.

        @param delay [Default: 1] positive real number. Sets delay in seconds
            between drawing animation frames.
        @param message [Default: ""] string. Text shown before the indicator
            (consider appeding a space as separator).
        @param num_dots [Default: 3] positive integer. Determines maximum
            number of dots drawn before clearing them with blanks.
        @param dot [Default: "."] string. Drawn up to num_dots times.
        @param blank [Default: " "] string. Placeholder for yet to be drawn
            dots. Must have same length as dot.
        """
        assert(len(dot) == len(blank))
        assert(num_dots > 0)
        # prevent spamming output with messages longer than a single line
        if get_terminal_size().columns < num_dots + len(message) + 1:
            sys.stdout.write(message)
            message = ""

        def _dots():
            while True: # loop animation
                for n in range(1, num_dots+1): # iterate animation frames
                    # return single frame
                    yield "{}{}{}".format(message, dot*n, blank*(num_dots-n))
        super().__init__(_dots(), delay)

class Percentage:
    """! Manages a ProgressIndicator with a predefined animation. """
    def __init__(self, delay=1, message="", percentage=0, use_bar=False,
            use_delayed_output=True, keep_output=True):
        """! initializes ProgressIndicator displaying a percentage.

        The percentage can be updated using the `set_progress` method.
        By default, this method spawns a new thread to print the animation.
        If use_delayed_output is set to False, the delay is ignored, and no
        new thread is spawned. The output is then updated in the main thread,
        whenever 'set_progress' is called.
        Start the animation by either using the `start`/`stop` functions or a
        `with` `as` block.

        @param delay [Default: 1] positive real number. Sets delay in seconds
            between drawing animation frames.
        @param message [Default: ""] string. Text shown before the indicator
            (consider appeding a space as separator).
        @param percentage [Default: 0] real number in [0, 1]. Initial
            percentage shown in the animation.
        @param use_bar [Default: False] bool. If True, adds a bar plotting
            the current progress.
        @param use_delayed_output [Default: True] bool. If False, delay is
            ignored and the animation is drawn in the main thread whenever
            set_progress() is called.
        @param keep_output [Default: True] bool. Set F if the last animation
            frame should be kept as a new line.
        """
        self._use_thread = use_delayed_output
        self._keep_output = keep_output
        # use lists to pass variables by reference
        self._progress = [percentage]
        width = get_terminal_size().columns
        if use_bar:
            # prevent spamming output with messages longer than a single line
            if width < len(message) + 12:
                sys.stdout.write(message)
                message = ""
            self._message = self._bar(self._progress, message, width)
        else:
            # prevent spamming output with messages longer than a single line
            if width < len(message) + 8:
                sys.stdout.write(message)
                message = ""
            self._message = [message]
        animator = self._perc(self._progress, self._message)
        self._indicator = ProgressIndicator(animator, delay)

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exception, value, trace):
        self.stop()

    def start(self):
        """! Start the animation. Must call stop() afterwards.
        
        By default, this method spawns a new thread. See option
        use_delayed_output in the class constructor.
        """
        if self._use_thread:
            self._indicator.start()

    def stop(self):
        """! Stops the animation. """
        if self._use_thread:
            self._indicator.show() # print manually to catch current progress
        if self._keep_output:
            sys.stdout.write("\n") # newline to keep output
        if self._use_thread:
            self._indicator.stop()
        else:
            sys.stdout.write("\033[K") # clear line

    def set_progress(self, percentage):
        """! Updates the percentage shown by the indicator.

        @param percentage real number. Must be in the interval [0, 1].
        """
        self._progress[0] = percentage
        if not self._indicator._enabled:
            self._indicator.show()

    @staticmethod
    def _perc(p, message):
        """! Method to create internal animator function.
        
        Uses lists to get parameters by reference.
        """
        while True:
            yield "{}{:6.2f}%".format(message[0], 100*p[0])

    class _bar:
        """! Uses [] to return a progress bar string if use_bar is enabled."""
        def __init__(self, percentage, message, width):
            self.p = percentage
            self.m = message
            # collected offset by percentage/extra characters from bar = 11
            self.w = width - len(message) - 11
        def __getitem__(self, _): # ignore key
            n = int(self.w * self.p[0])
            return self.m + "[" + "#" * n + " " * (self.w - n) + "] "

if __name__ == "__main__":
    print("This is only a usage example, and does not actually do anything.")
    # using start/stop
    p = Dots(message="waiting", delay=0.5)
    p.start()
    time.sleep(1.6)
    p.stop()
    # using with as block
    with Percentage(message="download 1 ", use_bar=True, delay=0.4) as p:
        for i in range(13):
            time.sleep(0.1467)
            p.set_progress((i+1)/13)
    with Percentage(message="download 2 ", use_delayed_output=False,
            keep_output=False) as p:
        for i in range(97):
            time.sleep(0.0367)
            p.set_progress((i+1)/97)
    # using with block ('as' is not usefull without Percentage)
    with Spinner(message="finish "):
        time.sleep(2)