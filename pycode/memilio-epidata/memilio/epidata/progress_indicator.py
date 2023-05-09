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
from shutil import get_terminal_size
from warnings import warn
from os import name as os_name


class ProgressIndicator:
    """! Print an animation to show that something is happening.

    Animations are rendered in a new thread, which is set up as deamon so that
    it stops on main thread exit.
    The methods Dots, Spinner and Percentage provide some default animations.
    Start the animation by either using the `start`/`stop` functions or a
    'with as' block, e.g
    ```
    with ProgressIndicator.Spinner():
        do_something...
    ```
    or
    ```
    with ProgressIndicator.Percentage() as indicator:
        for i in range(n):
            do_something...
            indicator.set_progress((i+1)/n)
    ```
    """

    _first_init = True
    _disabled = False

    def __init__(self, message, animator, delay, auto_adjust=False):
        """! Create a ProgressIndicator.

        @param message String. Shown in front of the animation without
            seperator. If it would not fit in a single line with the
            animation, it will be printed once in a new line above the
            animation instead. Must not contain `\\r` or `\\n`.
        @param animator Generator expression. The expression must loop, i.e.
            it should never return StopIteration. `next(animator)` must be a
            string of length < `os.get_terminal_size().columns`. The length
            should be constant, otherwise animations may leave artifacts (this
            can be worked around by prepending `"\\033[K"` to each string).
        @param delay Positive real number. Sets delay in seconds between
            drawing animation frames.
        @param auto_adjust [Default: False] Bool. Specify whether each frame
            of the animation should be forced to fit in a single line. Can be
            usefull for long, line filling animations.
        """
        assert (delay > 0)
        self._message = message
        self._animator = animator
        self._delay = delay
        self._thread = None
        self._thread_is_running = False  # whether _thread is running
        self._auto_adjust = auto_adjust
        self._space = 0  # free space in a line, set by _adjust_to_terminal_size
        self._frame = ""

        self._adjust_to_terminal_size()
        if ProgressIndicator._first_init:
            ProgressIndicator._first_init = False
            ProgressIndicator._console_setup()

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exception, value, trace):
        self.stop()

    @staticmethod
    def disable_indicators(disable):
        """! Globally prevents new indicators from starting.

        This does not affect currently running indicators.

        @param disable Bool. If True, no new indicators can be started.
            If False, resume default behaviour.
        """
        ProgressIndicator._disabled = disable

    @staticmethod
    def _console_setup():
        if os_name == 'nt':  # os name can be nt, posix, or java
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
                warn(msg, category=RuntimeWarning, stacklevel=2)

    def _adjust_to_terminal_size(self, reserve=0):
        """! Keep animation in a single line."""
        # use width - 1, since show() prepends a space
        width = get_terminal_size().columns - 1
        self._frame = self._frame[:(width)]
        if len(self._message) + len(self._frame) >= width - reserve:
            sys.stdout.write(self._message + "\n")
            self._message = ""
        self._space = width - len(self._frame) - len(self._message)

    def _advance(self):
        """! Advance the animation to the next frame. """
        self._frame = next(self._animator)
        if self._auto_adjust:
            self._adjust_to_terminal_size()

    def _render(self):
        """! Regularly update the animation. Do not call manually! """
        while self._thread_is_running:
            self.step()
            # wait before writing next frame
            time.sleep(self._delay)

    def set_message(self, message):
        """! Change the message displayed in front of the animation.
        @param message String. Shown in front of the animation without
            seperator. If it would not fit in a single line with the
            animation, it will be printed once in a new line above the
            animation instead. Must not contain `\\r` or `\\n`.
        """
        self._message = message
        self._adjust_to_terminal_size()

    def show(self):
        """! Print the animation without advancing it."""
        # prepend a space to fit the cursor
        sys.stdout.write(f" {self._message}{self._frame}\r")
        sys.stdout.flush()

    def step(self):
        """! Advance and print the animation. """
        self._advance()
        self.show()

    def start(self):
        """! Start the animation in a new thread.

        Must call stop() afterwards.
        """
        if not ProgressIndicator._disabled and not self._thread_is_running:
            self._thread_is_running = True
            # start new threat to render the animator in the background
            self._thread = threading.Thread(target=self._render)
            self._thread.daemon = True  # stops thread on main thread exit
            self._thread.start()

    def stop(self):
        """! Stop the animation and join the thread. """
        if self._thread_is_running:
            self._thread_is_running = False
            if self._thread and self._thread.is_alive():
                self._thread.join()
            sys.stdout.write("\033[K")  # clear line


class Spinner(ProgressIndicator):
    """! Subclass of ProgressIndicator with a predefined animation. """

    def __init__(self, delay=0.1, message=""):
        """! initializes a ProgressIndicator with a rotating line animation.

        This method spawns a new thread to print the animation.
        Start the animation by either using the `start`/`stop` functions or a
        'with' block.

        @param delay [Default: 0.1] positive real number. Sets delay in
            seconds between drawing animation frames.
        @param message [Default: ""] string. Text shown before the indicator.
        """
        def _spin():
            while True:  # loop animation
                yield from "|/-\\"
        super().__init__(message + " ", _spin(), delay)


class Dots(ProgressIndicator):
    """! Subclass of ProgressIndicator with a predefined animation. """

    def __init__(self, delay=1, message="", num_dots=3, dot=".", blank=" "):
        """! initializes ProgressIndicator with a 'dot, dot, dot' animation.

        This method spawns a new thread to print the animation.
        Start the animation by either using the `start`/`stop` functions or a
        'with' block.

        @param delay [Default: 1] positive real number. Sets delay in seconds
            between drawing animation frames.
        @param message [Default: ""] string. Text shown before the indicator.
        @param num_dots [Default: 3] positive integer. Determines maximum
            number of dots drawn before clearing them with blanks.
        @param dot [Default: "."] string. Drawn up to num_dots times.
        @param blank [Default: " "] string. Placeholder for yet to be drawn
            dots. Must have same length as dot.
        """
        assert (len(dot) == len(blank))
        assert (num_dots > 0)

        def _dots():
            while True:  # loop animation
                for n in range(1, num_dots+1):  # iterate animation frames
                    # return single frame
                    yield f"{dot*n}{blank*(num_dots-n)}"
        super().__init__(message + " ", _dots(), delay)


class Percentage(ProgressIndicator):
    """! Manages a ProgressIndicator with a predefined animation. """

    def __init__(self, delay=1, message="", percentage=0, use_bar=True,
                 keep_output=False):
        """! initializes ProgressIndicator displaying a percentage.

        The percentage can be updated using the `set_progress` method.
        By default, this method spawns a new thread to print the animation.
        Start the animation by either using the `start`/`stop` functions or a
        `with` `as` block.

        @param delay [Default: 1] non-negative real number. Sets delay in
            seconds between drawing animation frames. If delay is set to 0, no
            new thread is spawned. The output is then updated in the main
            thread, whenever 'set_progress' is called.
        @param message [Default: ""] string. Text shown before the indicator.
        @param percentage [Default: 0] real number in [0, 1]. Initial
            percentage shown in the animation.
        @param use_bar [Default: True] bool. If True, adds a progress bar
            visualizing the current percentage.
        @param keep_output [Default: False] bool. Whether the last animation
            frame should be kept as a new line  when stopping.
        """
        if delay == 0:
            self._use_thread = False
            delay = 1  # arbitrary, will not be used outside of init
        else:
            self._use_thread = True
        self._keep_output = keep_output
        self._use_bar = use_bar
        self._progress = percentage
        self._disabled = False

        def _perc():
            while True:
                yield f"{100*self._progress:6.2f}%"
        super().__init__(message + " ", _perc(), delay, use_bar)

    def _advance(self):
        """! Advance the animation to the next frame. """
        self._frame = next(self._animator)
        if self._use_bar:
            self._adjust_to_terminal_size(5)
            # prepend bar to frame
            self._frame = self._bar(self._space, self._progress) + self._frame

    def start(self):
        """! Start the animation. Must call stop() afterwards.

        If delay > 0, this method spawns a new thread.
        """
        self._disabled = ProgressIndicator._disabled
        if self._use_thread:
            super().start()

    def stop(self):
        """! Stops the animation. """
        if not self._disabled:
            if self._use_thread:
                self.step()
            if self._keep_output:
                sys.stdout.write("\n")  # newline to 'save' output
            else:
                sys.stdout.write("\033[K")  # clear line
        if self._use_thread:
            super().stop()

    def set_progress(self, percentage):
        """! Updates the percentage shown by the indicator. Steps if delay = 0.

        @param percentage real number. Must be in the interval [0, 1].
        """
        self._progress = percentage
        if not self._disabled and not self._thread_is_running:
            self.step()

    @staticmethod
    def _bar(width, percentage):
        """! Returns a progress bar.

        @param width Total width of the bar. Must be at least 4.
        @param percentage Float in [0,1].
        @return String of length width, visualizing percentage.
        """
        w = width - 3  # 3 == len("[] ")
        n = int(w * percentage)
        return "[" + "#" * n + " " * (w - n) + "] "
