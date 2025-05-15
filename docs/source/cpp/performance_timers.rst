Performance Timers
==================

Here we present MEmilio's own timing framework. 


Timer usage
-----------

In this section we present how to use the AutoTimer class. This is the preferred way of using the timing framework, as
the class takes care of running the timer, managing its lifetime for later evaluation, and ensuring thread safety with
OpenMP.

An AutoTimer starts when it is created, and stops when it is destroyed - which usually happens at the next closing
bracket :code:`}` or the next :code:`return`. This design, automating the starting and stopping of a timer, is
intentionally limiting, because it helps to avoid several issues or mistakes that can arise when manually running
timers. 

You can also try out the `code example <https://github.com/SciCompMod/memilio/tree/main/cpp/examples>`__.


Timing in executables
~~~~~~~~~~~~~~~~~~~~~

Let's say you have set up a Simulation, and want to measure how long it takes, without the setup. Then you can write 

.. code-block:: cpp

    #include "memilio/timer/auto_timer.h"

    int main() {
        Simulation sim; ... // setup

        {
            mio::timer::AutoTimer<"my simulation"> my_timer; // my_timer starts here
            sim.advance(t_max); // run the simulation
        } // my_timer stops here

        ... // evaluate results
    }

and will see a table printed at the end of your programm, that lists the time it took to :code:`advance` next to the
timer named "my simulation". That's it!

You can add more timers like this, but make sure you use unique names, otherwise the same timer will be reused, and the
measured times will be added together. The name of the timer object itself (here :code:`my_timer`) is not important, as
long as the compiler does not complain about it.


Timing in the library
~~~~~~~~~~~~~~~~~~~~~

Adding timers in the library is not much different to adding timers in main, but avoiding name collisions can be more
difficult. Hence, we use the optional second tempalte argument of AutoTimer to specify its scope, as shown in the
following examples.

To measure the time a class member function takes, add a timer like this:

.. code-block:: cpp

    #include "memilio/timer/auto_timer.h"

    namespace foo

    class Bar {
        void baz() {
            mio::timer::AutoTimer<"baz", "foo::Bar"> timer;
            
            ... // rest of the function
        }
    };
    
    } // namespace foo

Or, when timing a free function:

.. code-block:: cpp

    #include "memilio/timer/auto_timer.h"

    namespace foo {

    void bar {
        AutoTimer<"bar", "foo"> timer;
            
        ... // rest of the function
    }

    } // namespace foo

The first string given to AutoTimer is the timer's name, the second the scope it is in. They are used in combination
to identify the timer, similar to a map key, so they must be unique. This can be effectively guaranteed, if the name
matches the function and the scope contains all enclosing namespaces, like in the examples above.

If the containing function is used, a summary with timing results will be printed where both timers will show up as
:code:`foo::Bar::baz` and :code:`foo::bar`, respectively.


General Recommendations
-----------------------

- **Do not time every detail.**
  While accurate, the timers are not very precise, so if you want to know how much time one or a few instructions take,
  use a profiler like (g)perf or likwid. Also, adding too many timers will clutter the timing results.

- **Only time computationally intensive code.**
  Similar to the last point, avoid timing small functions like setters and getters, and reserve using timers for the
  main compute loops or functions. While the timers add only a little overhead, it will become measurable when used too
  often.

- **Time entire functions.**
  Adding scopes for timing parts of main is fine, but you should avoid segmenting functions, either with scopes for
  AutoTimer or with manually run timers. The reason for this is related less to timers and more to code design, because
  if you can segment the function into multiple distinct parts, it is probably doing too many things, and should be
  separated into smaller functions. Also, adding scope (and thus indents) for AutoTimer does make code slightly harder
  to read.


The Timing Framework
--------------------

The main goals of this timing framework are to provide time measuring capabilities with minimal runtime overhead and
without having to plan around them. This means that accessing, starting and stopping a timer should be as fast as
possible, while the interfaces of the classes or functions that are to be timed should not change. Additionally, the
timer should work in parallel environments.

The solution to this is AutoTimer, whose usage was already shown above. There are, of course, some drawbacks. For
example, NamedTimer (the class used by AutoTimer) cannot be instantiated dynamically, as their name (and scope) have to
be known at compile time. This also means that adding a lot of timers will impact the time it takes to compile the code,
though a couple hundred timers should only take around an additional second.

Classes and their responsibilities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this section, we describe the main components of the timing framework and how they interact. For more details on a
specific component, view its API documentation. 

- **BasicTimer**:
  The foundation of the timing framework. BasicTimer is a very simple class, that defines the methods start, stop,
  reset, and get_elapsed_time. These are used by all other classes in this framework. Uses a wall clock, so if compute
  resources are shared with other tasks, the timing results may be higher than expexted. In debug builds, it will log
  errors whenever a member function was used incorrectly, e.g. when start was called twice.

- **TimerRegistration**:
  This simple struct is used to keep track of timers and some additional information, but does not manage their storage.
  It consists of two strings for name and scope, a reference to a BasicTimer, and a thread id. The thread id specifies
  which thread the timer is used in, which could differ from the thread is is created by.

- **Printer**:
  A pure virtual class defining a print method to evaluate and output timing results via a list of TimerRegistrations.
  Implemented by TablePrinter and ListPrinter.

- **TimerRegistrar**:
  Keeps track of timers via a list of TimerRegistrations, and holds a Printer that can be used to display all
  registered timers after the end of main. Timers can be registered by passing a TimerRegistration to its add_timer
  method. Uses a singleton pattern to provide global access to the same object, that is, the only way to obtain a
  TimerRegistrar object is by using its get_instance method, which returns a reference to a static object. Importantly,
  this class does not manage or own timer objects, and there is intentionally no methods that retrive or delete
  TimerRegistrations.

- **NamedTimer**:
  Inherits from BasicTimer, with the main purpose of managing the lifetime, access, and registration of a timer.
  This is done using a singleton pattern, similar to TimerRegistrar, but the reference returned by get_instance is
  thread_local as well as static. The template parameters Name and Scope allow using more than one NamedTimer, since
  different template arguments define a different type. This effectively creates a global compile-time map, mapping a
  Name and Scope to a BasicTimer. Additionally, the NamedTimer registers itself automatically, and will only be
  destroyed after the TimerRegistrar.

- **AutoTimer**:
  Automates running an existing timer, by calling start in its constructor, and stop in its destructor. The timer used
  can be either specified via the Name and Scope template, fetching the corresponding NamedTimer internally, or by
  passing an lvalue reference to a BasicTimer.

Using NamedTimer and BasicTimer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Preferably, you should use AutoTimer where possible, as its limiting design helps to avoid common errors, for example
with parallel regions. But, if you have to, you can use a NamedTimer directly without any extra work:

.. code-block:: cpp

    #include "memilio/timer/named_timer.h"

    namespace foo

    class Bar {
        void baz() {
            mio::timer::NamedTimer<"baz", "foo::Bar">::get_instance().start();
            
            ... // rest of the function

            mio::timer::NamedTimer<"baz", "foo::Bar">::get_instance().stop();
        }
    };
    
    } // namespace foo

This will behave exactly like the AutoTimer in the example above, while also allowing you to use the reset or
get_elapsed_time methods defined by BasicTimer.

Last but not least, you can also use a BasicTimer directly. This means that you will have to manually take care of
the timer object, threading and evaluation. If you add such a BasicTimer to the TimerRegistrar, you will probably need
to disable the final timer summary, and call print manually. Of course, you can also make your own list of registrations
and use a Printer directly.
