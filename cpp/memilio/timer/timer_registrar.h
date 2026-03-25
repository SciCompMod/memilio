/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Rene Schmieding
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef MIO_TIMER_TIMER_REGISTRAR_H
#define MIO_TIMER_TIMER_REGISTRAR_H

#include "memilio/io/binary_serializer.h"
#include "memilio/io/default_serialize.h"
#include "memilio/timer/basic_timer.h"
#include "memilio/timer/registration.h"
#include "memilio/timer/table_printer.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/miompi.h"

#include <cassert>
#include <iostream>
#include <list>
#include <memory>
#include <mutex>

namespace mio
{
namespace timing
{

/// @brief TimerRegistrar is a singleton to keep track of and print timers. Does not manage storage.
class TimerRegistrar
{
public:
    /**
     * @brief Access the TimerRegistrar.
     *
     * TimerRegistrar acts as a singleton, using this method to access the instance.
     * This makes TimerRegistrar accessible everywhere, without having to manage a TimerRegistrar object or pass an
     * instance to any class or function that uses timers. In particular, it allows e.g. NamedTimer to register itself
     * automatically, while guaranteeing a correct destruction order.
     *
     * @return A reference to the TimerRegistrar instance.
     */
    static TimerRegistrar& get_instance()
    {
        static TimerRegistrar t;
        return t;
    }

    /**
     * @brief Add a new timer to the TimerRegistrar.
     *
     * Note that the TimerRegistrar does not manage any storage, it only keeps a list of the timers and their names.
     * Therefore, if you manually call add_timer, you have to make sure that the registered timer lives long enough
     * for e.g. the last print_timer call. You may want to disable_final_timer_summary, and look at NamedTimer as an
     * example that does guarantee correct object lifetime.  
     *
     * @param registration A new registration with a name, a reference to a timer, and the thread id it is created on.
     */
    void add_timer(TimerRegistration&& registration)
    {
        m_registration_lock.lock();
        m_register.emplace_back(std::move(registration));
        m_registration_lock.unlock();
    }

    /// @brief Returns a read only list of all TimerRegistration%s. Can be used to print or evaluate timers.
    const auto& get_register() const
    {
        return m_register;
    }

    /**
     * @brief Print all timers using a Printer.
     * 
     * By default, uses TablePrinter to write to stdout. The printer can be changed using the set_printer member.
     *
     * @param out Any output stream, defaults to std::cout.
     */
    void print_timers(std::ostream& out = std::cout) const
    {
        PRAGMA_OMP(single nowait)
        {
            m_printer->print(m_register, out);
        }
    }

    /**
     * @brief Print all timers gathered on the given MPI Comm using a Printer.
     * 
     * By default, uses TablePrinter to write to stdout. The printer can be changed using the set_printer member.
     *
     * @param comm An MPI communicator. Make sure this method is called by all ranks on comm!
     * @param out Any output stream, defaults to std::cout.
     */
    void print_gathered_timers(mpi::Comm comm, std::ostream& out = std::cout)
    {
        PRAGMA_OMP(single nowait)
        {
            std::list<BasicTimer> external_timers; // temporary storage for gathered registrations
            std::list<TimerRegistration> gathered_register; // combined registrations of all ranks
            gather_timers(comm, external_timers, gathered_register);
            if (mpi::is_root()) {
                m_printer->print(gathered_register, out);
            }
        }
    }

    /// @brief Prevent the TimerRegistrar from calling print_timers on exit from main.
    void disable_final_timer_summary()
    {
        m_print_on_death = false;
    }

    /**
     * @brief Replace the internal printer used for print_timers.
     *
     * To allow making a final timer summary on exit from main, the TimerRegistrar has to take ownership of the printer.
     * Usage may look as follows:
     * ```
     * auto printer = std::make_unique<mio::timing::TablePrinter>();
     * mio::timing::TimerRegistrar::get_instance().set_printer(std::move(printer));
     * ```
     * Note that after passing the printer using std::move, it will contain a nullptr and can not be used again.
     *
     * @param printer A unique pointer to something that implements Printer.
     */
    void set_printer(std::unique_ptr<Printer>&& printer)
    {
        m_printer.swap(printer);
        // old value of m_printer (now stored in printer) is deleted at end of scope
    }

    /**
     * @brief TimerRegistrar must not be copied or moved, use `TimerRegistrar::get_instance()` to access it.
     * @{
     */
    TimerRegistrar(TimerRegistrar&)  = delete;
    TimerRegistrar(TimerRegistrar&&) = delete;
    void operator=(TimerRegistrar&)  = delete;
    void operator=(TimerRegistrar&&) = delete;
    /** @} */

private:
    /// @brief Instead of constructing a TimerRegistrar, use its static method `TimerRegistrar::get_instance()`.
    TimerRegistrar() = default;

    /// @brief Specify a destructor to allow printing all timers after exit from main.
    ~TimerRegistrar()
    {
        if (m_print_on_death) {
            PRAGMA_OMP(single nowait)
            {
                // print_timers already uses a "single nowait", but the cout for titling the table needs one, too
                std::cout << "Final Timer Summary\n";
                print_timers();
            }
        }
    }

    /// @brief Gather timers from all ranks, using external_timers as temporary timer storage for gathered_register.
    void gather_timers(mpi::Comm comm, std::list<BasicTimer>& external_timers,
                       std::list<TimerRegistration>& gathered_register) const
    {
#ifndef MEMILIO_ENABLE_MPI
        mio::unused(comm, external_timers);
#endif
        if (mpi::is_root()) {
            std::ranges::transform(m_register, std::back_inserter(gathered_register), [](const TimerRegistration& r) {
                return TimerRegistration{r.name, r.scope, r.timer, r.thread_id, 0};
            });
        }
        if (comm == nullptr) {
            if (mpi::is_root()) {
                log_error("Got nullptr as MPI Comm. Only timers on root rank are gathered.");
            }
            return;
        }
#ifdef MEMILIO_ENABLE_MPI
        // name, scope, timer, thread_id
        using GatherableRegistration = std::tuple<std::string, std::string, BasicTimer, int>;
        int comm_size;
        MPI_Comm_size(comm, &comm_size);

        if (mpi::is_root()) {
            for (int snd_rank = 1; snd_rank < comm_size; snd_rank++) { // skip root rank!
                int bytes_size;
                MPI_Recv(&bytes_size, 1, MPI_INT, snd_rank, 0, comm, MPI_STATUS_IGNORE);
                ByteStream bytes(bytes_size);
                MPI_Recv(bytes.data(), bytes.data_size(), MPI_BYTE, snd_rank, 0, mpi::get_world(), MPI_STATUS_IGNORE);

                auto rec_register = deserialize_binary(bytes, Tag<std::vector<GatherableRegistration>>{});
                if (!rec_register) {
                    log_error("Error receiving ensemble results from rank {}.", snd_rank);
                }
                std::ranges::transform(
                    rec_register.value(), std::back_inserter(gathered_register), [&](const GatherableRegistration& r) {
                        const auto& [name, scope, timer, thread_id] = r;
                        external_timers.push_back(timer);
                        return TimerRegistration{name, scope, external_timers.back(), thread_id, snd_rank};
                    });
            }
        }
        else {
            std::vector<GatherableRegistration> snd_register;
            std::ranges::transform(m_register, std::back_inserter(snd_register), [](const TimerRegistration& r) {
                return GatherableRegistration{r.name, r.scope, r.timer, r.thread_id};
            });
            ByteStream bytes = serialize_binary(snd_register);
            int bytes_size   = int(bytes.data_size());
            MPI_Send(&bytes_size, 1, MPI_INT, 0, 0, comm);
            MPI_Send(bytes.data(), bytes.data_size(), MPI_BYTE, 0, 0, comm);
        }
#endif
    }

    std::unique_ptr<Printer> m_printer = std::make_unique<TablePrinter>(); ///< A printer to visualize all timers.
    bool m_print_on_death              = true; ///< Whether to call m_printer during the destructor.
    std::list<TimerRegistration> m_register; ///< List that allows access to timers without having their name.
    std::mutex m_registration_lock; ///< Lock to safeguard m_register against concurrent writes.
};

} // namespace timing
} // namespace mio

#endif // MIO_TIMER_TIMER_REGISTRAR_H
