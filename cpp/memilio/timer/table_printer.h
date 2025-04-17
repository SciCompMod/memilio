#ifndef MIO_TIMER_TABLE_PRINTER_H
#define MIO_TIMER_TABLE_PRINTER_H

#include "memilio/timer/definitions.h"
#include "memilio/timer/registration.h"
#include "spdlog/fmt/bundled/base.h"
#include "spdlog/fmt/bundled/format.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <iostream>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace mio
{
namespace timing
{

namespace details
{

template <class T>
class Table
{
    struct Row {
        std::string name;
        std::vector<T> values;
    };

public:
    Table(const std::string& name, const std::vector<std::string>& column_names,
          const std::vector<std::string>& row_names = {})
        : m_name(name)
        , m_column_names(column_names)
        , m_rows()
    {
        m_rows.reserve(row_names.size());
        for (auto& rn : row_names) {
            add_row(rn);
        }
    }

    void add_row(const std::string& name, const std::vector<T>& values)
    {
        assert(values.size() == m_column_names.size());
        m_rows.push_back({name, values});
    }

    void add_row(const std::string& name)
    {
        add_row(name, std::vector<T>(m_column_names.size()));
    }

    T& operator()(size_t row, size_t col)
    {
        assert(row < m_rows.size());
        assert(col < m_rows[row].values.size());
        return m_rows[row].values[col];
    }

    std::string& get_row_name(size_t row)
    {
        assert(row < m_rows.size());
        return m_rows[row].name;
    }

    std::string& get_col_name(size_t col)
    {
        assert(col < column_names.size());
        return m_column_names[col];
    }

    std::string& get_name()
    {
        return m_name;
    }

    size_t rows() const
    {
        return m_rows.size();
    }

    size_t cols() const
    {
        return m_column_names.size();
    }

private:
    std::string m_name;
    std::vector<std::string> m_column_names;
    std::vector<Row> m_rows;
};

} // namespace details

class TablePrinter : public Printer
{
public:
    void set_time_format(std::string format_string)
    {
        m_time_format = format_string;
    }

    void print(const std::list<TimerRegistration>& timer_register, std::ostream& out = std::cout) override
    {
        auto [table, is_multithreaded, name_width, max_val] = create_table(timer_register);
        // calculate column_widths
        const size_t time_width   = fmt::format(m_time_format, max_val).size();
        const size_t thread_width = table.get_col_name(table.cols() - 1).size();
        std::vector<size_t> col_widths(table.cols() + 1);
        // note that col_width is offset by 1 relative to table, as uses index 0 for the first value,
        // but here we want to print the timer (or row) name first
        col_widths[0] = name_width;
        for (size_t col = 0; col < table.cols(); col++) {
            col_widths[col + 1] = std::max(time_width, table.get_col_name(col).size());
        }
        if (is_multithreaded) { // reduce size for "#Threads" column
            col_widths[table.cols()] = thread_width;
        }
        // vertical decorations
        const std::string separator = " | ";
        const std::string border_l  = "| ";
        const std::string border_r  = " |\n";
        // horizontal decoration, as a lambda
        const auto draw_hline = [&]() {
            for (auto& w : col_widths) {
                out << "+" + std::string(w + 2, '-');
            }
            out << "+\n";
        };

        draw_hline();
        // table header, prints table and column names
        out << border_l << std::setw(col_widths[0]) << std::left << table.get_name();
        for (size_t col = 0; col < table.cols(); col++) {
            out << separator << std::setw(col_widths[col + 1]) << std::left << table.get_col_name(col);
        }
        out << border_r;
        draw_hline();
        // print table content, only adding statistics when mulithreaded
        for (size_t row = 0; row < table.rows(); row++) {
            out << border_l << std::setw(col_widths[0]) << std::left << table.get_row_name(row);
            out << separator << std::setw(col_widths[1]) << fmt::format(m_time_format, table(row, 0));
            if (is_multithreaded) {
                const int& num_threads = static_cast<int>(table(row, 4));
                for (int col = 1; col < 4; col++) {
                    if (num_threads == 1) {
                        out << separator << std::setw(col_widths[col + 1]) << "-- ";
                    }
                    else {
                        out << separator << std::setw(col_widths[col + 1])
                            << fmt::format(m_time_format, table(row, col));
                    }
                }
                out << separator << std::setw(thread_width) << std::right << num_threads;
            }
            out << border_r;
        }
        draw_hline();
    }

private:
    inline static std::tuple<details::Table<double>, bool, size_t, double>
    create_table(const std::list<TimerRegistration>& timer_register)
    {
        std::vector<std::string> rows; // list of all timer names
        std::map<std::string, size_t> row_to_index; // map from name to index, used to fill table
        bool is_multithreaded = false; // keep track of whether a thread > 0 exists
        // map rows from thread 0 first, so the order of timers (somewhat) corresponds to their call order
        for (const auto& [name, _, thread] : timer_register) {
            if (thread == 0) {
                if (row_to_index.emplace(name, rows.size()).second) {
                    rows.push_back(name);
                }
            }
            else {
                is_multithreaded = true;
            }
        }
        // make a second pass to add timers from other threads
        // this does nothing, if all timers are used on thread 0 at least once
        for (auto& [name, _, thread] : timer_register) {
            if (thread != 0) {
                if (row_to_index.emplace(name, rows.size()).second) {
                    rows.push_back(name);
                }
            }
        }
        // create a table, omitting statistics if not multithreaded
        details::Table<double> table =
            is_multithreaded
                ? details::Table<double>("Timers", {"Elapsed Time", "Min", "Max", "Average", "#Threads"}, rows)
                : details::Table<double>("Timers", {"Elapsed Time"}, rows);
        size_t name_width = 6; // width of the name column, defaults to length of "Timers".
        // shorthands for the columns, so there are fewer mystery numbers
        const int elapsed = 0, min = 1, max = 2, avg = 3, num = 4;
        // accumulate elapsed time and gather statistics in the table
        // averages are calculated later, using finished values from elapsed and num
        for (auto& [name, timer, thread] : timer_register) {
            const auto row  = row_to_index[name];
            const auto time = time_in_seconds(timer.get_elapsed_time());
            table(row, elapsed) += time;
            if (is_multithreaded) {
                if (table(row, min) == 0.0) {
                    // we need to overwrite the first min value, otherwise it would stay at 0
                    table(row, min) = time;
                }
                else {
                    table(row, min) = std::min(table(row, min), time);
                }
                table(row, max) = std::max(table(row, max), time);
                table(row, num) += 1;
            }
        }
        double max_val = 0.0; // maximum table value to calculate formatted time width.
        for (size_t row = 0; row < table.rows(); row++) {
            max_val    = std::max(max_val, table(row, elapsed));
            name_width = std::max(name_width, table.get_row_name(row).size());
            if (is_multithreaded) {
                table(row, avg) = table(row, elapsed) / table(row, num);
            }
        }
        return {table, is_multithreaded, name_width, max_val};
    }

    std::string m_time_format = "{:e}";
};

} // namespace timing
} // namespace mio

#endif // MIO_TIMER_TABLE_PRINTER_H
