/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: René Schmieding
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
#include "matchers.h"

#ifdef MEMILIO_HAS_JSONCPP

void Json::PrintTo(const Value& json, std::ostream* os)
{
    if (json.isObject()) {
        // move opening bracket into its own line
        *os << "\n";
    }
    const static auto js_writer = [] {
        StreamWriterBuilder swb;
        swb["indentation"] = "  ";
        return std::unique_ptr<StreamWriter>(swb.newStreamWriter());
    }();
    js_writer->write(json, os);
}

std::string json_type_to_string(Json::ValueType t)
{
    using namespace Json;
    switch (t) {
    case nullValue:
        return "Null";
    case intValue:
        return "Int";
    case uintValue:
        return "UInt";
    case realValue:
        return "Real";
    case stringValue:
        return "String";
    case booleanValue:
        return "Boolean";
    case arrayValue:
        return "Array";
    case objectValue:
        return "Object";
    default:
        assert(false && "Unreachable");
        return "";
    }
}

#endif // MEMILIO_HAS_JSONCPP
