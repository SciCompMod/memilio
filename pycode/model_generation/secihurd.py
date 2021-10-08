import generator as gen

project_path = ""

if project_path == "":
    from subprocess import check_output
    project_path = check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + "/cpp/"

header_file = "epidemiology/secir/secihurd.h"
source_file = "examples/secihurd.cpp"

### compartments and names ###
s, e, c, i, h, u, r, d = gen.compartments("Susceptible Exposed Carrier Infected Hospitalized ICU Recovered Dead")
gen.name("secihurd_example")
gen.namespace("gen")

### parameters ###
d0          = gen.parameter("StartDay",                        "double", 0)
season      = gen.parameter("Seasonality",                     "double", 0)
tinc        = gen.parameter("IncubationTime",                  "double", 5.2)
tinfmild    = gen.parameter("InfectiousTimeMild",              "double", 6)
tserint     = gen.parameter("SerialInterval",                  "double", 4.2)
thosp2home  = gen.parameter("HospitalizedToHomeTime",          "double", 1)
thome2hosp  = gen.parameter("HomeToHospitalizedTime",          "double", 5)
thosp2icu   = gen.parameter("HospitalizedToICUTime",           "double", 2)
ticu2home   = gen.parameter("ICUToHomeTime",                   "double", 8)
ticu2death  = gen.parameter("ICUToDeathTime",                  "double", 5)
inf_prob    = gen.parameter("InfectionProbabilityFromContact", "double", 0.05)
carr_infec  = gen.parameter("RelativeCarrierInfectability",    "double", 1)
alpha       = gen.parameter("AsymptoticCasesPerInfectious",    "double", 0.09)
beta        = gen.parameter("RiskOfInfectionFromSympomatic",   "double", 0.25)
rho         = gen.parameter("HospitalizedCasesPerInfectious",  "double", 0.2)
theta       = gen.parameter("ICUCasesPerHospitalized",         "double", 0.25)
delta       = gen.parameter("DeathsPerHospitalized",           "double", 0.3)
tnt_cap     = gen.parameter("TestAndTraceCapacity",            "double", "double(std::numeric_limits<double>::max())")
icu_cap     = gen.parameter("ICUCapacity",                     "double", 0)
tinfasym    = gen.parameter("InfectiousTimeAsymptomatic",      "double", 1/ (0.5 / (5.2 - 4.2)) + 0.5 * 6)

### declare some variables ###
icu_op, tnt_req, risk_symp = gen.variables("icu_occupancy test_and_trace_required risk_from_symptomatic")
n, divN, season_val, cont_freq_eff = gen.variables("N divN season_val cont_freq_eff")
prob_hosp2icu, prob_hosp2dead= gen.variables("prob_hosp2icu prob_hosp2dead")
dummy_R2, dummy_R3, dummy_S = gen.variables("dummy_R2 dummy_R3 dummy_S")

t = gen.time()

from sympy import Function # TDOD: make model version of this
cos = Function(r"epi::smoother_cosine")
sin = Function("sin")
mod = Function(r"std::fmod")

eq = gen.equation

eq(icu_op, "=", 0)
eq(tnt_req, "=", 0)
eq(dummy_R3, "=", 0.5 / (tinc - tserint))
eq(tnt_req, "=", (1 - alpha) * dummy_R3 * c)
eq(icu_op, "=", u)
eq(s, "=", 0)
eq(e, "=", 0)
eq(dummy_R2, "=", 1 / (2 * tserint - tinc))
eq(dummy_R3, "=", 0.5 / (tinc - tserint))
eq(risk_symp, "=", cos(tnt_req, tnt_cap, tnt_cap * 5, beta, 0))
eq(season_val, "=", (1 + season * sin(3.141592653589793 * (mod((d0 + t), 365) / 182.5 + 0.5))))
eq(cont_freq_eff, "=", season_val * 1)
eq(n, "=", s+e+c+i+h+u+r)
eq(divN, "=", 1/n)
eq(dummy_S, "=", s * cont_freq_eff * divN * inf_prob * (carr_infec * c + risk_symp * i))
eq(s, "-=", dummy_S)
eq(e, "+=", dummy_S)
eq(prob_hosp2icu, "=", cos(icu_op, 0.9 * icu_cap, icu_cap, theta, 0))
eq(prob_hosp2dead, "=", theta - prob_hosp2icu)
eq(e, "-=", dummy_R2 * e)
eq(c, "=", dummy_R2 * e - ((1 - alpha) * dummy_R3 + alpha / tinfasym) * c)
eq(i, "=", (1 - alpha) * dummy_R3 * c - ((1 - rho) / tinfmild + rho / thome2hosp) * i)
eq(h, "=", rho / thome2hosp * i - ((1 - theta) / thosp2home + theta / thosp2icu) * h)
eq(u, "=", -((1 - delta) / ticu2home + delta / ticu2death) * u)
eq(u, "+=", prob_hosp2icu / thosp2icu * h)
eq(r, "=", alpha / tinfasym * c + (1-rho) / tinfmild * i + (1-theta) / thosp2home * h + (1 - delta) / ticu2home * u)
eq(d, "=", delta / ticu2death * u)
eq(d, "+=", prob_hosp2dead / thosp2icu * h)

gen.print_header(project_path + header_file)
gen.print_source(header_file, 0, 50, 0.1, [99760, 100, 50, 50, 20, 10, 10, 0], project_path + source_file)