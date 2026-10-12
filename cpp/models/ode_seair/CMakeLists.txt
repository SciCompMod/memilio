memilio_add_model(ode_seair
    infection_state.h
    model.h
    model.cpp
    parameters.h
)
target_link_libraries(ode_seair PUBLIC AD::AD)
