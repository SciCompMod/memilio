function Parameters($form) {
    let listeners = [];

    function buildParameters() {
        let serialized = $form.serializeArray();
        let parameters = {};
        serialized.forEach(s => {
            parameters[s.name] = s.value;
        });
        return parameters;
    }

    $form
        .on('change', 'input', function(event){
            event.preventDefault();
            event.stopPropagation();
            
            let parameters = buildParameters();

            listeners
                .forEach(l => {
                    l(serialized);
                });
        });

    return {
        onchange: function(listener) {
            listeners.push(listener);
        },
        getParameters: function() {
            return buildParameters();
        }
    }
}