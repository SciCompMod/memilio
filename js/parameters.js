function Parameters($form) {
    let listeners = [];


    $form
        .submit(function(event){
            event.preventDefault();
            event.stopPropagation();
            
            let serialized = $form.serializeArray();

            listeners
                .forEach(l => {
                    l(serialized);
                });
        });

    return {
        onsubmit: function(listener) {
            listeners.push(listener);
        }
    }
}