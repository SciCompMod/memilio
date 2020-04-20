export { Parameters };

class Parameters {
  constructor($form) {
    this.listeners = [];
    this.$form = $form;
    let self = this;

    $form.on('change', 'input', function (event) {
      event.preventDefault();
      event.stopPropagation();

      let parameters = self.buildParameters();

      self.listeners.forEach((l) => {
        l(parameters);
      });
    });
  }

  buildParameters() {
    let serialized = this.$form.serializeArray();
    let parameters = {};
    serialized.forEach((s) => {
      parameters[s.name] = parseFloat(s.value);
    });
    return parameters;
  }

  onchange(listener) {
    this.listeners.push(listener);
  }
  getParameters() {
    return this.buildParameters();
  }
}
