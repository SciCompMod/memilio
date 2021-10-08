from collections import OrderedDict
import generator.Model

class ModelManager:
    models = OrderedDict()

    @classmethod
    def get_active(self):
        if not self.models: # if models is empty
            self.models[0] = generator.Model(0)
        return next(reversed(self.models.values()))

    @classmethod
    def get_all(self):
        return self.models.values() if self.models else None

    @classmethod
    def set_active(self, model):
        self.models[model.num] = model
        self.models.move_to_end(model.num)

    @classmethod
    def add_model(self):
        num = 0
        while self.models.get(num):
            num += 1
        model = generator.Model(num)
        self.models[num] = model
        return model

    @classmethod
    def remove_model(self, model):
        self.models.pop(model.num)