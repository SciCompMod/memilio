import pickle


class Run_Saver:

    __scores = None
    __history = None
    __predtest = None
    __testY = None
    __y_scaler = None
    __x_mlp_scaler = None
    __calculation_time = None
    __n_epochs = None
    __n_batch_size = None
    __name = None
    __eta = None
    __loss_func = None
    __metrics_funcs = None
    __running_status = -1
    __dataset_name = None
    __ml_model_name = None


    def __init__(self, n_epochs, n_batch_size, name, eta, loss_func, metrics_funcs, dataset_name, ml_model_name, running_status=1):
        self.__n_epochs = n_epochs
        self.__n_batch_size = n_batch_size
        self.__name = name
        self.__eta = eta
        self.__loss_func = loss_func
        self.__metrics_funcs = metrics_funcs
        self.__dataset_name = dataset_name
        self.__ml_model_name = ml_model_name
        self.__running_status = running_status

    def setResults(self, scores, history, predtest, testY, y_scaler, x_mlp_scaler, calculation_time):
        self.__scores = scores
        self.__history = history
        self.__predtest = predtest
        self.__testY = testY
        self.__y_scaler = y_scaler
        self.__x_mlp_scaler = x_mlp_scaler
        self.__calculation_time = calculation_time
        self.__running_status = 0

    def getScores(self):
        return self.__scores

    def getHistory(self):
        return self.__history

    def getPredtest(self):
        return self.__predtest

    def getTestY(self):
        return self.__testY

    def getYScaler(self):
        return self.__y_scaler

    def getXMLPScaler(self):
        return self.__x_mlp_scaler

    def getCalculation_time(self):
        return self.__calculation_time

    def getN_epochs(self):
        return self.__n_epochs

    def getN_batch_size(self):
        return self.__n_batch_size

    def getName(self):
        return self.__name

    def getEta(self):
        return self.__eta

    def getLoss_func(self):
        return self.__loss_func

    def getMetrics_funcs(self):
        return self.__metrics_funcs

    def getRunning_status(self):
        return self.__running_status

    def getDataset_name(self):
        return self.__dataset_name

    def getMl_model_name(self):
        return self.__ml_model_name


    def save(self, file_name):
        if(len(file_name) > 4 and file_name[len(file_name)-4:] != ".run"):
            file_name += ".run"

        with open(file_name, "wb") as f:
            pickle.dump(self,f,pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def load(file_name):
        with open(file_name, "rb") as f:
            return pickle.load(f)
