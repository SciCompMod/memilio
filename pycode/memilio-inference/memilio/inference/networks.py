import tensorflow as tf
from tensorflow.keras.layers import LSTM, Bidirectional, Dense
from tensorflow.keras.models import Sequential

from bayesflow import default_settings as defaults

from bayesflow.helper_networks import MultiConv1D


class TwoLevelSequenceNetwork(tf.keras.Model):
    """Adjustment of SequensceNetwork from BayesFlow for TwoLevelModel
    """

    def __init__(
        self, summary_dim=10, num_conv_layers=2, lstm_units=128, bidirectional=False, conv_settings=None, **kwargs
    ):
        """Creates a stack of inception-like layers followed by an LSTM network, with the idea
        of learning vector representations from multivariate time series data.

        Parameters
        ----------
        summary_dim     : int, optional, default: 10
            The number of learned summary statistics.
        num_conv_layers : int, optional, default: 2
            The number of convolutional layers to use.
        lstm_units      : int, optional, default: 128
            The number of hidden LSTM units.
        conv_settings   : dict or None, optional, default: None
            The arguments passed to the `MultiConv1D` internal networks. If `None`,
            defaults will be used from `default_settings`. If a dictionary is provided,
            it should contain the following keys:
            - layer_args      (dict) : arguments for `tf.keras.layers.Conv1D` without kernel_size
            - min_kernel_size (int)  : the minimum kernel size (>= 1)
            - max_kernel_size (int)  : the maximum kernel size
        bidirectional   : bool, optional, default: False
            Indicates whether the involved LSTM network is bidirectional (forward and backward in time)
            or unidirectional (forward in time). Defaults to False, but may increase performance.
        **kwargs        : dict
            Optional keyword arguments passed to the __init__() method of tf.keras.Model
        """

        super().__init__(**kwargs)

        # Take care of None conv_settings
        if conv_settings is None:
            conv_settings = defaults.DEFAULT_SETTING_MULTI_CONV

        self.net = Sequential([MultiConv1D(conv_settings)
                              for _ in range(num_conv_layers)])

        self.lstm = Bidirectional(
            LSTM(lstm_units)) if bidirectional else LSTM(lstm_units)
        self.out_layer = Dense(summary_dim, activation="linear")
        self.summary_dim = summary_dim

    def call(self, x, **kwargs):
        """Performs a forward pass through the network by first passing `x` through the sequence of
        multi-convolutional layers and then applying the LSTM network.

        Parameters
        ----------
        x : tf.Tensor
            Input of shape (batch_size, n_groups, n_time_steps, n_time_series)

        Returns
        -------
        out : tf.Tensor
            Output of shape (batch_size, n_groups, summary_dim)
        """

        # Reshape: (batch_size, n_groups, timepoints, values) -> (batch_size * n_groups, timepoints, values)
        batch_size, n_groups, timepoints, values = x.shape
        out = tf.reshape(x, (-1, timepoints, values))

        out = self.net(out, **kwargs)
        out = self.lstm(out, **kwargs)
        out = self.out_layer(out, **kwargs)

        # Reshape: (batch_size * n_groups, summary_dim) -> (batch_size, n_groups, summary_dim)
        out = tf.reshape(out, (batch_size, n_groups, -1))

        return out
