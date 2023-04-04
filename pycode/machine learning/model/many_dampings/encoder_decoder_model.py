import torch
import torch.nn as nn

from torch.utils.data import Dataset, DataLoader


class StoreItemDataset(Dataset):
    def __init__(self, cat_columns=[],
                 num_columns=[],
                 embed_vector_size=None, decoder_input=True,
                 ohe_cat_columns=False):
        super().__init__()
        self.sequence_data = None
        self.cat_columns = cat_columns
        self.num_columns = num_columns
        self.cat_classes = {}
        self.cat_embed_shape = []
        self.cat_embed_vector_size = embed_vector_size if embed_vector_size is not None else {}
        self.pass_decoder_input = decoder_input
        self.ohe_cat_columns = ohe_cat_columns
        self.cat_columns_to_decoder = False

    def get_embedding_shape(self):
        return self.cat_embed_shape

    def load_sequence_data(self, processed_data):
        self.sequence_data = processed_data

    def process_cat_columns(self, column_map=None):
        column_map = column_map if column_map is not None else {}
        for col in self.cat_columns:
            self.sequence_data[col] = self.sequence_data[col].astype(
                'category')
            if col in column_map:
                self.sequence_data[col] = self.sequence_data[col].cat.set_categories(
                    column_map[col]).fillna('#NA#')
            else:
                self.sequence_data[col].cat.add_categories(
                    '#NA#', inplace=True)
            self.cat_embed_shape.append(
                (len(self.sequence_data[col].cat.categories),
                 self.cat_embed_vector_size.get(col, 50)))

    def __len__(self):
        return len(self.sequence_data)

    def __getitem__(self, idx):
        row = self.sequence_data.iloc[[idx]]
        x_inputs = [torch.tensor(
            row['x_sequence'].values[0],
            dtype=torch.float32)]
        y = torch.tensor(row['y_sequence'].values[0], dtype=torch.float32)
        if self.pass_decoder_input:
            decoder_input = torch.tensor(
                row['y_sequence'].values[0][:, 1:],
                dtype=torch.float32)
        if len(self.num_columns) > 0:
            for col in self.num_columns:
                num_tensor = torch.tensor(
                    [row[col].values[0]],
                    dtype=torch.float32)
                x_inputs[0] = torch.cat((x_inputs[0], num_tensor.repeat(
                    x_inputs[0].size(0)).unsqueeze(1)), axis=1)
                decoder_input = torch.cat(
                    (decoder_input, num_tensor.repeat(
                        decoder_input.size(0)).unsqueeze(1)),
                    axis=1)
        if len(self.cat_columns) > 0:
            if self.ohe_cat_columns:
                for ci, (num_classes, _) in enumerate(self.cat_embed_shape):
                    col_tensor = torch.zeros(num_classes, dtype=torch.float32)
                    col_tensor[row[self.cat_columns[ci]
                                   ].cat.codes.values[0]] = 1.0
                    col_tensor_x = col_tensor.repeat(x_inputs[0].size(0), 1)
                    x_inputs[0] = torch.cat(
                        (x_inputs[0], col_tensor_x), axis=1)
                    if self.pass_decoder_input and self.cat_columns_to_decoder:
                        col_tensor_y = col_tensor.repeat(
                            decoder_input.size(0), 1)
                        decoder_input = torch.cat(
                            (decoder_input, col_tensor_y), axis=1)
            else:
                cat_tensor = torch.tensor(
                    [row[col].cat.codes.values[0] for col in self.cat_columns],
                    dtype=torch.long
                )
                x_inputs.append(cat_tensor)
        if self.pass_decoder_input:
            x_inputs.append(decoder_input)
            y = torch.tensor(
                row['y_sequence'].values[0][:, 0],
                dtype=torch.float32)
        if len(x_inputs) > 1:
            return tuple(x_inputs), y
        return x_inputs[0], y


class RNNEncoder(nn.Module):
    def __init__(
            self, rnn_num_layers=1, input_feature_len=1, sequence_len=168,
            hidden_size=100, bidirectional=False, device='cpu',
            rnn_dropout=0.2):
        super().__init__()
        self.sequence_len = sequence_len
        self.hidden_size = hidden_size
        self.input_feature_len = input_feature_len
        self.num_layers = rnn_num_layers
        self.rnn_directions = 2 if bidirectional else 1
        self.gru = nn.GRU(
            num_layers=rnn_num_layers,
            input_size=input_feature_len,
            hidden_size=hidden_size,
            batch_first=True,
            bidirectional=bidirectional,
            dropout=rnn_dropout
        )
        self.device = device

    def forward(self, input_seq):
        ht = torch.zeros(
            self.num_layers * self.rnn_directions, input_seq.size(0),
            self.hidden_size, device=self.device)
        if input_seq.ndim < 3:
            input_seq.unsqueeze_(2)
        gru_out, hidden = self.gru(input_seq, ht)
        print(gru_out.shape)
        print(hidden.shape)
        if self.rnn_directions * self.num_layers > 1:
            num_layers = self.rnn_directions * self.num_layers
            if self.rnn_directions > 1:
                gru_out = gru_out.view(
                    input_seq.size(0),
                    self.sequence_len, self.rnn_directions, self.hidden_size)
                gru_out = torch.sum(gru_out, axis=2)
            hidden = hidden.view(
                self.num_layers, self.rnn_directions, input_seq.size(0),
                self.hidden_size)
            if self.num_layers > 0:
                hidden = hidden[-1]
            else:
                hidden = hidden.squeeze(0)
            hidden = hidden.sum(axis=0)
        else:
            hidden.squeeze_(0)
        return gru_out, hidden


class DecoderCell(nn.Module):
    def __init__(self, input_feature_len, hidden_size, dropout=0.2):
        super().__init__()
        self.decoder_rnn_cell = nn.GRUCell(
            input_size=input_feature_len,
            hidden_size=hidden_size,
        )
        self.out = nn.Linear(hidden_size, 1)
        self.attention = False
        self.dropout = nn.Dropout(dropout)

    def forward(self, prev_hidden, y):
        rnn_hidden = self.decoder_rnn_cell(y, prev_hidden)
        output = self.out(rnn_hidden)
        return output, self.dropout(rnn_hidden)


class EncoderDecoderWrapper(nn.Module):
    def __init__(self, encoder, decoder_cell, output_size=3,
                 teacher_forcing=0.3, sequence_len=336, decoder_input=True,
                 device='cpu'):
        super().__init__()
        self.encoder = encoder
        self.decoder_cell = decoder_cell
        self.output_size = output_size
        self.teacher_forcing = teacher_forcing
        self.sequence_length = sequence_len
        self.decoder_input = decoder_input
        self.device = device

    def forward(self, xb, yb=None):
        if self.decoder_input:
            decoder_input = xb[-1]
            input_seq = xb[0]
            if len(xb) > 2:
                encoder_output, encoder_hidden = self.encoder(
                    input_seq, *xb[1: -1])
            else:
                encoder_output, encoder_hidden = self.encoder(input_seq)
        else:
            if type(xb) is list and len(xb) > 1:
                input_seq = xb[0]
                encoder_output, encoder_hidden = self.encoder(*xb)
            else:
                input_seq = xb
                encoder_output, encoder_hidden = self.encoder(input_seq)
        prev_hidden = encoder_hidden
        outputs = torch.zeros(
            input_seq.size(0),
            self.output_size, device=self.device)
        y_prev = input_seq[:, -1, 0].unsqueeze(1)
        for i in range(self.output_size):
            step_decoder_input = torch.cat(
                (y_prev, decoder_input[:, i]), axis=1)
            if (yb is not None) and (i > 0) and (torch.rand(1) < self.teacher_forcing):
                step_decoder_input = torch.cat(
                    (yb[:, i].unsqueeze(1), decoder_input[:, i]), axis=1)
            rnn_output, prev_hidden = self.decoder_cell(
                prev_hidden, step_decoder_input)
            y_prev = rnn_output
            outputs[:, i] = rnn_output.squeeze(1)
        return outputs
