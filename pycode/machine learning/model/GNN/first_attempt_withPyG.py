from torch_geometric.nn import GCNConv
import torch
import json
import numpy as np
from torch_geometric.data import Data
import os
import pickle
from torch_geometric.loader import DataLoader
import torch.nn.functional as F
import argparse
from torch_geometric.nn import global_mean_pool


path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(
    os.path.dirname(
        os.path.realpath(os.path.dirname(os.path.realpath(path)))),
    'data_GNN_3pop_nodamp_30days_1k')

file = open(os.path.join(path_data, 'data_secir_age_groups.pickle'), 'rb')
data_secir = pickle.load(file)

new_inputs = np.asarray(data_secir['inputs']).reshape(1000, 3, 240)
new_labels = np.asarray(data_secir['labels']).reshape(1000, 3, 1200)

edge_index = torch.tensor([[0, 1],
                           [1, 0],
                           [0, 2],
                           [2, 0],
                           [1, 2],
                           [2, 1]], dtype=torch.long)


dataset_x = torch.tensor(new_inputs, dtype=torch.float)
dataset_y = torch.tensor(new_labels, dtype=torch.float)
edge_attributes = [693, 345, 62, 22, 1103, 22]
data = Data(
    x=dataset_x[0], y=dataset_y[0], edge_index=edge_index.t().contiguous(),
    edge_attr=edge_attributes)

train_dataset_x = dataset_x[:800]
train_loader_x = DataLoader(train_dataset_x, batch_size=32, shuffle=True)

test_dataset_x = dataset_x[800:]
test_loader_x = DataLoader(test_dataset_x,  batch_size=32, shuffle=True)

train_dataset_y = dataset_y[:800]
train_loader_y = DataLoader(train_dataset_y, batch_size=32, shuffle=True)


test_dataset_y = dataset_y[800:]
test_loader_x = DataLoader(test_dataset_y)

# loader = DataLoader(dataset_x, batch_size=32, shuffle=True)

output_shape = data.y.shape[1]


class GCN(torch.nn.Module):
    def __init__(self, hidden_channels):
        super().__init__()
        torch.manual_seed(1234567)
        self.conv1 = GCNConv(data.num_features, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, output_shape)
        self.LinearLayer = torch.nn.Linear(output_shape, output_shape)

    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index)
        x = x.relu()
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.conv2(x, edge_index)

        x = self.LinearLayer(x, edge_index)
        return x


# class GCN(torch.nn.Module):
#     def __init__(self, hidden_channels):
#         super().__init__()
#         torch.manual_seed(1234567)
#         self.conv1 = GCNConv(data.num_features, hidden_channels)
#         self.conv2 = GCNConv(hidden_channels, output_shape)
#         self.LinearLayer = torch.nn.Linear(output_shape, output_shape)

#     def forward(self, x, edge_index, batch):
#         x = self.conv1(x, edge_index)
#         x = x.relu()
#         x = F.dropout(x, p=0.5, training=self.training)
#         x = self.conv2(x, edge_index)

#         # 2. Readout layer
#         x = global_mean_pool(x, batch)  # [batch_size, hidden_channels]


#         x = self.LinearLayer(x)
#         return x


model = GCN(hidden_channels=16)
print(model)


optimizer = torch.optim.Adam(model.parameters(), lr=0.005, weight_decay=5e-4)


def train():
    model.train()

    # Iterate in batches over the training dataset.
    for data, labels in zip(train_loader_x, train_loader_y):
        data_ = Data(x=data, y=labels, edge_index=edge_index.t().contiguous(),
                     edge_attr=edge_attributes)
        # Perform a single forward pass.
        out = model(data_.x, data_.edge_index, data_.batch)
        loss = F.mse_loss(out, data_.y)  # Compute the loss.
        loss.backward()  # Derive gradients.
        optimizer.step()  # Update parameters based on gradients.
        optimizer.zero_grad()  # Clear gradients.


def test(loader_x, loader_y):
    model.eval()

    mape = []
    # Iterate in batches over the training/test dataset.
    for data, labels in zip(loader_x, loader_y):
        data_ = Data(x=data, y=labels, edge_index=edge_index.t().contiguous(),
                     edge_attr=edge_attributes)
        out = model(data_, data_.edge_index, data_.batch)
        pred = out.argmax(dim=1)  # Use the class with highest probability.
        # Check against ground-truth labels.)
        mape_run = abs((pred-data_.y)/data_.y)
        mape.append(mape_run)
    return (1/len(test_loader_x))*sum(mape)


for epoch in range(1, 171):
    train()
    train_acc = test(train_loader_x, train_loader_y)

    print(f'Epoch: {epoch:03d}, Train Acc: {train_acc:.4f}')


###########
def train():
    model.train()
    optimizer.zero_grad()
    output = model(data.x, data.edge_index)
    loss_train = F.mse_loss(output, data.y)
    loss_train.backward(retain_graph=True)
    optimizer.step()
    return output, loss_train


# tensor with node features in size [number of nodes, number of node features]
# e.g. x = torch.tensor([[population node 1 ], [population node 2], [population node 3]], dtype=torch.float)
# x = [data[]]

# data = Data(x=x, edge_index=edge_index.t().contiguous())


###############  without  batches ####################
dataset_x = torch.tensor(new_inputs, dtype=torch.float)
dataset_y = torch.tensor(new_labels, dtype=torch.float)
edge_attributes = [693, 345, 62, 22, 1103, 22]
data = Data(
    x=dataset_x[0], y=dataset_y[0], edge_index=edge_index.t().contiguous(),
    edge_attr=edge_attributes)

train_dataset_x = dataset_x[:800]
train_loader_x = DataLoader(train_dataset_x)

test_dataset_x = dataset_x[800:]
test_loader_x = DataLoader(test_dataset_x)

train_dataset_y = dataset_y[:800]
train_loader_y = DataLoader(train_dataset_y)


test_dataset_y = dataset_y[800:]
test_loader_x = DataLoader(test_dataset_y)

# loader = DataLoader(dataset_x, batch_size=32, shuffle=True)

output_shape = data.y.shape[1]


class GCN(torch.nn.Module):
    def __init__(self, hidden_channels):
        super().__init__()
        torch.manual_seed(1234567)
        self.conv1 = GCNConv(data.num_features, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, output_shape)
        self.LinearLayer = torch.nn.Linear(output_shape, output_shape)

    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index)
        x = x.relu()
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.conv2(x, edge_index)

        x = self.LinearLayer(x)

        return x


model = GCN(hidden_channels=16)
print(model)


optimizer = torch.optim.Adam(model.parameters(), lr=0.005, weight_decay=5e-4)


def train():
    model.train()

    # Iterate in batches over the training dataset.
    for data, labels in zip(train_loader_x, train_loader_y):
        data_ = Data(x=data, y=labels, edge_index=edge_index.t().contiguous(),
                     edge_attr=edge_attributes)
        # Perform a single forward pass.
        out = model(data_.x.reshape(3, 240), data_.edge_index)
        loss = F.mse_loss(out, data_.y)  # Compute the loss.
        loss.backward()  # Derive gradients.
        optimizer.step()  # Update parameters based on gradients.
        optimizer.zero_grad()  # Clear gradients.


def test(loader_x, loader_y):
    model.eval()

    mape = []
    # Iterate in batches over the training/test dataset.
    for data, labels in zip(loader_x, loader_y):
        data_ = Data(x=data, y=labels, edge_index=edge_index.t().contiguous(),
                     edge_attr=edge_attributes)
        out = model(data_, data_.edge_index)
        # pred = out.argmax(dim=1)  # Use the class with highest probability.
        pred = out
        # Check against ground-truth labels.)
        mape_run = abs((pred-labels)/data.y)
        mape.append(mape_run)
    return (1/len(test_loader_x))*sum(mape)
