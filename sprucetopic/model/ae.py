
import torch; torch.manual_seed(0)
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader
from torch.utils.data import Dataset

class Stacklayers(nn.Module):
	def __init__(self,input_size,layers):
		super(Stacklayers, self).__init__()
		self.layers = nn.ModuleList()
		self.input_size = input_size
		for next_l in layers:
			self.layers.append(nn.Linear(self.input_size,next_l))
			self.layers.append(self.get_activation())
			nn.BatchNorm1d(next_l)
			self.input_size = next_l
		
	def forward(self, input_data):
		for layer in self.layers:
			input_data = layer(input_data)
		return input_data

	def get_activation(self):
		return nn.ReLU()


class Decoder(nn.Module):
	def __init__(self, input_dims, latent_dims):
		super(Decoder, self).__init__()
		self.linear1 = nn.Linear(latent_dims, 128)
		self.linear2 = nn.Linear(128, input_dims)
	def forward(self, z):
		output = F.relu(self.linear1(z))
		x_hat = F.relu(self.linear2(output))
		return x_hat


class Autoencoder(nn.Module):
	def __init__(self,input_dims,out_dims,latent_dims,layers):
		super(Autoencoder, self).__init__()
		self.encoder = Stacklayers(input_dims,layers)
		self.decoder = Decoder(out_dims,latent_dims)
	def forward(self, x):
		z = self.encoder(x)
		return self.decoder(z)


class TabularDataset(Dataset):
	def __init__(self, x,y):
		self.X = x
		self.y = y	
	def __len__(self):
		return len(self.X)
	def __getitem__(self, idx):
		return [self.X[idx], self.y[idx]]

def load_data(x,y,device,batch_size):
	x = x.astype('float32')
	x = torch.from_numpy(x).to(device)
	return DataLoader(TabularDataset(x,y), batch_size=batch_size, shuffle=True)

def train(autoencoder, data,device, epochs,lr):
	opt = torch.optim.Adam(autoencoder.parameters(),lr=lr)
	mse = nn.MSELoss()
	loss_values =[]
	for epoch in range(epochs):
		loss = 0
		for indx,(x,y) in enumerate(data):
			x = x.to(device)
			opt.zero_grad()
			x_hat = autoencoder(x)
			train_loss = mse(x_hat,x)
			train_loss.backward()
			opt.step()
			loss += train_loss.item()
		if epoch % 30 == 0:  
			print('====> Epoch: {} Average loss: {:.4f}'.format(epoch, loss/len(data)))
		loss_values.append(loss/len(data))
	return loss_values
