from dataloader import DynamicsDataset
from torch.utils.data import DataLoader
import torch
import torch.nn as nn
from networks import SimpleDynamicsNetwork
from tqdm import tqdm

if __name__ == '__main__':
    BATCH_SIZE = int(4e5)
    EPOCHS = 30

    TRAIN_RATIO = 0.95

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print('Device being used: ', device)

    network = SimpleDynamicsNetwork().to(device)

    #workin on dataset
    dataset = DynamicsDataset('./dataset/')

    TRAIN_SIZE = int(TRAIN_RATIO*len(dataset))
    TEST_SIZE = len(dataset) - TRAIN_SIZE
    print(f'Training Size: ', TRAIN_SIZE)
    print(f'Testing Size: ', TEST_SIZE)

    train_dataset, test_dataset = torch.utils.data.random_split(dataset, [TRAIN_SIZE, TEST_SIZE])

    dataloader = DataLoader(train_dataset,shuffle=True,num_workers=4,batch_size=BATCH_SIZE,pin_memory=False)

    adam = torch.optim.Adam(network.parameters(),lr=1e-2)
    scheduler = torch.optim.lr_scheduler.MultiStepLR(adam, [20]*1, gamma=0.1)
    loss = nn.MSELoss(reduction='mean')


    x_test, labels_test = test_dataset[:]
    x_test = x_test.to(device)
    labels_test = labels_test.to(device)
    for epoch in range(EPOCHS):
        print('EPOCH: ',epoch)
        train_error = 0
        for idx,batch in enumerate(tqdm(dataloader)):
            adam.zero_grad()

            x, labels = batch
            x = x.to(device)
            labels = labels.to(device)

            out = network(x)

            mse = loss(labels,out)
            train_error += mse

            mse.backward()
            adam.step()
        scheduler.step()
        print(f'Training Error: {train_error}')

        out = network(x)
        mse = loss(labels,out)
        print(f'Testing Error: {mse}')
    torch.save(network, './models/hcw_simple_model.pt')
